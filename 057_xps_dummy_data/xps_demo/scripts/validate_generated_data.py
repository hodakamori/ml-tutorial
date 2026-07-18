"""生成した XPS ダミーデータの自動検証スクリプト。

generate_xps_data.py で作られたデータが仕様どおりかを確認する。
プレビュー画像（preview/*.png）も併せて生成する。

チェック項目（プラン Phase 6 / 完成条件に対応）:
    * スペクトルが 50 個存在する
    * 各スペクトルが 201 点ある
    * 結合エネルギーが単調に減少する
    * 欠損値・無限値がない
    * 強度が 0 以上である
    * 真のピーク面積比の合計が 1 である
    * 表面処理の強さと酸素含有ピーク（C-O, O-C=O）の平均面積比に期待した傾向がある
    * 同一条件内にも十分なばらつきがある
    * 同じ乱数シードで同じデータを再生成できる（再現性）

使い方:
    python validate_generated_data.py [--dir PATH] [--no-plots]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class Checker:
    """PASS/FAIL を集計する簡易チェッカ。"""

    def __init__(self) -> None:
        self.n_pass = 0
        self.n_fail = 0

    def check(self, name: str, ok: bool, detail: str = "") -> None:
        mark = "PASS" if ok else "FAIL"
        if ok:
            self.n_pass += 1
        else:
            self.n_fail += 1
        line = f"[{mark}] {name}"
        if detail:
            line += f"  ({detail})"
        print(line)

    def summary(self) -> bool:
        print("-" * 60)
        print(f"合計: {self.n_pass} PASS / {self.n_fail} FAIL")
        return self.n_fail == 0


def load_config_seed(cfg_path: Path) -> int:
    import yaml
    with open(cfg_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)["random_seed"]


def validate(root: Path, make_plots: bool) -> bool:
    ck = Checker()
    spectra_dir = root / "spectra"
    gt_dir = root / "ground_truth"

    meta = pd.read_csv(root / "metadata.csv")
    peak_df = pd.read_csv(gt_dir / "true_peak_parameters.csv")

    # --- スペクトル数 ---
    spec_files = sorted(spectra_dir.glob("PET_*.csv"))
    ck.check("スペクトルが 50 個存在する", len(spec_files) == 50,
             f"found {len(spec_files)}")

    # --- 各スペクトルの中身 ---
    all_201 = True
    all_monotonic = True
    all_finite = True
    all_nonneg = True
    for fp in spec_files:
        df = pd.read_csv(fp)
        if len(df) != 201:
            all_201 = False
        be = df["binding_energy_eV"].to_numpy()
        if not np.all(np.diff(be) < 0):
            all_monotonic = False
        if not np.all(np.isfinite(df.to_numpy())):
            all_finite = False
        if (df["intensity_cps"] < 0).any():
            all_nonneg = False
    ck.check("各スペクトルが 201 点ある", all_201)
    ck.check("結合エネルギーが単調に減少する", all_monotonic)
    ck.check("欠損値・無限値がない", all_finite)
    ck.check("強度が 0 以上である", all_nonneg)

    # --- 真のピーク面積比の合計が 1 ---
    frac_sum = peak_df.groupby("sample_id")["area_fraction"].sum()
    ck.check("真のピーク面積比の合計が 1", np.allclose(frac_sum.values, 1.0, atol=1e-6),
             f"min={frac_sum.min():.6f}, max={frac_sum.max():.6f}")

    # --- 処理強度と酸素含有ピーク面積比の傾向 ---
    merged = peak_df.merge(meta[["sample_id", "treatment_group", "power_W"]],
                           on="sample_id")
    oxygen = merged[merged["component"].isin(["C-O", "O-C=O"])]
    grp_ox = (oxygen.groupby("treatment_group")["area_fraction"].sum()
              .groupby(level=0).mean())
    # 群平均の酸素含有比
    ox_by_grp = (oxygen.groupby(["treatment_group", "sample_id"])["area_fraction"]
                 .sum().groupby("treatment_group").mean())
    order = ["G1", "G2", "G3", "G4", "G5"]
    ox_vals = [ox_by_grp[g] for g in order]
    increasing = all(ox_vals[i] < ox_by_grp[order[i + 1]] for i in range(len(order) - 1))
    ck.check("処理が強いほど酸素含有ピーク比が平均的に増加",
             increasing,
             " -> ".join(f"{g}:{v:.3f}" for g, v in zip(order, ox_vals)))

    # C-C/C-H は逆に減少
    cc = merged[merged["component"] == "C-C_C-H"]
    cc_by_grp = cc.groupby("treatment_group")["area_fraction"].mean()
    cc_vals = [cc_by_grp[g] for g in order]
    decreasing = all(cc_vals[i] > cc_vals[i + 1] for i in range(len(order) - 1))
    ck.check("処理が強いほど C-C/C-H 比が平均的に減少", decreasing,
             " -> ".join(f"{g}:{v:.3f}" for g, v in zip(order, cc_vals)))

    # --- 同一条件内のばらつき ---
    within_std = (oxygen.groupby(["treatment_group", "sample_id"])["area_fraction"]
                  .sum().groupby("treatment_group").std())
    ck.check("同一条件内にも十分なばらつきがある",
             (within_std > 0.003).all(),
             f"min within-group std={within_std.min():.4f}")

    # --- 再現性 ---
    seed = load_config_seed(root / "config" / "generation_config.yaml")
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        gen_script = root / "scripts" / "generate_xps_data.py"
        subprocess.run(
            [sys.executable, str(gen_script), "--outdir", str(tmp_path)],
            check=True, capture_output=True,
        )
        df_a = pd.read_csv(spectra_dir / "PET_001.csv")
        df_b = pd.read_csv(tmp_path / "spectra" / "PET_001.csv")
        reproducible = np.allclose(df_a.to_numpy(), df_b.to_numpy())
        # 真値ファイルも一致するか
        gt_a = pd.read_csv(gt_dir / "true_peak_parameters.csv")
        gt_b = pd.read_csv(tmp_path / "ground_truth" / "true_peak_parameters.csv")
        gt_ok = gt_a.equals(gt_b)
    ck.check(f"同じ乱数シード({seed})で同じデータを再生成できる",
             reproducible and gt_ok)

    # --- プレビュー画像 ---
    if make_plots:
        make_preview(root, meta, ox_by_grp, order)
        print("preview/*.png を生成しました。")

    return ck.summary()


def make_preview(root: Path, meta: pd.DataFrame, ox_by_grp, order) -> None:
    """代表スペクトル・異常スペクトル・群別面積比のプレビュー画像を作る。"""
    spectra_dir = root / "spectra"
    prev = root / "preview"
    prev.mkdir(exist_ok=True)
    anomaly = pd.read_csv(root / "ground_truth" / "anomaly_labels.csv")

    # 1) 各群 1 試料ずつの正常スペクトル
    fig, ax = plt.subplots(figsize=(7, 5))
    reps = {g: meta[meta["treatment_group"] == g]["sample_id"].iloc[0] for g in order}
    for g in order:
        sid = reps[g]
        df = pd.read_csv(spectra_dir / f"{sid}.csv")
        ax.plot(df["binding_energy_eV"], df["intensity_cps"], label=f"{g} ({sid})")
    ax.set_xlabel("Binding energy [eV]")
    ax.set_ylabel("Intensity [cps]")
    ax.set_title("Example C1s spectra by treatment group")
    ax.invert_xaxis()  # XPS 慣習: 高 BE を左に
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(prev / "example_normal_spectra.png", dpi=120)
    plt.close(fig)

    # 2) 異常スペクトル
    fig, ax = plt.subplots(figsize=(7, 5))
    for _, r in anomaly.iterrows():
        df = pd.read_csv(spectra_dir / f"{r['sample_id']}.csv")
        ax.plot(df["binding_energy_eV"], df["intensity_cps"],
                label=f"{r['sample_id']} ({r['anomaly_type']})")
    ax.set_xlabel("Binding energy [eV]")
    ax.set_ylabel("Intensity [cps]")
    ax.set_title("Anomalous spectra")
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(prev / "example_anomaly_spectra.png", dpi=120)
    plt.close(fig)

    # 3) 群別の真の面積比（積み上げ）
    peak_df = pd.read_csv(root / "ground_truth" / "true_peak_parameters.csv")
    m = peak_df.merge(meta[["sample_id", "treatment_group"]], on="sample_id")
    pivot = (m.groupby(["treatment_group", "component"])["area_fraction"].mean()
             .unstack())
    comp_order = ["C-C_C-H", "C-O", "O-C=O"]
    pivot = pivot.loc[order, comp_order]
    fig, ax = plt.subplots(figsize=(7, 5))
    bottom = np.zeros(len(order))
    colors = ["#4C72B0", "#DD8452", "#55A868"]
    for comp, c in zip(comp_order, colors):
        ax.bar(order, pivot[comp].values, bottom=bottom, label=comp, color=c)
        bottom += pivot[comp].values
    ax.set_ylabel("Mean true area fraction")
    ax.set_xlabel("Treatment group")
    ax.set_title("True area fraction by treatment group")
    ax.legend()
    fig.tight_layout()
    fig.savefig(prev / "true_area_fraction_by_group.png", dpi=120)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--dir", type=Path, default=root)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    ok = validate(args.dir, make_plots=not args.no_plots)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
