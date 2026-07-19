"""PET フィルム C1s XPS ダミースペクトル生成スクリプト。

社内講習会「生成AIを使ったMI」デモ3向けの、合成（架空）XPS データセットを作成する。
生の測定データから MI 用特徴量を生成するデモの入力として設計している。

生成物:
    * spectra/PET_001.csv ... PET_050.csv : 解析対象の個別スペクトル
        （binding_energy_eV, intensity_cps の 2 列）
    * metadata.csv                        : 解析コードから参照してよいメタ情報
    * ground_truth/true_peak_parameters.csv     : 各ピークの真値
    * ground_truth/true_spectrum_parameters.csv : スペクトル全体の真値
    * ground_truth/anomaly_labels.csv           : 異常データの正解ラベル

真値（ground_truth/）は生成に用いた真のパラメータであり、
デモ中のフィッティングコードからは参照させない。フィッティング結果と後から
定量比較するための答え合わせ用データである。

生成条件は config/generation_config.yaml に分離してあり、乱数シードにより完全に
再現できる。数値は「それらしい」オーダーで作った完全な架空データであり実測値ではない。

使い方:
    python generate_xps_data.py [--config PATH] [--outdir PATH]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# ピーク形状: 面積で規格化した pseudo-Voigt 関数
# ---------------------------------------------------------------------------
def pseudo_voigt(x: np.ndarray, center: float, fwhm: float, eta: float,
                 area: float) -> np.ndarray:
    """面積規格化 pseudo-Voigt。

    Gauss 成分と Lorentz 成分をそれぞれ面積 1 に規格化し、eta で線形結合する。
    したがって返り値を x で積分するとおよそ ``area`` になる。面積を直接
    パラメータに与えることで、生成時の真の面積とフィッティング面積を比較しやすくする。

    Args:
        x: 結合エネルギー配列 [eV]
        center: ピーク中心 [eV]
        fwhm: 半値全幅 [eV]（Gauss/Lorentz 共通の FWHM）
        eta: Lorentzian 混合率（0=純 Gauss, 1=純 Lorentz）
        area: ピーク面積
    """
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    gauss = np.exp(-((x - center) ** 2) / (2.0 * sigma ** 2)) / (
        sigma * math.sqrt(2.0 * math.pi)
    )
    gamma = fwhm / 2.0
    lorentz = (gamma / math.pi) / ((x - center) ** 2 + gamma ** 2)
    return area * (eta * lorentz + (1.0 - eta) * gauss)


def shirley_like_background(energy: np.ndarray, peak_signal: np.ndarray,
                            ratio: float) -> np.ndarray:
    """Shirley 型に近い背景。

    低結合エネルギー側ほど背景が高くなる（光電子の非弾性散乱による）挙動を、
    ピーク信号の累積積分で近似する。energy は降順（高 -> 低）で渡される想定。

    ステップの高さは「ピーク最大高さ × ratio」とする。累積和そのものを使うと
    段差がピーク面積のオーダー（= 面積/刻み幅）になり、低 BE 側が非物理的に
    持ち上がってしまうため、0->1 に正規化した累積割合でスケールする。
    """
    # 高結合エネルギー側からの累積和（降順配列なので index 増加 = 低 BE 方向）
    cumulative = np.cumsum(peak_signal)
    total = cumulative[-1]
    peak_max = peak_signal.max()
    if total <= 0 or peak_max <= 0:
        return np.zeros_like(energy)
    step_height = ratio * peak_max
    return step_height * (cumulative / total)


# ---------------------------------------------------------------------------
# 設定読み込みとメタデータ生成
# ---------------------------------------------------------------------------
def load_config(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_metadata(cfg: dict) -> pd.DataFrame:
    """50 試料分のメタデータ（sample_id, 処理条件, 反復番号）を作る。"""
    rows = []
    sample_no = 0
    for grp in cfg["treatment_groups"]:
        for rep in range(1, grp["n_replicates"] + 1):
            sample_no += 1
            rows.append(
                dict(
                    sample_id=f"PET_{sample_no:03d}",
                    treatment_group=grp["group"],
                    treatment=grp["treatment"],
                    power_W=grp["power_W"],
                    time_s=grp["time_s"],
                    replicate=rep,
                )
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 真のパラメータ生成
# ---------------------------------------------------------------------------
def sample_true_parameters(cfg: dict, meta: pd.DataFrame,
                           rng: np.random.Generator,
                           anomaly_map: dict) -> tuple[list, list]:
    """試料ごとに真のピーク・スペクトルパラメータを決定する。

    Returns:
        (peak_records, spectrum_records) の 2 リスト。
    """
    peaks_cfg = cfg["peaks"]
    ps = cfg["peak_shape"]
    grp_lookup = {g["group"]: g for g in cfg["treatment_groups"]}

    peak_records = []
    spectrum_records = []

    for _, row in meta.iterrows():
        sid = row["sample_id"]
        grp = grp_lookup[row["treatment_group"]]
        anomaly = anomaly_map.get(sid)

        # --- 面積比: 群の中心値まわりに Dirichlet ばらつき ---
        center_frac = np.asarray(grp["area_fraction_center"], dtype=float)
        conc = cfg["area_fraction_concentration"]
        frac = rng.dirichlet(center_frac * conc)  # 合計 1 に正規化済み

        # --- 全ピーク面積（試料ごとに ±~18% 変動） ---
        total_area = cfg["total_area"]["base"] * (
            1.0 + rng.normal(0.0, cfg["total_area"]["variation_frac"])
        )
        total_area = max(total_area, cfg["total_area"]["base"] * 0.3)

        # 低信号異常: 総面積を約 30% に
        if anomaly and anomaly["type"] == "low_signal":
            total_area *= cfg["anomaly_params"]["low_signal_factor"]

        # --- 共通エネルギーシフト（帯電） ---
        if anomaly and anomaly["type"] == "charging_shift":
            lo, hi = cfg["charging_shift"]["anomaly_shift_eV"]
            mag = rng.uniform(lo, hi)
            energy_shift = mag * rng.choice([-1.0, 1.0])
        else:
            energy_shift = rng.normal(0.0, cfg["charging_shift"]["normal_sigma_eV"])

        # --- 3 ピークに共通させる FWHM のベース（弱い相関のため） ---
        fwhm_lo, fwhm_hi = ps["fwhm_eV"]
        fwhm_base = rng.uniform(fwhm_lo, fwhm_hi)
        corr = ps["fwhm_correlation"]

        jit_lo, jit_hi = ps["center_jitter_eV"]
        center_sigma = rng.uniform(jit_lo, jit_hi)

        for i, pk in enumerate(peaks_cfg):
            # 位置: 基準値 + 帯電シフト + 個別ゆらぎ
            center = pk["center_eV"] + energy_shift + rng.normal(0.0, center_sigma)
            # FWHM: 共通ベースと個別値を corr で混合（弱い相関）
            fwhm_indiv = rng.uniform(fwhm_lo, fwhm_hi)
            fwhm = corr * fwhm_base + (1.0 - corr) * fwhm_indiv
            eta = rng.uniform(*ps["eta"])
            area = total_area * frac[i]
            peak_records.append(
                dict(
                    sample_id=sid,
                    component=pk["name"],
                    center_eV=round(center, 4),
                    area=round(area, 3),
                    area_fraction=round(float(frac[i]), 5),
                    fwhm_eV=round(fwhm, 4),
                    eta=round(eta, 4),
                )
            )

        # --- バックグラウンド係数 ---
        bg = cfg["background"]
        offset = rng.uniform(*bg["offset"])
        slope = rng.uniform(*bg["linear_slope"])
        shirley = rng.uniform(*bg["shirley_ratio"])
        quad = rng.uniform(*bg["quadratic"])

        # --- ノイズ強度 ---
        nz = cfg["noise"]
        count_scale = rng.uniform(*nz["count_scale"])
        gauss_sigma = rng.uniform(*nz["gaussian_sigma"])
        if anomaly and anomaly["type"] == "high_noise":
            factor = cfg["anomaly_params"]["high_noise_factor"]
            count_scale *= factor
            gauss_sigma *= factor

        spectrum_records.append(
            dict(
                sample_id=sid,
                total_peak_area=round(total_area, 3),
                energy_shift_eV=round(float(energy_shift), 4),
                bg_offset=round(offset, 3),
                bg_linear_slope=round(slope, 4),
                bg_shirley_ratio=round(shirley, 5),
                bg_quadratic=round(quad, 5),
                noise_count_scale=round(count_scale, 4),
                noise_gaussian_sigma=round(gauss_sigma, 3),
                area_frac_CC_CH=round(float(frac[0]), 5),
                area_frac_CO=round(float(frac[1]), 5),
                area_frac_OCO=round(float(frac[2]), 5),
            )
        )

    return peak_records, spectrum_records


# ---------------------------------------------------------------------------
# スペクトル合成
# ---------------------------------------------------------------------------
def synthesize_spectrum(energy: np.ndarray, sid: str, peak_df: pd.DataFrame,
                        spec_row: pd.Series, cfg: dict,
                        anomaly: dict | None,
                        rng: np.random.Generator) -> np.ndarray:
    """1 試料分の観測強度を合成する。"""
    # 1) 3 個のピーク
    peak_signal = np.zeros_like(energy)
    for _, pk in peak_df.iterrows():
        peak_signal += pseudo_voigt(
            energy, pk["center_eV"], pk["fwhm_eV"], pk["eta"], pk["area"]
        )

    # 未モデル化ショルダー異常: 285.5 eV 付近に弱い第 4 成分を追加
    if anomaly and anomaly["type"] == "unmodeled_shoulder":
        ap = cfg["anomaly_params"]
        shoulder_area = spec_row["total_peak_area"] * ap["shoulder_area_frac"]
        peak_signal += pseudo_voigt(
            energy,
            ap["shoulder_center_eV"] + spec_row["energy_shift_eV"],
            ap["shoulder_fwhm_eV"],
            0.3,
            shoulder_area,
        )

    # 2) バックグラウンド（定数 + 一次 + Shirley 型 + 二次）
    e0 = energy.max()
    de = energy - e0  # 高 BE 端を基準とした相対エネルギー
    background = (
        spec_row["bg_offset"]
        + spec_row["bg_linear_slope"] * de
        + spec_row["bg_quadratic"] * de ** 2
    )
    background += shirley_like_background(energy, peak_signal,
                                          spec_row["bg_shirley_ratio"])

    intensity = peak_signal + background

    # 3) ベースラインドリフト（ごく緩やかな正弦的うねり）
    drift_amp = rng.uniform(5.0, 20.0)
    drift_phase = rng.uniform(0, 2 * math.pi)
    drift = drift_amp * np.sin(2 * math.pi * (de / (e0 - energy.min())) + drift_phase)
    intensity = intensity + drift

    # 4) ノイズ（カウントノイズ + Gaussian ノイズ）
    count_noise = rng.normal(0.0, 1.0, size=energy.shape) * (
        spec_row["noise_count_scale"] * np.sqrt(np.clip(intensity, 0, None))
    )
    gauss_noise = rng.normal(0.0, spec_row["noise_gaussian_sigma"], size=energy.shape)
    intensity = intensity + count_noise + gauss_noise

    # 5) スパイク異常: 数点に大きな外れ値
    if anomaly and anomaly["type"] == "spike":
        ap = cfg["anomaly_params"]
        idx = rng.choice(len(energy), size=ap["spike_count"], replace=False)
        amp = rng.uniform(*ap["spike_amplitude"], size=ap["spike_count"])
        intensity[idx] += amp

    # 6) 負の強度を 0 以上に補正
    intensity = np.clip(intensity, 0.0, None)
    return intensity


# ---------------------------------------------------------------------------
# メイン
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve().parent
    root = here.parent
    parser.add_argument("--config", type=Path,
                        default=root / "config" / "generation_config.yaml")
    parser.add_argument("--outdir", type=Path, default=root)
    args = parser.parse_args()

    cfg = load_config(args.config)
    rng = np.random.default_rng(cfg["random_seed"])

    outdir = args.outdir
    spectra_dir = outdir / "spectra"
    gt_dir = outdir / "ground_truth"
    spectra_dir.mkdir(parents=True, exist_ok=True)
    gt_dir.mkdir(parents=True, exist_ok=True)

    # 結合エネルギー軸（降順）
    n = cfg["energy"]["n_points"]
    energy = np.round(
        np.linspace(cfg["energy"]["start_eV"], cfg["energy"]["stop_eV"], n), 2
    )

    # 異常マップ（sample_index=1-based -> PET_xxx）
    anomaly_map = {
        f"PET_{a['sample_index']:03d}": a for a in cfg["anomalies"]
    }

    # メタデータ
    meta = build_metadata(cfg)

    # 真のパラメータ
    peak_records, spec_records = sample_true_parameters(cfg, meta, rng, anomaly_map)
    peak_df_all = pd.DataFrame(peak_records)
    spec_df_all = pd.DataFrame(spec_records)

    # スペクトル合成 & 保存
    for _, mrow in meta.iterrows():
        sid = mrow["sample_id"]
        pk = peak_df_all[peak_df_all["sample_id"] == sid]
        srow = spec_df_all[spec_df_all["sample_id"] == sid].iloc[0]
        anomaly = anomaly_map.get(sid)
        intensity = synthesize_spectrum(energy, sid, pk, srow, cfg, anomaly, rng)
        spec_csv = pd.DataFrame(
            {"binding_energy_eV": energy, "intensity_cps": np.round(intensity, 2)}
        )
        spec_csv.to_csv(spectra_dir / f"{sid}.csv", index=False)

    # メタデータ・真値の保存
    meta.to_csv(outdir / "metadata.csv", index=False)
    peak_df_all.to_csv(gt_dir / "true_peak_parameters.csv", index=False)
    spec_df_all.to_csv(gt_dir / "true_spectrum_parameters.csv", index=False)

    anomaly_df = pd.DataFrame(
        [
            dict(sample_id=f"PET_{a['sample_index']:03d}",
                 anomaly_type=a["type"], description=a["description"])
            for a in cfg["anomalies"]
        ]
    ).sort_values("sample_id")
    anomaly_df.to_csv(gt_dir / "anomaly_labels.csv", index=False)

    print(f"生成完了: {len(meta)} スペクトル -> {spectra_dir}")
    print(f"  metadata.csv, ground_truth/*.csv を出力しました。")


if __name__ == "__main__":
    main()
