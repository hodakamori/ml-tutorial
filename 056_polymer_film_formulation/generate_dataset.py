"""ポリマーフィルム配合物の仮想データセット生成スクリプト。

社内講習会「生成AIを使ったMI」デモ1（整形済み実験データから Lasso 回帰モデルを
作成する）向けの、合成（架空）データセットを作成する。

特徴:
    * 全 500 レコード
    * ポリマーフィルムの配合物（ベースポリマー + 添加剤）が対象
    * 物質名は実在の材料名を使用し、PubChem 等から SMILES を取得できることを想定
      （デモではその SMILES からフィンガープリントを計算して説明変数にする）
    * 配合温度・配合時間などのプロセス条件を含む
    * 目的変数は 引張強度・透明度・ガスバリア性（酸素透過度）
    * 説明変数を多数用意している。目的変数に効くもの・効かないもの（ノイズ）が混在し、
      Lasso による変数選択（不要な係数をゼロにする挙動）を体験できるようにしてある。
    * 一定割合の欠損値を含む

出力: polymer_film_formulation.csv

注意: 数値は材料科学的に「それらしい」オーダーで生成した完全な架空データであり、
      実測値ではない。教材デモ専用。
"""

import numpy as np
import pandas as pd

N_RECORDS = 500
SEED = 20260712

rng = np.random.default_rng(SEED)


# ---------------------------------------------------------------------------
# ベースポリマー: 代表的な物性のベースライン値・標準加工温度域・代表的な分子量域
#   tensile : 引張強度 [MPa]
#   transp  : 透明度 [%]
#   logbar  : 酸素透過度の常用対数 log10(cc*mm/(m^2*day*atm))  小さいほど高バリア
#   temp    : 標準的な配合(溶融混練)温度域 [degC]
#   mw      : 重量平均分子量 Mw の代表域 [g/mol]
# ---------------------------------------------------------------------------
BASE_POLYMERS = {
    "Polyethylene":                 dict(tensile=15, transp=86, logbar=np.log10(1900), temp=(150, 200), mw=(80_000, 250_000)),
    "Polypropylene":                dict(tensile=33, transp=88, logbar=np.log10(1400), temp=(180, 230), mw=(150_000, 400_000)),
    "Polystyrene":                  dict(tensile=45, transp=90, logbar=np.log10(2500), temp=(190, 240), mw=(150_000, 350_000)),
    "Poly(vinyl chloride)":         dict(tensile=52, transp=80, logbar=np.log10(120),  temp=(160, 200), mw=(60_000, 130_000)),
    "Poly(methyl methacrylate)":    dict(tensile=66, transp=92, logbar=np.log10(400),  temp=(200, 240), mw=(80_000, 200_000)),
    "Polycarbonate":                dict(tensile=64, transp=89, logbar=np.log10(1000), temp=(240, 290), mw=(25_000, 60_000)),
    "Poly(lactic acid)":            dict(tensile=55, transp=90, logbar=np.log10(550),  temp=(170, 210), mw=(80_000, 200_000)),
    "Poly(ethylene terephthalate)": dict(tensile=80, transp=89, logbar=np.log10(50),   temp=(250, 290), mw=(30_000, 80_000)),
    "Poly(vinyl alcohol)":          dict(tensile=70, transp=92, logbar=np.log10(3),    temp=(180, 220), mw=(30_000, 130_000)),
    "Nylon 6":                      dict(tensile=76, transp=85, logbar=np.log10(30),   temp=(230, 270), mw=(20_000, 40_000)),
}
# 出現頻度の重み（汎用樹脂を多めに）
POLYMER_WEIGHTS = np.array([1.6, 1.6, 1.2, 1.1, 1.0, 0.9, 1.2, 1.3, 0.8, 0.9])


# ---------------------------------------------------------------------------
# 主添加剤: カテゴリ・代表的な配合量域(wt%)・各物性への 1wt% あたりの効果
#   d_tensile : 引張強度への寄与 [MPa / wt%]
#   d_transp  : 透明度への寄与   [% / wt%]
#   d_logbar  : log10(酸素透過度) への寄与 [/ wt%]  （負 = バリア向上）
# ---------------------------------------------------------------------------
ADDITIVES = {
    # 可塑剤: 柔軟化で引張強度を下げ、バリアはやや悪化
    "Dioctyl phthalate":                dict(cat="plasticizer", load=(3, 30),  d_tensile=-0.70, d_transp=+0.05, d_logbar=+0.010),
    "Dibutyl phthalate":                dict(cat="plasticizer", load=(3, 25),  d_tensile=-0.75, d_transp=+0.05, d_logbar=+0.012),
    "Acetyl tributyl citrate":          dict(cat="plasticizer", load=(3, 25),  d_tensile=-0.60, d_transp=+0.06, d_logbar=+0.008),
    "Triacetin":                        dict(cat="plasticizer", load=(2, 20),  d_tensile=-0.55, d_transp=+0.04, d_logbar=+0.010),
    "Glycerol":                         dict(cat="plasticizer", load=(2, 18),  d_tensile=-0.65, d_transp=-0.10, d_logbar=+0.015),
    # 充填剤/フィラー: 透明度を大きく下げ、迂回効果でバリア向上
    "Titanium dioxide":                 dict(cat="filler",      load=(1, 20),  d_tensile=+0.30, d_transp=-2.20, d_logbar=-0.020),
    "Calcium carbonate":                dict(cat="filler",      load=(5, 30),  d_tensile=+0.15, d_transp=-1.60, d_logbar=-0.015),
    "Talc":                             dict(cat="filler",      load=(3, 25),  d_tensile=+0.25, d_transp=-1.80, d_logbar=-0.025),
    # 酸化防止剤: 物性への影響は小さい
    "Butylated hydroxytoluene":         dict(cat="antioxidant", load=(0.1, 1.0), d_tensile=+0.05, d_transp=-0.10, d_logbar=0.000),
    "Irganox 1010":                     dict(cat="antioxidant", load=(0.1, 1.0), d_tensile=+0.08, d_transp=-0.08, d_logbar=0.000),
    # 紫外線吸収剤: 透明度をやや低下
    "Benzophenone":                     dict(cat="uv_absorber", load=(0.2, 3.0), d_tensile=-0.05, d_transp=-0.40, d_logbar=+0.002),
    "2-Hydroxy-4-methoxybenzophenone":  dict(cat="uv_absorber", load=(0.2, 3.0), d_tensile=-0.04, d_transp=-0.35, d_logbar=+0.002),
    # 滑剤/スリップ剤: ほぼ中立
    "Stearic acid":                     dict(cat="lubricant",   load=(0.2, 3.0), d_tensile=-0.10, d_transp=-0.05, d_logbar=+0.003),
    "Erucamide":                        dict(cat="lubricant",   load=(0.1, 2.0), d_tensile=-0.08, d_transp=-0.05, d_logbar=+0.004),
}
ADDITIVE_NAMES = list(ADDITIVES.keys())


# ---------------------------------------------------------------------------
# 第2添加剤（相容化剤・造核剤など）。引張強度に弱く効く。"None" は無添加。
#   d_tensile : 引張強度への寄与 [MPa / wt%]
# ---------------------------------------------------------------------------
SECOND_ADDITIVES = {
    "None":              dict(load=(0.0, 0.0), d_tensile=0.00),
    "Maleic anhydride":  dict(load=(0.3, 3.0), d_tensile=+0.90),   # 相容化剤
    "Sorbitol":          dict(load=(0.1, 2.0), d_tensile=+0.30),   # 造核・透明化剤
    "Pentaerythritol":   dict(load=(0.1, 2.0), d_tensile=+0.20),
    "Zinc stearate":     dict(load=(0.1, 2.0), d_tensile=+0.05),   # 潤滑・安定剤（ほぼ中立）
}
SECOND_ADDITIVE_NAMES = list(SECOND_ADDITIVES.keys())
SECOND_ADDITIVE_WEIGHTS = np.array([0.40, 0.18, 0.16, 0.14, 0.12])  # 4割は無添加

# ノイズ用のカテゴリ列（目的変数に一切影響しない）
OPERATORS = ["Op-A", "Op-B", "Op-C", "Op-D", "Op-E", "Op-F"]
LINES = ["Line-1", "Line-2", "Line-3"]


def sample_records(n):
    polymer_names = list(BASE_POLYMERS.keys())
    p = POLYMER_WEIGHTS / POLYMER_WEIGHTS.sum()
    p2 = SECOND_ADDITIVE_WEIGHTS / SECOND_ADDITIVE_WEIGHTS.sum()

    rows = []
    for i in range(n):
        poly = rng.choice(polymer_names, p=p)
        pp = BASE_POLYMERS[poly]

        add = rng.choice(ADDITIVE_NAMES)
        ap = ADDITIVES[add]
        loading = float(rng.uniform(*ap["load"]))

        sec = rng.choice(SECOND_ADDITIVE_NAMES, p=p2)
        sp = SECOND_ADDITIVES[sec]
        sec_loading = 0.0 if sec == "None" else float(rng.uniform(*sp["load"]))

        # ベースポリマーの分子量 Mw [g/mol]
        mw_lo, mw_hi = pp["mw"]
        mw = float(rng.uniform(mw_lo, mw_hi))
        mw_norm = (mw - mw_lo) / (mw_hi - mw_lo)         # 0 .. 1

        # 配合温度 [degC]: 標準域の中央付近を狙いつつ、上下にばらつく
        t_lo, t_hi = pp["temp"]
        t_mid = 0.5 * (t_lo + t_hi)
        t_width = 0.5 * (t_hi - t_lo)
        temperature = float(rng.normal(t_mid, t_width * 0.55))

        # 各種プロセス条件
        comp_time = float(rng.uniform(2, 30))            # 配合(混練)時間 [min]
        screw_speed = float(rng.uniform(50, 400))        # スクリュー回転数 [rpm]
        cooling_rate = float(rng.uniform(1, 100))        # 冷却速度 [degC/min]
        draw_ratio = float(rng.uniform(1.0, 5.0))        # 延伸倍率 [-]
        filler_size = float(rng.uniform(0.1, 20.0))      # フィラー平均粒径 [um]
        moisture = float(rng.uniform(50, 5000))          # 残留水分 [ppm]
        film_thickness = float(rng.uniform(10, 300))     # フィルム厚 [um]（透過"係数"には効かない）

        # ノイズ用（目的変数に無関係な列）
        operator = rng.choice(OPERATORS)
        line = rng.choice(LINES)
        lab_humidity = float(rng.uniform(20, 70))        # 試験室湿度 [%]
        ambient_temp = float(rng.uniform(18, 30))        # 室温 [degC]
        qc_batch = f"Q-{int(rng.integers(1, 60)):03d}"
        noise_index = float(rng.uniform(0, 1))

        # --- 温度が最適域からずれると分散不良・熱劣化で物性低下 ---
        dev = (temperature - t_mid) / t_width            # -約2 .. +約2
        temp_penalty_tensile = -6.0 * dev**2
        temp_penalty_transp = -3.0 * dev**2

        # --- 混練時間: 長いほど分散が進み物性向上（頭打ち）、過剰では微減 ---
        time_factor = 1.0 - np.exp(-comp_time / 8.0)     # 0 .. ~1
        over_mix = max(0.0, comp_time - 24.0)
        time_tensile = 7.0 * time_factor - 0.4 * over_mix
        time_transp = 3.5 * time_factor

        # --- 主添加剤の配合量効果 ---
        add_tensile = ap["d_tensile"] * loading
        add_transp = ap["d_transp"] * loading
        add_logbar = ap["d_logbar"] * loading

        # --- 追加の説明変数の効果 ---
        mw_tensile = 4.0 * (mw_norm - 0.5)               # 高分子量ほど強い
        sec_tensile = sp["d_tensile"] * sec_loading      # 第2添加剤
        screw_tensile = 0.010 * (screw_speed - 200.0)    # 分散性（弱い）
        screw_transp = 0.006 * (screw_speed - 200.0)
        draw_tensile = 3.0 * (draw_ratio - 1.0)          # 延伸配向で強化
        draw_logbar = -0.06 * (draw_ratio - 1.0)         # 配向でバリア向上
        cool_transp = 0.030 * (cooling_rate - 30.0)      # 急冷=低結晶=高透明
        cool_logbar = +0.0020 * (cooling_rate - 30.0)    # 低結晶=バリア悪化
        # フィラー使用時のみ、粒径が小さいほど分散良好で透明・強度に寄与
        is_filler = ap["cat"] == "filler"
        filler_tensile = (3.0 - filler_size) * 0.30 if is_filler else 0.0
        filler_transp = (3.0 - filler_size) * 0.50 if is_filler else 0.0
        # 親水性ポリマーのみ残留水分でバリア悪化
        moisture_logbar = 5e-5 * moisture if poly in ("Poly(vinyl alcohol)", "Nylon 6") else 0.0

        # 引張強度 [MPa]
        tensile = (pp["tensile"] + add_tensile + temp_penalty_tensile + time_tensile
                   + mw_tensile + sec_tensile + screw_tensile + draw_tensile
                   + filler_tensile + rng.normal(0, 2.5))
        tensile = max(1.0, tensile)

        # 透明度 [%]
        transparency = (pp["transp"] + add_transp + temp_penalty_transp + time_transp
                        + screw_transp + cool_transp + filler_transp + rng.normal(0, 1.8))
        transparency = float(np.clip(transparency, 1.0, 99.0))

        # ガスバリア性: 酸素透過度 [cc*mm/(m^2*day*atm)]（小さいほど高バリア）
        logbar = (pp["logbar"] + add_logbar + 0.05 * max(0.0, dev)
                  + draw_logbar + cool_logbar + moisture_logbar + rng.normal(0, 0.10))
        oxygen_permeability = float(10.0 ** logbar)

        rows.append({
            "sample_id": f"PF-{i + 1:04d}",
            # --- 物質名（SMILES 取得対象）---
            "base_polymer": poly,
            "additive": add,
            "second_additive": sec,
            # --- 配合条件 ---
            "additive_loading_wt%": round(loading, 2),
            "second_additive_loading_wt%": round(sec_loading, 2),
            "base_polymer_mw_g_per_mol": round(mw, 0),
            "filler_particle_size_um": round(filler_size, 2),
            # --- プロセス条件 ---
            "compounding_temperature_C": round(temperature, 1),
            "compounding_time_min": round(comp_time, 1),
            "screw_speed_rpm": round(screw_speed, 0),
            "cooling_rate_C_per_min": round(cooling_rate, 1),
            "draw_ratio": round(draw_ratio, 2),
            "annealing_moisture_ppm": round(moisture, 0),
            "film_thickness_um": round(film_thickness, 1),
            # --- 目的変数に無関係なノイズ列 ---
            "operator": operator,
            "production_line": line,
            "lab_humidity_percent": round(lab_humidity, 1),
            "ambient_temperature_C": round(ambient_temp, 1),
            "qc_batch": qc_batch,
            "noise_index": round(noise_index, 4),
            # --- 目的変数 ---
            "tensile_strength_MPa": round(tensile, 2),
            "transparency_percent": round(transparency, 2),
            "oxygen_permeability_cc_mm_m2_day_atm": round(oxygen_permeability, 3),
        })
    return pd.DataFrame(rows)


def inject_missing(df):
    """一部の列に欠損値を導入する（デモで前処理を体験させるため）。

    物質名（base_polymer）は SMILES 取得に必須なので欠損させない。
    additive はまれに欠損、数値の説明変数・目的変数に数%の欠損を入れる。
    """
    n = len(df)

    def mask(frac):
        return rng.random(n) < frac

    for col, frac in {
        "additive": 0.02,
        "additive_loading_wt%": 0.06,
        "second_additive_loading_wt%": 0.05,
        "base_polymer_mw_g_per_mol": 0.05,
        "filler_particle_size_um": 0.06,
        "compounding_temperature_C": 0.05,
        "compounding_time_min": 0.07,
        "screw_speed_rpm": 0.05,
        "cooling_rate_C_per_min": 0.05,
        "draw_ratio": 0.05,
        "annealing_moisture_ppm": 0.06,
        "film_thickness_um": 0.04,
        "lab_humidity_percent": 0.04,
        "ambient_temperature_C": 0.04,
        "tensile_strength_MPa": 0.08,
        "transparency_percent": 0.07,
        "oxygen_permeability_cc_mm_m2_day_atm": 0.07,
    }.items():
        df.loc[mask(frac), col] = np.nan
    return df


def main():
    df = sample_records(N_RECORDS)
    df = inject_missing(df)

    out = "polymer_film_formulation.csv"
    df.to_csv(out, index=False)

    n = len(df)
    print(f"wrote {out}: {n} records, {df.shape[1]} columns")
    print("\n--- columns ---")
    print(list(df.columns))
    print("\n--- missing value ratio ---")
    print((df.isna().mean() * 100).round(1).astype(str) + " %")
    print("\n--- head ---")
    print(df.head(5).to_string(index=False))


if __name__ == "__main__":
    main()
