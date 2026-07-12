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
# ベースポリマー: 代表的な物性のベースライン値と代表的な加工温度域
#   tensile : 引張強度 [MPa]
#   transp  : 透明度 [%]
#   logbar  : 酸素透過度の常用対数 log10(cc*mm/(m^2*day*atm))  小さいほど高バリア
#   temp    : 標準的な配合(溶融混練)温度域 [degC]
# ---------------------------------------------------------------------------
BASE_POLYMERS = {
    "Polyethylene":                 dict(tensile=15, transp=86, logbar=np.log10(1900), temp=(150, 200)),
    "Polypropylene":                dict(tensile=33, transp=88, logbar=np.log10(1400), temp=(180, 230)),
    "Polystyrene":                  dict(tensile=45, transp=90, logbar=np.log10(2500), temp=(190, 240)),
    "Poly(vinyl chloride)":         dict(tensile=52, transp=80, logbar=np.log10(120),  temp=(160, 200)),
    "Poly(methyl methacrylate)":    dict(tensile=66, transp=92, logbar=np.log10(400),  temp=(200, 240)),
    "Polycarbonate":                dict(tensile=64, transp=89, logbar=np.log10(1000), temp=(240, 290)),
    "Poly(lactic acid)":            dict(tensile=55, transp=90, logbar=np.log10(550),  temp=(170, 210)),
    "Poly(ethylene terephthalate)": dict(tensile=80, transp=89, logbar=np.log10(50),   temp=(250, 290)),
    "Poly(vinyl alcohol)":          dict(tensile=70, transp=92, logbar=np.log10(3),    temp=(180, 220)),
    "Nylon 6":                      dict(tensile=76, transp=85, logbar=np.log10(30),   temp=(230, 270)),
}
# 出現頻度の重み（汎用樹脂を多めに）
POLYMER_WEIGHTS = np.array([1.6, 1.6, 1.2, 1.1, 1.0, 0.9, 1.2, 1.3, 0.8, 0.9])


# ---------------------------------------------------------------------------
# 添加剤: カテゴリ・代表的な配合量域(wt%)・各物性への 1wt% あたりの効果
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


def sample_records(n):
    polymer_names = list(BASE_POLYMERS.keys())
    p = POLYMER_WEIGHTS / POLYMER_WEIGHTS.sum()

    rows = []
    for i in range(n):
        poly = rng.choice(polymer_names, p=p)
        pp = BASE_POLYMERS[poly]

        add = rng.choice(ADDITIVE_NAMES)
        ap = ADDITIVES[add]
        loading = float(rng.uniform(*ap["load"]))

        # 配合温度: 標準域の中央付近を狙いつつ、上下にばらつく
        t_lo, t_hi = pp["temp"]
        t_mid = 0.5 * (t_lo + t_hi)
        t_width = 0.5 * (t_hi - t_lo)
        temperature = float(rng.normal(t_mid, t_width * 0.55))

        # 配合(混練)時間 [min]
        comp_time = float(rng.uniform(2, 30))

        # --- 温度が最適域からずれると分散不良・熱劣化で物性低下 ---
        dev = (temperature - t_mid) / t_width          # -約2 .. +約2
        temp_penalty_tensile = -6.0 * dev**2
        temp_penalty_transp = -3.0 * dev**2

        # --- 混練時間: 長いほど分散が進み物性向上（頭打ち）、過剰では微減 ---
        time_factor = 1.0 - np.exp(-comp_time / 8.0)   # 0 .. ~1
        over_mix = max(0.0, comp_time - 24.0)           # 過混練
        time_tensile = 7.0 * time_factor - 0.4 * over_mix
        time_transp = 3.5 * time_factor

        # --- 添加剤の配合量効果 ---
        add_tensile = ap["d_tensile"] * loading
        add_transp = ap["d_transp"] * loading
        add_logbar = ap["d_logbar"] * loading

        # 引張強度 [MPa]
        tensile = (pp["tensile"] + add_tensile + temp_penalty_tensile
                   + time_tensile + rng.normal(0, 2.5))
        tensile = max(1.0, tensile)

        # 透明度 [%]
        transparency = (pp["transp"] + add_transp + temp_penalty_transp
                        + time_transp + rng.normal(0, 1.8))
        transparency = float(np.clip(transparency, 1.0, 99.0))

        # ガスバリア性: 酸素透過度 [cc*mm/(m^2*day*atm)]（小さいほど高バリア）
        logbar = pp["logbar"] + add_logbar + 0.05 * max(0.0, dev) + rng.normal(0, 0.10)
        oxygen_permeability = float(10.0 ** logbar)

        rows.append({
            "sample_id": f"PF-{i + 1:04d}",
            "base_polymer": poly,
            "additive": add,
            "additive_loading_wt%": round(loading, 2),
            "compounding_temperature_C": round(temperature, 1),
            "compounding_time_min": round(comp_time, 1),
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

    df.loc[mask(0.02), "additive"] = np.nan
    df.loc[mask(0.06), "additive_loading_wt%"] = np.nan
    df.loc[mask(0.05), "compounding_temperature_C"] = np.nan
    df.loc[mask(0.07), "compounding_time_min"] = np.nan
    df.loc[mask(0.08), "tensile_strength_MPa"] = np.nan
    df.loc[mask(0.07), "transparency_percent"] = np.nan
    df.loc[mask(0.07), "oxygen_permeability_cc_mm_m2_day_atm"] = np.nan
    return df


def main():
    df = sample_records(N_RECORDS)
    df = inject_missing(df)

    out = "polymer_film_formulation.csv"
    df.to_csv(out, index=False)

    n = len(df)
    print(f"wrote {out}: {n} records, {df.shape[1]} columns")
    print("\n--- missing value ratio ---")
    print((df.isna().mean() * 100).round(1).astype(str) + " %")
    print("\n--- head ---")
    print(df.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
