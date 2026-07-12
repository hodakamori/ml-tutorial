# ポリマーフィルム配合物 仮想データセット（デモ1用）

社内講習会「生成AIを使ったMI」**デモ1**（整形済み実験データから Lasso 回帰モデルを
作成し、専用MIツール相当の解析を行う）で使う、**架空（合成）** のデータセットです。

- ファイル: `polymer_film_formulation.csv`
- レコード数: **500** / カラム数: **24**
- 生成スクリプト: `generate_dataset.py`（`python generate_dataset.py` で再生成可能・乱数シード固定）
- 参考: `smiles_reference.csv`（物質名→SMILES の対応表。デモの答え合わせ用）

> ⚠️ 数値は材料科学的に「それらしい」オーダーで作った**完全な架空データ**であり、実測値では
> ありません。教材デモ専用です。

## 想定するデモの流れ

1. CSV を読み込み、行数・欠損値・分布を確認する（整形済み実験データ）。
2. `base_polymer` / `additive` / `second_additive` の**物質名から SMILES を取得**する（PubChem 等）。
3. SMILES から**フィンガープリント**（例: Morgan/RDKit FP）を計算し、説明変数にする。
4. プロセス条件・配合条件と結合して特徴量行列を作る。
5. 欠損値を処理し、**Lasso 回帰**で引張強度・透明度・ガスバリア性を予測、重要因子を確認する。

## 列の説明（データ辞書）

説明変数は多めに用意してあり、**目的変数に効く列と効かない（ノイズ）列が混在**しています。
Lasso が不要な特徴量の係数をゼロに縮小し、変数選択する挙動を体験できます。

### ID
| 列名 | 単位 | 説明 |
| --- | --- | --- |
| `sample_id` | – | サンプル識別子（`PF-0001` 形式） |

### 物質名（SMILES 取得の対象）
| 列名 | 単位 | 説明 |
| --- | --- | --- |
| `base_polymer` | – | ベースポリマー名（10種） |
| `additive` | – | 主添加剤名（14種） |
| `second_additive` | – | 第2添加剤名（相容化剤・造核剤など。`None`=無添加を含む） |

### 配合条件（説明変数）
| 列名 | 単位 | 効き | 説明 |
| --- | --- | --- | --- |
| `additive_loading_wt%` | wt% | ○ | 主添加剤の配合量 |
| `second_additive_loading_wt%` | wt% | △ | 第2添加剤の配合量（引張強度に弱く効く） |
| `base_polymer_mw_g_per_mol` | g/mol | ○ | ベースポリマーの重量平均分子量 Mw（高いほど強い） |
| `filler_particle_size_um` | µm | △ | フィラー平均粒径（充填剤使用時のみ効く。小さいほど分散良好） |

### プロセス条件（説明変数）
| 列名 | 単位 | 効き | 説明 |
| --- | --- | --- | --- |
| `compounding_temperature_C` | ℃ | ○ | 配合（溶融混練）温度。最適域を外れると物性低下 |
| `compounding_time_min` | min | ○ | 配合（混練）時間。長いほど分散向上（頭打ち） |
| `screw_speed_rpm` | rpm | △ | スクリュー回転数（分散性に弱く効く） |
| `cooling_rate_C_per_min` | ℃/min | ○ | 冷却速度。急冷=低結晶=高透明・低バリア |
| `draw_ratio` | – | ○ | 延伸倍率。配向で引張強度とバリアが向上 |
| `annealing_moisture_ppm` | ppm | △ | 残留水分（親水性ポリマーのバリアに効く） |
| `film_thickness_um` | µm | ✕ | フィルム厚（透過"係数"には効かない＝ほぼノイズ） |

### 目的変数に無関係なノイズ列（変数選択の練習用）
| 列名 | 単位 | 効き | 説明 |
| --- | --- | --- | --- |
| `operator` | – | ✕ | 作業者ID |
| `production_line` | – | ✕ | 製造ライン |
| `lab_humidity_percent` | % | ✕ | 試験室湿度 |
| `ambient_temperature_C` | ℃ | ✕ | 室温 |
| `qc_batch` | – | ✕ | 品質管理バッチ番号 |
| `noise_index` | – | ✕ | 意味のない乱数値 |

### 目的変数
| 列名 | 単位 | 説明 |
| --- | --- | --- |
| `tensile_strength_MPa` | MPa | 引張強度 |
| `transparency_percent` | % | 透明度（0–100） |
| `oxygen_permeability_cc_mm_m2_day_atm` | cc·mm/(m²·day·atm) | ガスバリア性 = 酸素透過度。**小さいほど高バリア** |

（効き列: ○=主要因子 / △=弱い・条件付き / ✕=目的変数に無関係）

## 欠損値について

前処理を体験できるよう、多くの列に欠損値（`NaN`）を含みます（各列おおむね数〜10%程度）。
`base_polymer` は SMILES 取得に必須なため欠損させていません。`additive` はまれに欠損します。

## 収録している物質

- **ベースポリマー（10種）**: Polyethylene, Polypropylene, Polystyrene, Poly(vinyl chloride),
  Poly(methyl methacrylate), Polycarbonate, Poly(lactic acid), Poly(ethylene terephthalate),
  Poly(vinyl alcohol), Nylon 6
- **主添加剤（14種）**: 可塑剤（Dioctyl phthalate, Dibutyl phthalate, Acetyl tributyl citrate,
  Triacetin, Glycerol）、充填剤（Titanium dioxide, Calcium carbonate, Talc）、
  酸化防止剤（Butylated hydroxytoluene, Irganox 1010）、
  紫外線吸収剤（Benzophenone, 2-Hydroxy-4-methoxybenzophenone）、
  滑剤（Stearic acid, Erucamide）
- **第2添加剤（5種）**: None（無添加）, Maleic anhydride, Sorbitol, Pentaerythritol, Zinc stearate

いずれも実在の材料・化合物名で、PubChem 等から SMILES を取得できることを想定しています。
ポリマー名は代表構造（繰り返し単位/モノマー）の SMILES に解決されるのが一般的です。
`smiles_reference.csv` に対応表を用意しています（デモで自力取得した結果の照合用）。

## 参考: Lasso での学習可能性（乱数シード固定時）

物質名のフィンガープリント代わりにカテゴリを one-hot したうえで Lasso 回帰にかけると、
5分割交差検証の R² はおおよそ 引張強度 0.94 / 透明度 0.83 / 酸素透過度 0.97。
`noise_index` や `film_thickness_um` などの無関係列は係数がほぼ 0 に縮小されます。

## データに埋め込んだ傾向（ネタバレ）

- 引張強度・透明度・ガスバリア性は主に**ベースポリマー種**で決まる。
- **可塑剤**は引張強度を下げる。**充填剤**は透明度を大きく下げ、迂回効果でバリアを向上させる。
- 配合温度が各ポリマーの**最適域から外れる**と、分散不良・熱劣化で引張強度・透明度が低下する。
- **混練時間**・**スクリュー回転数**が大きいほど分散が進み物性が向上する（頭打ち・過混練で微減）。
- **分子量**が大きいほど引張強度が高い。**延伸倍率**が大きいほど引張強度・バリアが向上する。
- **冷却速度**が速いほど結晶化が抑えられ、透明度は上がるがバリアは下がる。
- `operator` / `production_line` / `lab_humidity_percent` / `ambient_temperature_C` /
  `qc_batch` / `noise_index` / `film_thickness_um` は目的変数に**影響しない**（Lasso で係数≈0）。
