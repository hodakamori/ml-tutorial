# ポリマーフィルム配合物 仮想データセット（デモ1用）

社内講習会「生成AIを使ったMI」**デモ1**（整形済み実験データから Lasso 回帰モデルを
作成し、専用MIツール相当の解析を行う）で使う、**架空（合成）** のデータセットです。

- ファイル: `polymer_film_formulation.csv`
- レコード数: **500**
- 生成スクリプト: `generate_dataset.py`（`python generate_dataset.py` で再生成可能・乱数シード固定）
- 参考: `smiles_reference.csv`（物質名→SMILES の対応表。デモの答え合わせ用）

> ⚠️ 数値は材料科学的に「それらしい」オーダーで作った**完全な架空データ**であり、実測値では
> ありません。教材デモ専用です。

## 想定するデモの流れ

1. CSV を読み込み、行数・欠損値・分布を確認する（整形済み実験データ）。
2. `base_polymer` と `additive` の**物質名から SMILES を取得**する（PubChem 等）。
3. SMILES から**フィンガープリント**（例: Morgan/RDKit FP）を計算し、説明変数にする。
4. プロセス条件（配合温度・時間・添加剤配合量）と結合して特徴量行列を作る。
5. 欠損値を処理し、**Lasso 回帰**で引張強度・透明度・ガスバリア性を予測、重要因子を確認する。

## 列の説明（データ辞書）

| 列名 | 種別 | 単位 | 説明 |
| --- | --- | --- | --- |
| `sample_id` | ID | – | サンプル識別子（`PF-0001` 形式） |
| `base_polymer` | 説明変数（物質名） | – | ベースポリマー名。SMILES 取得の対象 |
| `additive` | 説明変数（物質名） | – | 添加剤名。SMILES 取得の対象 |
| `additive_loading_wt%` | 説明変数 | wt% | 添加剤の配合量 |
| `compounding_temperature_C` | 説明変数（プロセス） | ℃ | 配合（溶融混練）温度 |
| `compounding_time_min` | 説明変数（プロセス） | min | 配合（混練）時間 |
| `tensile_strength_MPa` | 目的変数 | MPa | 引張強度 |
| `transparency_percent` | 目的変数 | % | 透明度（全光線透過率イメージ、0–100） |
| `oxygen_permeability_cc_mm_m2_day_atm` | 目的変数 | cc·mm/(m²·day·atm) | ガスバリア性 = 酸素透過度。**小さいほど高バリア** |

## 欠損値について

前処理を体験できるよう、一部の列に欠損値（`NaN`）を含みます（各列おおむね数%程度）。
`base_polymer` は SMILES 取得に必須なため欠損させていません。`additive` はまれに欠損します。

## 収録している物質

- **ベースポリマー（10種）**: Polyethylene, Polypropylene, Polystyrene, Poly(vinyl chloride),
  Poly(methyl methacrylate), Polycarbonate, Poly(lactic acid), Poly(ethylene terephthalate),
  Poly(vinyl alcohol), Nylon 6
- **添加剤（14種）**: 可塑剤（Dioctyl phthalate, Dibutyl phthalate, Acetyl tributyl citrate,
  Triacetin, Glycerol）、充填剤（Titanium dioxide, Calcium carbonate, Talc）、
  酸化防止剤（Butylated hydroxytoluene, Irganox 1010）、
  紫外線吸収剤（Benzophenone, 2-Hydroxy-4-methoxybenzophenone）、
  滑剤（Stearic acid, Erucamide）

いずれも実在の材料・化合物名で、PubChem 等から SMILES を取得できることを想定しています。
ポリマー名は代表構造（繰り返し単位/モノマー）の SMILES に解決されるのが一般的です。
`smiles_reference.csv` に対応表を用意しています（デモで自力取得した結果の照合用）。

## データに埋め込んだ傾向（ネタバレ）

Lasso で意味のある係数が出るよう、次のような関係を（ノイズ付きで）埋め込んでいます。

- 引張強度・透明度・ガスバリア性は主に**ベースポリマー種**で決まる。
- **可塑剤**は引張強度を下げる。**充填剤**は透明度を大きく下げ、迂回効果でバリアを向上させる。
- 配合温度が各ポリマーの**最適域から外れる**と、分散不良・熱劣化で引張強度・透明度が低下する。
- **混練時間**が長いほど分散が進み物性が向上する（頭打ち、過混練で微減）。
