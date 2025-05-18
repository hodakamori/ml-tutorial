import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator as rfg

import pandas as pd
import numpy as np
import glob
import pickle
import os
from typing import Sequence, Literal
from pathlib import Path
# from autogluon.tabular import TabularPredictor
# from sklearn.model_selection import train_test_split


def smiles_to_bits(smiles):
    N_BITS, RADIUS = 2048, 2
    gen = rfg.GetMorganGenerator(radius=RADIUS, fpSize=N_BITS)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.full(N_BITS, np.nan, dtype="float32")
    fp = gen.GetFingerprint(mol)
    arr = np.zeros(N_BITS, dtype="uint8")
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def _load_vec(path: Path, model: str) -> np.ndarray:
    """path の pkl から model のベクトルを float32 1-D ndarray で返す"""
    with path.open("rb") as f:
        emb_dict = pickle.load(f)
    vec = emb_dict["emb"][model]
    if hasattr(vec, "detach"):  # torch.Tensor の場合
        vec = vec.detach().cpu().numpy()
    return np.asarray(vec, dtype="float32").flatten()


def calc_nnp_emb_merge(
    df: pd.DataFrame,
    model: str,
    emb_dir="pyg/embeddings",
    id_dtype: Literal["str", "int"] = "str",  # id の型を合わせるため
) -> pd.DataFrame:
    emb_paths = sorted(Path(emb_dir).glob("*.pkl"))
    if not emb_paths:
        raise FileNotFoundError(f"No *.pkl found in {emb_dir}")

    first_vec = _load_vec(emb_paths[0], model)
    dim = first_vec.size
    emb_cols = [f"{model}_{i}" for i in range(dim)]

    ids, vecs = [], []
    for p in emb_paths:
        vec = _load_vec(p, model)
        if vec.size != dim:
            continue
        ids.append(p.stem)
        vecs.append(vec)

    emb_mat = np.vstack(vecs)
    emb_df = pd.DataFrame(emb_mat, columns=emb_cols, dtype="float32")
    emb_df.insert(0, "id", ids)

    if id_dtype == "int":
        emb_df["id"] = emb_df["id"].astype(int)
    else:
        df["id"] = df["id"].astype(str)
        emb_df["id"] = emb_df["id"].astype(str)

    df_out = (
        df.reset_index(drop=False)
        .merge(emb_df, on="id", how="left", validate="one_to_one")
        .sort_values("index")
        .drop(columns="index")
        .reset_index(drop=True)
    )
    return df_out


df = (
    pd.read_csv("freesolv.csv")
    .rename(
        columns={
            "# compound id (and file prefix)": "id",
            " SMILES": "smiles",
            " experimental value (kcal/mol)": "dG",
        }
    )
    .loc[:, ["id", "smiles", "dG"]]
)

# bit_mat = np.vstack(
#     df["smiles"].map(smiles_to_bits).to_numpy()
# )  # shape = (n_samples, 2048)

# bit_cols = [f"morgan_{i}" for i in range(N_BITS)]
# bits_df = pd.DataFrame(bit_mat, columns=bit_cols, dtype="uint8")

# df = pd.concat([df.reset_index(drop=True), bits_df], axis=1, copy=False)


# df = calc_nnp_emb_merge(df, model="dimenetpp", emb_dir="pyg/embeddings", id_dtype="str")
# df.to_csv("freesolv_dimenetpp.csv", index=False)
# print(df)

mf = pd.read_csv("molformer/molformer_embed.csv")
df = pd.merge(df, mf, left_on="smiles", right_on="SMILES", how="left")
df = df.drop(columns=["SMILES"])
df.to_csv("freesolv_molformer.csv", index=False)

print(df)
# drop_cols = ["id", "smiles"]
# data = df.drop(columns=drop_cols, errors="ignore")


# train_df, test_df = train_test_split(data, test_size=0.2, random_state=42, shuffle=True)

# predictor = TabularPredictor(label="dG", eval_metric="root_mean_squared_error").fit(
#     train_data=train_df,
#     presets="best_quality",
#     time_limit=None,
# )

# perf = predictor.evaluate(test_df)  # RMSE, MAE などが表示される
# print(perf)

# pred = predictor.predict(test_df.drop(columns=["dG"]))
# print("予測値 (kcal/mol) 上位 5 行:\n", pred.head())
