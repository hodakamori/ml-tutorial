import pandas as pd
from typing import Literal
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import root_mean_squared_error
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.linear_model import Lasso, Ridge
from sklearn.svm import SVR
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
import pickle
import glob
from sklearn.preprocessing import StandardScaler


def split_ood_dG(
    df, criteria: Literal["upper", "lower"], threshold: float
) -> pd.DataFrame:
    threshold_value = df["dG"].quantile(threshold)
    if criteria == "upper":
        test = df[df["dG"] > threshold_value]
        train = df[df["dG"] <= threshold_value]
    elif criteria == "lower":
        test = df[df["dG"] < threshold_value]
        train = df[df["dG"] >= threshold_value]
    else:
        raise ValueError("criteria must be 'upper' or 'lower'")
    return train, test


def load_freesolv():
    """Load the FreeSolv dataset."""
    df = pd.read_csv("freesolv.csv")
    df = df.rename(
        columns={
            "# compound id (and file prefix)": "id",
            " SMILES": "smiles",
            " experimental value (kcal/mol)": "dG",
        }
    )
    df = df[["id", "smiles", "dG"]]
    return df


def calc_morganfp(df):
    """Calculate the RDKit fingerprints for the FreeSolv dataset."""

    def calc_fp(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return fp

    df["fp"] = df["smiles"].apply(calc_fp)
    return df


def calc_nnp_emb(df, model):
    df_out = df.copy()
    df_out["fp"] = None
    df_out["fp"] = df_out["fp"].astype("object")

    def _load_fp(row):
        pkl_files = glob.glob(f"pyg/embeddings/{row['id']}.pkl")
        if not pkl_files:
            return None
        with open(pkl_files[0], "rb") as f:
            emb = pickle.load(f)
        vec = emb["emb"][model]
        if hasattr(vec, "detach"):
            vec = vec.detach().cpu().numpy()
        return vec.flatten()

    df_out["fp"] = df_out.apply(_load_fp, axis=1)
    return df_out


def train(df, model):
    """Train the model on the FreeSolv dataset."""

    X = df["fp"].tolist()
    y = df["dG"].tolist()
    model.fit(X, y)
    return model


def evaluate(model, df_train, df_test):
    """Evaluate the model on the FreeSolv dataset."""
    X_train = df_train["fp"].tolist()
    y_train = df_train["dG"].tolist()
    X_test = df_test["fp"].tolist()
    y_test = df_test["dG"].tolist()
    rmse_train = root_mean_squared_error(y_train, model.predict(X_train))
    rmse_test = root_mean_squared_error(y_test, model.predict(X_test))
    print(f"Train RMSE: {rmse_train:.4f}")
    print(f"Test RMSE: {rmse_test:.4f}")

    return rmse_train, rmse_test


def calc_descriptors(df_train, df_test, desc_name):
    if desc_name == "morganfp":
        df_train = calc_morganfp(df_train)
        df_test = calc_morganfp(df_test)
    elif desc_name == "schnet_emb":
        df_train = calc_nnp_emb(df_train, "schnet")
        df_test = calc_nnp_emb(df_test, "schnet")
    elif desc_name == "dimenet_emb":
        df_train = calc_nnp_emb(df_train, "dimenet")
        df_test = calc_nnp_emb(df_test, "dimenet")
    elif desc_name == "dimenetpp_emb":
        df_train = calc_nnp_emb(df_train, "dimenetpp")
        df_test = calc_nnp_emb(df_test, "dimenetpp")
    else:
        raise ValueError("desc_name must be 'morganfp' or 'mordredfp'")
    df_train = df_train.dropna(subset=["fp"])
    df_test = df_test.dropna(subset=["fp"])
    sc = StandardScaler()
    X_train = df_train["fp"].tolist()
    X_test = df_test["fp"].tolist()
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)
    df_train["fp"] = list(X_train)
    df_test["fp"] = list(X_test)
    return df_train, df_test


models = {
    "lasso": Lasso(alpha=0.1, random_state=42),
    "ridge": Ridge(alpha=0.1, random_state=42),
    "rf": RandomForestRegressor(n_estimators=100, random_state=42),
    "xgb": XGBRegressor(n_estimators=100, random_state=42),
    "lgbm": LGBMRegressor(random_state=42, verbosity=-1),
}


df = load_freesolv()
df_train, df_test = split_ood_dG(df, "upper", 0.9)
print(f"Train size: {len(df_train)}")
print(f"Test size: {len(df_test)}")
results = []
for desc_name in ["morganfp", "schnet_emb", "dimenet_emb", "dimenetpp_emb"]:
    df_train, df_test = calc_descriptors(df_train, df_test, desc_name)
    for name, model in models.items():
        model = train(df_train, model)
        rmse_train, rmse_test = evaluate(model, df_train, df_test)
        results.append(
            {
                "desc_name": desc_name,
                "model_name": name,
                "rmse_train": rmse_train,
                "rmse_test": rmse_test,
            }
        )
df_results = pd.DataFrame(results)
df_results.to_csv("results.csv", index=False)
