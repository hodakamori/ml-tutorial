# schnetsmiles.py
from typing import Tuple

import torch
from torch_geometric.data import Data
from torch_geometric.nn import radius_graph

from rdkit import Chem
from rdkit.Chem import AllChem


def _embed_mol(mol: Chem.Mol, max_attempts: int = 10) -> Chem.Mol:
    """RDKit で 3D 座標を付与し UFF 最適化する（失敗時は例外）"""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0  # 再現性のため
    for _ in range(max_attempts):
        if AllChem.EmbedMolecule(mol, params) == 0:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            return mol
    raise RuntimeError("3D embedding failed.")


def smiles_to_schnet_data(
    smiles: str,
    cutoff: float = 5.0,
    loop: bool = False,
    device: torch.device | str = "cpu",
) -> Data:
    """
    SMILES から SchNet 用 Data を作成するユーティリティ

    Parameters
    ----------
    smiles : str
        入力 SMILES 文字列
    cutoff : float, default 5.0
        radius_graph のカットオフ [Å]
    loop : bool, default False
        自己ループを含めるかどうか
    device : torch.device or str
        返却する Data が載るデバイス

    Returns
    -------
    data : torch_geometric.data.Data
        attr:
            z (LongTensor)   : [n_nodes] 原子番号
            pos (FloatTensor): [n_nodes, 3] 3 D 座標 (Å)
            edge_index       : [2, n_edges] 隣接行列
    """
    # --- 1. RDKit で分子を生成し 3D 埋め込み ------------------------
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = _embed_mol(mol)
    conf = mol.GetConformer()

    # --- 2. 原子番号と座標を取得 ------------------------------------
    atomic_numbers = []
    positions = []
    for atom in mol.GetAtoms():
        atomic_numbers.append(atom.GetAtomicNum())
        pos = conf.GetAtomPosition(atom.GetIdx())
        positions.append([pos.x, pos.y, pos.z])

    z = torch.tensor(atomic_numbers, dtype=torch.long, device=device)
    pos = torch.tensor(positions, dtype=torch.float, device=device)

    # --- 3. radius_graph でエッジを張る -----------------------------
    edge_index = radius_graph(pos, r=cutoff, loop=loop)

    # --- 4. Data オブジェクトへ -------------------------------------
    data = Data(z=z, pos=pos, edge_index=edge_index)
    return data


# ---------------------- 使い方例 ------------------------------------
if __name__ == "__main__":
    from torch_geometric.nn.models import SchNet

    smiles = "c1ccccc1O"  # フェノール
    cutoff = 5.0
    data = smiles_to_schnet_data(smiles, cutoff=cutoff)

    model = SchNet(
        hidden_channels=128,
        num_filters=128,
        num_interactions=3,
        num_gaussians=50,
        cutoff=cutoff,
    )

    # 分子全体のエネルギー (スカラ) を出力する例
    data = data.to("cpu")  # CUDA を使うなら "cuda"
    out = model(data.z, data.pos, data.batch if "batch" in data else None)
    print(out)  # tensor([E])
