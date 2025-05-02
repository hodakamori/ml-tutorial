import argparse
import torch
import numpy as np

from torch_geometric.data import Data
from torch_geometric.nn import radius_graph
from torch_geometric.nn.models import SchNet, DimeNet, DimeNetPlusPlus
from torch_geometric.datasets import QM9  # QM9 dataset import
from rdkit import Chem
from rdkit.Chem import AllChem


def _embed_mol(mol: Chem.Mol, max_attempts: int = 10) -> Chem.Mol:
    """Add 3D coordinates with RDKit and optimize with UFF (exception on failure)"""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0  # reproducibility
    for _ in range(max_attempts):
        if AllChem.EmbedMolecule(mol, params) == 0:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            return mol
    raise RuntimeError("3D embedding failed.")


def smiles_to_data(
    smiles: str,
    cutoff: float = 5.0,
    loop: bool = False,
    device: torch.device | str = "cpu",
) -> Data:
    """
    Create Data(z, pos, edge_index) from SMILES.
    DimeNet series automatically generates radius_graph → RBF → angle features → triplets internally.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = _embed_mol(mol)
    conf = mol.GetConformer()

    zs, pos = [], []
    for atom in mol.GetAtoms():
        zs.append(atom.GetAtomicNum())
        p = conf.GetAtomPosition(atom.GetIdx())
        pos.append([p.x, p.y, p.z])

    z = torch.tensor(zs, dtype=torch.long, device=device)
    pos = torch.tensor(pos, dtype=torch.float, device=device)
    edge_index = radius_graph(pos, r=cutoff, loop=loop)

    return Data(z=z, pos=pos, edge_index=edge_index)


def main():
    parser = argparse.ArgumentParser(
        description="SMILES → GNN inference (SchNet / DimeNet / DimeNet++)"
    )
    parser.add_argument("smiles", help="Input SMILES")
    parser.add_argument(
        "--model",
        choices=["schnet", "dimenet", "dimenetpp"],
        default="schnet",
        help="Model to use",
    )
    parser.add_argument(
        "--cutoff", type=float, default=5.0, help="radius_graph cutoff [Å]"
    )
    parser.add_argument(
        "--device", type=str, default="cpu", help="CPU or CUDA (e.g. cuda:0)"
    )
    args = parser.parse_args()

    # 1) SMILES → Data
    data = smiles_to_data(
        args.smiles, cutoff=args.cutoff, loop=False, device=args.device
    )

    # 2) Prepare QM9 dataset (for pretrained model)
    root_qm9 = "./data/QM9"
    dataset = QM9(root_qm9)

    # 3) Model initialization (load all pretrained models)
    if args.model == "schnet":
        model, *_ = SchNet.from_qm9_pretrained(
            root=root_qm9,
            dataset=dataset,
            target=7,  # Target index (example: 7 = Egap)
        )
    elif args.model == "dimenet":
        model, *_ = DimeNet.from_qm9_pretrained(
            root=root_qm9,
            dataset=dataset,
            target=7,
        )
    elif args.model == "dimenetpp":
        model, *_ = DimeNetPlusPlus.from_qm9_pretrained(
            root=root_qm9,
            dataset=dataset,
            target=7,
        )
    else:
        raise ValueError(f"Unknown model: {args.model}")

    # 4) Inference
    # Get batch information
    batch = data.batch if hasattr(data, "batch") else None
    # Input (z, pos, batch) for all models
    out = model(data.z, data.pos, batch)

    print(f"{args.model} prediction:", out.detach().cpu())

if __name__ == "__main__":
    main()
