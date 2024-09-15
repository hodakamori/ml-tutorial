from rdkit import Chem
from rdkit.Chem import AllChem
import torch
import numpy as np
from graphormer.tasks.is2re import LMDBDataset, PBCDataset, AtomDataset, KeywordDataset, NestedDictionaryDataset
from graphormer.models.graphormer_3d import Graphormer3D

mol = Chem.MolFromSmiles("CCCOC")
mol = AllChem.AddHs(mol)
AllChem.EmbedMolecule(mol)
pos = torch.tensor(mol.GetConformer(0).GetPositions()).float()
atoms = torch.tensor([atom.GetAtomicNum() for atom in mol.GetAtoms()]).long()
tags = torch.ones_like(atoms)
real_mask = torch.ones_like(atoms)

class Args:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

args = Args(
    layers=6,
    blocks=1,
    embed_dim=768,
    ffn_embed_dim=768,
    attention_heads=32,
    dropout=0.1,
    attention_dropout=0.1,
    activation_dropout=0.1,
    input_dropout=0.1,
    num_kernel=128,
    node_loss_weight=1.0,
    min_node_loss_weight=1.0
)

model = Graphormer3D.build_model(
    task="is2re",
    args=args
)

energy, node_output, node_target_mask = model(
    atoms = atoms.unsqueeze(0),
    tags = tags.unsqueeze(0),
    pos = pos.unsqueeze(0),
    real_mask = real_mask.unsqueeze(0)
)
print(energy)