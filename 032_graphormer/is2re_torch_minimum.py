import numpy as np
from graphormer.tasks.is2re import LMDBDataset, PBCDataset, AtomDataset, KeywordDataset, NestedDictionaryDataset
from graphormer.models.graphormer_3d import Graphormer3D

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

lmdb_dataset = LMDBDataset(db_path="data.lmdb")
pbc_dataset = PBCDataset(lmdb_dataset)

atoms = AtomDataset(pbc_dataset, "atoms")
tags = KeywordDataset(pbc_dataset, "tags")
real_mask = KeywordDataset(pbc_dataset, "real_mask")
pos = KeywordDataset(pbc_dataset, "pos")
relaxed_energy = KeywordDataset(pbc_dataset, "relaxed_energy", is_scalar=True)

dataset = NestedDictionaryDataset(
    {
        "net_input": {
            "pos": pos,
            "atoms": atoms,
            "tags": tags,
            "real_mask": real_mask,
        },
        "targets": {
            "relaxed_energy": relaxed_energy,
        },
    },
    sizes=[np.zeros(len(atoms))],
)
data = dataset[0]
energy, node_output, node_target_mask = model(
    atoms = data["net_input.atoms"].unsqueeze(0),
    tags = data["net_input.tags"].unsqueeze(0),
    pos = data["net_input.pos"].unsqueeze(0),
    real_mask = data["net_input.real_mask"].unsqueeze(0)
)
print(energy)