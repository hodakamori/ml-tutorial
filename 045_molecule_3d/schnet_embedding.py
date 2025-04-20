# schnet_features.py
from typing import Literal, Optional

import torch
from torch_geometric.nn import global_add_pool, global_mean_pool
from torch_scatter import scatter_add, scatter_mean
from torch_geometric.nn.models import SchNet as _SchNetBase


class SchNetFeatureExtractor(_SchNetBase):
    """
    SchNet の forward を改変し，最後の線形層に入る前の
    hidden state h を返すラッパー．

    Parameters
    ----------
    pool : {"add", "mean", None}, default "mean"
        - "add"  : グラフごとに加算 (global_add_pool)
        - "mean" : グラフごとに平均 (global_mean_pool)
        - None   : プーリングせず各原子ベクトルそのまま返す
    """

    def __init__(self, *args, pool: Literal["add", "mean", None] = "mean", **kwargs):
        super().__init__(*args, **kwargs)
        self._feat_pool = pool

    @torch.no_grad()
    def forward(self, z, pos, batch: Optional[torch.Tensor] = None):
        batch = torch.zeros_like(z) if batch is None else batch
        h = self.embedding(z)  # [num_atoms, hidden_channels]

        edge_index, edge_weight = self.interaction_graph(pos, batch)
        edge_attr = self.distance_expansion(edge_weight)

        # ------------- interaction blocks -------------
        for interaction in self.interactions:
            h = h + interaction(h, edge_index, edge_weight, edge_attr)
        # ----------------------------------------------

        if self._feat_pool is None:
            return h  # [num_atoms, hidden_channels]

        if self._feat_pool == "add":
            return scatter_add(h, batch, dim=0)  # [num_graphs, hidden_channels]
        else:  # "mean"
            return scatter_mean(h, batch, dim=0)  # [num_graphs, hidden_channels]


from schnet import smiles_to_schnet_data

smiles = "CC(=O)O"  # 酢酸
data = smiles_to_schnet_data(smiles)

extractor = SchNetFeatureExtractor(
    hidden_channels=128,
    num_filters=128,
    num_interactions=6,
    num_gaussians=50,
    cutoff=5.0,
    pool="mean",  # グラフ全体の 128 次元特徴
).eval()

feat = extractor(data.z, data.pos)  # shape: [1, 128]
print(feat)
