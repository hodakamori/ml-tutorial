import torch
import torch.nn as nn
from torch_geometric.data import Data
from torch_geometric.datasets import QM9
from torch_geometric.nn.models import SchNet, DimeNet, DimeNetPlusPlus
from torch_geometric.nn import radius_graph, global_add_pool, global_mean_pool
from typing import Literal

# -------------------------------------------------
# Utility
# -------------------------------------------------

def smiles_to_data(
    smiles: str,
    cutoff: float = 5.0,
    loop: bool = False,
    device: str = "cpu",
) -> Data:
    """Wrap external SMILES→Data conversion function."""
    from pyg_gnn import smiles_to_data as _builder

    return _builder(smiles, cutoff=cutoff, loop=loop, device=device)


# -------------------------------------------------
# Base class for pooling helpers
# -------------------------------------------------

def _maybe_pool(x: torch.Tensor, batch: torch.Tensor, mode: Literal["add", "mean", None]):
    if mode == "add":
        return global_add_pool(x, batch)
    if mode == "mean":
        return global_mean_pool(x, batch)
    return x


# -------------------------------------------------
# SchNet extractor
# -------------------------------------------------

class SchNetEmbeddingExtractor(nn.Module):
    """Return atom‑level embedding h0 from a QM9‑pretrained SchNet."""

    def __init__(self, model: SchNet, pool: Literal["add", "mean", None] = None):
        super().__init__()
        self.model = model.eval()
        self.pool = pool

    @classmethod
    def from_qm9_pretrained(
        cls,
        root: str,
        dataset: QM9,
        target: int,
        pool: Literal["add", "mean", None] = None,
    ) -> "SchNetEmbeddingExtractor":
        model, *_ = SchNet.from_qm9_pretrained(root=root, dataset=dataset, target=target)
        return cls(model, pool)

    @torch.no_grad()
    def forward(self, data: Data):
        z = data.z
        batch = getattr(data, "batch", None)
        if batch is None:
            batch = torch.zeros(z.size(0), dtype=torch.long, device=z.device)
        h0 = self.model.embedding(z)
        return _maybe_pool(h0, batch, self.pool)


# -------------------------------------------------
# DimeNet / DimeNet++ extractor (emb call)
# -------------------------------------------------

class _BaseDimeNetExtractor(nn.Module):
    """Common logic for DimeNet/DimeNet++ embedding extraction.

    * The `emb` method outputs *edge*-level embeddings (|E|, hidden).
    * Convert them to *node*-level embeddings (|V|, hidden) by
      averaging over target node `j`.
    """

    def __init__(self, model, pool):
        super().__init__()
        self.model = model.eval()
        self.pool = pool

    @torch.no_grad()
    def forward(self, data: Data):
        from torch_scatter import scatter_mean

        z, pos = data.z, data.pos
        batch = getattr(data, "batch", None)
        if batch is None:
            batch = torch.zeros(z.size(0), dtype=torch.long, device=z.device)

        edge_index = radius_graph(
            pos,
            r=self.model.cutoff,
            batch=batch,
            max_num_neighbors=self.model.max_num_neighbors,
        )
        i, j = edge_index  # source -> target (as in DimeNet paper)
        dist = (pos[i] - pos[j]).pow(2).sum(dim=-1).sqrt()
        rbf = self.model.rbf(dist)

        # Edge‑level embedding
        e_emb = self.model.emb(z, rbf, i, j)  # (|E|, hidden)

        # Convert to node‑level by averaging over incoming edges (target j)
        v_emb = scatter_mean(e_emb, j, dim=0, dim_size=z.size(0))  # (|V|, hidden)

        return _maybe_pool(v_emb, batch, self.pool)


class DimeNetEmbeddingExtractor(_BaseDimeNetExtractor):
    """QM9‑pretrained DimeNet embedding extractor."""

    @classmethod
    def from_qm9_pretrained(
        cls,
        root: str,
        dataset: QM9,
        target: int,
        pool: Literal["add", "mean", None] = None,
    ) -> "DimeNetEmbeddingExtractor":
        model, *_ = DimeNet.from_qm9_pretrained(root=root, dataset=dataset, target=target)
        return cls(model, pool)


class DimeNetPPEmbeddingExtractor(_BaseDimeNetExtractor):
    """QM9‑pretrained DimeNet++ embedding extractor."""

    @classmethod
    def from_qm9_pretrained(
        cls,
        root: str,
        dataset: QM9,
        target: int,
        pool: Literal["add", "mean", None] = None,
    ) -> "DimeNetPPEmbeddingExtractor":
        model, *_ = DimeNetPlusPlus.from_qm9_pretrained(
            root=root, dataset=dataset, target=target
        )
        return cls(model, pool)


# -------------------------------------------------
# Example usage
# -------------------------------------------------

if __name__ == "__main__":
    root_qm9 = "./data/QM9"
    dataset = QM9(root_qm9)
    data = smiles_to_data("CC(=O)O")

    sne = SchNetEmbeddingExtractor.from_qm9_pretrained(root_qm9, dataset, 7, pool="mean")
    print("SchNet h0:", sne(data).shape)

    dne = DimeNetEmbeddingExtractor.from_qm9_pretrained(root_qm9, dataset, 7, pool="mean")
    print("DimeNet pooled:", dne(data).shape)

    dpp = DimeNetPPEmbeddingExtractor.from_qm9_pretrained(root_qm9, dataset, 7, pool="mean")
    print("DimeNet++ pooled:", dpp(data).shape)
