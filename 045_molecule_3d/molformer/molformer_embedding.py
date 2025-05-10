import types, sys, importlib.machinery

if 'apex' not in sys.modules:
    apex_pkg = types.ModuleType('apex')
    apex_pkg.__spec__  = importlib.machinery.ModuleSpec('apex', loader=None, is_package=True)
    apex_pkg.__path__  = []
    sys.modules['apex'] = apex_pkg

    amp_mod = types.ModuleType('apex.amp')
    amp_mod.__spec__ = importlib.machinery.ModuleSpec('apex.amp', loader=None, is_package=False)
    sys.modules['apex.amp'] = amp_mod

    optim_mod = types.ModuleType('apex.optimizers')
    optim_mod.__spec__ = importlib.machinery.ModuleSpec('apex.optimizers', loader=None, is_package=False)
    apex_pkg.optimizers = optim_mod
    sys.modules['apex.optimizers'] = optim_mod

from argparse import Namespace
import yaml
from tokenizer.tokenizer import MolTranBertTokenizer
from train_pubchem_light import LightningModule
import pandas as pd
import torch
from fast_transformers.masking import LengthMask as LM


def batch_split(data, batch_size=64):
    i = 0
    while i < len(data):
        yield data[i:min(i+batch_size, len(data))]
        i += batch_size

def embed(model, smiles, tokenizer, batch_size=64):
    model.eval()
    embeddings = []
    for batch in batch_split(smiles, batch_size=batch_size):
        batch_enc = tokenizer.batch_encode_plus(batch, padding=True, add_special_tokens=True)
        idx, mask = torch.tensor(batch_enc['input_ids']), torch.tensor(batch_enc['attention_mask'])
        with torch.no_grad():
            token_embeddings = model.blocks(model.tok_emb(idx), length_mask=LM(mask.sum(-1)))
        # average pooling over tokens
        input_mask_expanded = mask.unsqueeze(-1).expand(token_embeddings.size()).float()
        sum_embeddings = torch.sum(token_embeddings * input_mask_expanded, 1)
        sum_mask = torch.clamp(input_mask_expanded.sum(1), min=1e-9)
        embedding = sum_embeddings / sum_mask
        embeddings.append(embedding.detach().cpu())
    return torch.cat(embeddings)


with open('../../data/Pretrained MoLFormer/hparams.yaml', 'r') as f:
    config = Namespace(**yaml.safe_load(f))

tokenizer = MolTranBertTokenizer('bert_vocab.txt')
ckpt = '../../data/Pretrained MoLFormer/checkpoints/N-Step-Checkpoint_3_30000.ckpt'
lm = LightningModule(config, tokenizer.vocab).load_from_checkpoint(ckpt, config=config, vocab=tokenizer.vocab)

df = pd.read_csv("freesolv.csv")
smiles = df[" SMILES"]
X = embed(lm, smiles, tokenizer).numpy()
X_df = pd.DataFrame(X, columns=[f"molformer_{i}" for i in range(X.shape[1])])
X_df["SMILES"] = smiles
X_df.to_csv("molformer_embed.csv", index=False)