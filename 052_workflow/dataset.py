from dataclasses import dataclass, field
from typing import List, Optional, Literal
from rdkit import Chem
from rdkit.Chem import AllChem
import psi4
import numpy as np
import pandas as pd
import pickle
import os
import resp

@dataclass
class DFTConfig:
    basis: Literal["6-31g", "sto-3g"]
    func: Literal["scf", "b3lyp"]
    memory: int = 500 # in unit of MB
    output_prefix: str = "output"

@dataclass
class Molecule:
    name: str # propane
    smiles: str # CCC
    elements: Optional[List[str]] = None # ["C", "C", "C", "H", ...]
    positions: Optional[List[List[float]]] = None # [[0, 0.1, 0.2], ...]
    dipole: Optional[List[float]] = None # [0, 0, 0.1]
    dipole_magnitude: Optional[float] = None # 0.1
    partial_charge: Optional[List[float]] =None # [0.1, -0.1, ...]
    net_charge: Optional[int] = 0
    spin: Optional[int] = 1
    dft_config: DFTConfig = field(default_factory=lambda: DFTConfig(basis="6-31g", func="scf"))

    def __post_init__(self):
        rdmol = Chem.MolFromSmiles(self.smiles)
        rdmol = Chem.AddHs(rdmol)
        self.elements = [atom.GetSymbol() for atom in rdmol.GetAtoms()]
        AllChem.EmbedMolecule(rdmol)
        self.positions = rdmol.GetConformer(0).GetPositions()

    def _to_psi4geom(self):
        psi4geom_input = []
        psi4geom_input.append(f"{self.net_charge} {self.spin}")
        for elem, pos in zip(self.elements, self.positions):
            psi4geom_input.append(f"{elem} {pos[0]} {pos[1]} {pos[2]}")
        psi4geom = psi4.geometry("\n".join(psi4geom_input))
        return psi4geom

    def calc_dipole(self):
        psi4.core.set_output_file(f"{self.dft_config.output_prefix}_dipole.dat", False)
        psi4.set_memory(f"{self.dft_config.memory} MB")

        molecule = self._to_psi4geom()

        opt_energy, opt_wfn, history = psi4.optimize(
            f"{self.dft_config.func}/{self.dft_config.basis}", return_wfn=True, return_history=True
        )
        psi4.oeprop(opt_wfn, "DIPOLE")
        dipole = psi4.variable("SCF DIPOLE")
        self.dipole = dipole
        self.dipole_magnitude = np.linalg.norm(dipole)

    def calc_resp_charge(self):
        psi4.core.set_output_file(f"{self.dft_config.output_prefix}_resp.dat", False)
        psi4.set_memory(f"{self.dft_config.memory} MB")

        molecule = self._to_psi4geom()
        opt_energy, opt_wfn, history = psi4.optimize(
            f"{self.dft_config.func}/{self.dft_config.basis}", return_wfn=True, return_history=True
        )
        mol = opt_wfn.molecule()
        options = {
            "VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0],
            "VDW_POINT_DENSITY": 1.0,
            "RESP_A": 0.0005,
            "RESP_B": 0.1,
        }

        # Call for first stage fit
        charges = resp.resp([mol], options)
        self.partial_charge = charges[1]

    def to_dict(self):
        return {
            "SMILES": self.smiles,
            "dipole": self.dipole,
            "dipole_magnitude": self.dipole_magnitude,
            "partial_charge": self.partial_charge
        }
    
    def to_pickle(self):
        with open(f"{self.name}.pk", "wb") as f:
            pickle.dump(self, f)


dataset = []
for name, smiles in zip(["water", "methanol", "methane"], ["O", "CO", "C"]):
    if os.path.exists(f"{name}.pk"):
        with open(f"{name}.pk", "rb") as f:
            mol = pickle.load(f)
    else:
        mol = Molecule(name=name, smiles=smiles)
    if mol.dipole is None:
        mol.calc_dipole()
    if mol.partial_charge is None:
        mol.calc_resp_charge()
    mol.to_pickle()
    # aws s3 cp {name}.pk s3://***/{name}.pk
    dataset.append(mol.to_dict())

df = pd.DataFrame(dataset)
df["partial_charge_max"] = df["partial_charge"].apply(max)
df["partial_charge_min"] = df["partial_charge"].apply(min)
df["partial_charge_delta"] = df["partial_charge_max"] - df["partial_charge_min"]

print(df)