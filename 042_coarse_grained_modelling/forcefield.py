from typing import Any, Dict

import openmm
import openmm.app as app
from openmm import unit


class CGForceField:
    def __init__(self, param_dict: Dict[str, Any]):
        self.param_dict = param_dict

    @classmethod
    def from_dict(cls, param_dict: Dict[str, Any]):
        return cls(param_dict)

    def create_system(self, topology: app.Topology) -> openmm.System:
        system = openmm.System()
        atoms = list(topology.atoms())

        bond_force = openmm.HarmonicBondForce()
        system.addForce(bond_force)

        nonbonded_force = openmm.NonbondedForce()
        system.addForce(nonbonded_force)

        nonbonded_params = self.param_dict.get("NONBONDED", {})
        for i, atom in enumerate(atoms):
            res_label = atom.residue.name
            nb_data = nonbonded_params[res_label]
            mass_value = nb_data["mass"]
            mass_amu = mass_value * unit.amu
            system.addParticle(mass_amu)

            sigma_nm = nb_data["sigma"]  # nm
            epsilon_kj = nb_data["epsilon"]  # kJ/mol
            charge_e = nb_data["charge"]  # elementary charge

            nonbonded_force.addParticle(
                charge_e,
                sigma_nm * unit.nanometer,
                epsilon_kj * unit.kilojoule_per_mole,
            )

        bond_params = self.param_dict.get("BOND", {})
        for bond in topology.bonds():
            atom1, atom2 = bond
            i1 = atom1.index
            i2 = atom2.index
            label1 = atom1.residue.name
            label2 = atom2.residue.name

            key1 = (label1, label2)
            key2 = (label2, label1)

            bp = None
            if key1 in bond_params:
                bp = bond_params[key1]
            elif key2 in bond_params:
                bp = bond_params[key2]

            if bp is not None:
                r0 = bp["r0"]  # nm
                k_val = bp["k"]  # kJ/mol/nm^2
            else:
                raise ValueError(f"Unknown bond: {label1}-{label2}")

            bond_force.addBond(
                i1,
                i2,
                r0 * unit.nanometer,
                k_val * unit.kilojoule_per_mole / (unit.nanometer**2),
            )

        # -------------------
        # 4) Add angle or torsion force (if present)
        # -------------------

        # -------------------
        # 5) Set nonbonded interaction
        # -------------------
        nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

        return system
