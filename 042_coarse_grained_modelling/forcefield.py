from collections import defaultdict
from typing import Any, Dict, List, Tuple

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

        angle_force = openmm.HarmonicAngleForce()
        system.addForce(angle_force)

        torsion_force = openmm.PeriodicTorsionForce()
        system.addForce(torsion_force)

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
        # --- 3) AngleForce (自動生成された angles を参照)
        # param_dict["ANGLE"] は { (resA, resB, resC): {"theta0":..., "k":...}, ... }
        angle_params = self.param_dict.get("ANGLE", {})

        angles = self._get_angles(topology)  # List[(i,j,k)]

        # atoms[i].residue.name, ...
        for i, j, k in angles:
            resA = atoms[i].residue.name
            resB = atoms[j].residue.name
            resC = atoms[k].residue.name

            key1 = (resA, resB, resC)
            key2 = None
            # もし対称を考えるなら別の書き換え (ここではあえてkey2は使わない)

            ap = None
            if key1 in angle_params:
                ap = angle_params[key1]
            # elif key2 in angle_params:
            #    ap = angle_params[key2]
            if ap is not None:
                theta0 = ap["theta0"]  # rad
                k_ang = ap["k"]  # kJ/(mol * rad^2)
                angle_force.addAngle(
                    i,
                    j,
                    k,
                    theta0 * unit.radian,  # or unit.degrees if you prefer
                    k_ang * unit.kilojoule_per_mole / (unit.radian**2),
                )

        # --- 4) TorsionForce (自動生成された torsions を参照)
        # param_dict["TORSION"] は { (resA,resB,resC,resD): {"periodicity":..., "phase":..., "k":...}, ...}
        torsion_params = self.param_dict.get("TORSION", {})

        tors = self._get_torsions(topology)
        for i, j, k, l in tors:
            resA = atoms[i].residue.name
            resB = atoms[j].residue.name
            resC = atoms[k].residue.name
            resD = atoms[l].residue.name

            key1 = (resA, resB, resC, resD)
            tp = None
            if key1 in torsion_params:
                tp = torsion_params[key1]
            # elif ... (if you want symmetrical checks)

            if tp is not None:
                periodicity = tp["periodicity"]
                phase = tp["phase"]  # rad
                k_tor = tp["k"]  # kJ/mol

                torsion_force.addTorsion(
                    i,
                    j,
                    k,
                    l,
                    periodicity,
                    phase * unit.radian,
                    k_tor * unit.kilojoule_per_mole,
                )

        # -------------------
        # 5) Set nonbonded interaction
        # -------------------
        nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

        return system

    def _build_adjacency_list(self, topology: app.Topology) -> Dict[int, set]:
        adjacency = defaultdict(set)
        for bond in topology.bonds():
            a1, a2 = bond  # app.TopologyAtom
            i1 = a1.index
            i2 = a2.index
            adjacency[i1].add(i2)
            adjacency[i2].add(i1)
        return adjacency

    def _get_angles(self, topology: app.Topology) -> List[Tuple[int, int, int]]:
        adjacency = self._build_adjacency_list(topology)
        angles = []
        for j in adjacency:
            neighbors_j = list(adjacency[j])
            for i in neighbors_j:
                if i < j:
                    for k in neighbors_j:
                        if k > j and k != i:
                            angles.append((i, j, k))
        return angles

    def _get_torsions(self, topology: app.Topology) -> List[Tuple[int, int, int, int]]:
        adjacency = self._build_adjacency_list(topology)
        angles = self._get_angles(topology)  # (i, j, k)

        torsions = []
        for i, j, k in angles:
            # k の隣接を見る
            for l in adjacency[k]:
                if l > k and l not in (i, j):
                    torsions.append((i, j, k, l))
        return torsions
