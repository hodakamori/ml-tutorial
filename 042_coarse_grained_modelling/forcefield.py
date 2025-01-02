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

        # -- Forces
        bond_force = openmm.HarmonicBondForce()
        system.addForce(bond_force)

        angle_force = openmm.HarmonicAngleForce()
        system.addForce(angle_force)

        torsion_force = openmm.PeriodicTorsionForce()
        system.addForce(torsion_force)

        nonbonded_force = openmm.NonbondedForce()
        system.addForce(nonbonded_force)

        # 1) nonbond params: mass, LJ, charge
        nonbonded_params = self.param_dict.get("NONBONDED", {})
        for i, atom in enumerate(atoms):
            res_label = atom.residue.name
            nb_data = nonbonded_params[res_label]
            mass_value = nb_data["mass"]  # amu
            sigma_nm = nb_data["sigma"]  # nm
            epsilon_kj = nb_data["epsilon"]  # kJ/mol
            charge_e = nb_data["charge"]  # e

            # mass
            system.addParticle(mass_value * unit.amu)
            # nonbond (LJ + coulomb)
            nonbonded_force.addParticle(
                charge_e,
                sigma_nm * unit.nanometer,
                epsilon_kj * unit.kilojoule_per_mole,
            )
            adjacency = defaultdict(set)

            # bond一覧を取得
            bond_list = list(topology.bonds())
            # adjacency list
            for bond in bond_list:
                a1, a2 = bond
                i1 = a1.index
                i2 = a2.index
                adjacency[i1].add(i2)
                adjacency[i2].add(i1)

            # 1-2: 完全除外
            for atom1, atom2 in bond_list:
                i1, i2 = atom1.index, atom2.index
                # addException( p1, p2, chargeProd, sigma, epsilon, replace=False )
                # 完全除外 => chargeProd=0.0, epsilon=0.0
                nonbonded_force.addException(i1, i2, 0.0, 1.0, 0.0, replace=True)

            angle_list = []
            for j in adjacency:
                neigh = list(adjacency[j])
                for i in neigh:
                    if i < j:
                        for k in neigh:
                            if k > j and k != i:
                                angle_list.append((i, j, k))
                                # 1-3 => exclude
                                nonbonded_force.addException(
                                    i, k, 0.0, 1.0, 0.0, replace=True
                                )
        nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        nonbonded_force.setCutoffDistance(1.1)

        # 2) BondForce
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
            else:
                raise ValueError(f"Unknown bond: {label1}-{label2}")

            r0 = bp["r0"]  # nm
            k_val = bp["k"]  # kJ/(mol·nm^2)
            bond_force.addBond(
                i1,
                i2,
                r0 * unit.nanometer,
                k_val * unit.kilojoule_per_mole / (unit.nanometer**2),
            )
            # system.addConstraint(i1, i2, r0 * unit.nanometer)

        angle_params = self.param_dict.get("ANGLE", {})
        angles = self._get_angles(topology)  # (i,j,k)

        for i, j, k in angles:
            resA = atoms[i].residue.name
            resB = atoms[j].residue.name
            resC = atoms[k].residue.name

            key1 = (resA, resB, resC)
            if key1 in angle_params:
                ap = angle_params[key1]
                theta0_rad = ap["theta0"]
                k_rad = ap["k"]

                angle_force.addAngle(
                    i,
                    j,
                    k,
                    theta0_rad * unit.radian,
                    k_rad * unit.kilojoule_per_mole / (unit.radian**2),
                )

        torsion_params = self.param_dict.get("TORSION", {})
        tors_list = self._get_torsions(topology)

        for i, j, k, l in tors_list:
            resA = atoms[i].residue.name
            resB = atoms[j].residue.name
            resC = atoms[k].residue.name
            resD = atoms[l].residue.name

            key1 = (resA, resB, resC, resD)
            if key1 in torsion_params:
                tp = torsion_params[key1]
                periodicity = tp["periodicity"]
                phase_rad = tp["phase"]  # already in rad
                k_tor = tp["k"]  # kJ/mol

                torsion_force.addTorsion(
                    i,
                    j,
                    k,
                    l,
                    periodicity,
                    phase_rad * unit.radian,
                    k_tor * unit.kilojoule_per_mole,
                )

        # 6) set PBC
        box_dims = topology.getUnitCellDimensions()
        if box_dims is not None:
            lx, ly, lz = box_dims
            system.setDefaultPeriodicBoxVectors(
                openmm.Vec3(lx, 0, 0) * unit.nanometer,
                openmm.Vec3(0, ly, 0) * unit.nanometer,
                openmm.Vec3(0, 0, lz) * unit.nanometer,
            )
            system.usesPeriodicBoundaryConditions = True

        return system

    def _build_adjacency_list(self, topology: app.Topology) -> Dict[int, set]:
        adjacency = defaultdict(set)
        for bond in topology.bonds():
            a1, a2 = bond
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
            for l in adjacency[k]:
                if l > k and l not in (i, j):
                    torsions.append((i, j, k, l))
        return torsions
