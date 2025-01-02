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
        """
        全体の流れ:
         - system, forces を作成
         - nonbond, bond, angle, torsion の設定
         - 1-2,1-3 exclusionを高速に行う
         - param_dict["ANGLE"] は (already rad, kJ/(mol·rad^2)) など
         - param_dict["TORSION"] は 単一 or 複数辞書 (list)
        """
        system = openmm.System()
        atoms = list(topology.atoms())

        # (1) Forces
        bond_force = openmm.HarmonicBondForce()
        angle_force = openmm.HarmonicAngleForce()
        torsion_force = openmm.PeriodicTorsionForce()
        nonbonded_force = openmm.NonbondedForce()

        system.addForce(bond_force)
        system.addForce(angle_force)
        system.addForce(torsion_force)
        system.addForce(nonbonded_force)

        # (2) Nonbonded particles
        nonbonded_params = self.param_dict.get("NONBONDED", {})
        for i, atom in enumerate(atoms):
            res_label = atom.residue.name
            nb_data = nonbonded_params[res_label]
            mass_value = nb_data["mass"]  # amu
            sigma_nm = nb_data["sigma"]  # nm
            epsilon_kj = nb_data["epsilon"]  # kJ/mol
            charge_e = nb_data["charge"]  # e

            system.addParticle(mass_value * unit.amu)
            nonbonded_force.addParticle(
                charge_e,
                sigma_nm * unit.nanometer,
                epsilon_kj * unit.kilojoule_per_mole,
            )

        # (3) build adjacency for 1-2,1-3 exclusion
        adjacency = defaultdict(set)
        bond_list = list(topology.bonds())
        for bond in bond_list:
            a1, a2 = bond
            i1, i2 = a1.index, a2.index
            adjacency[i1].add(i2)
            adjacency[i2].add(i1)

        # 1-2, 1-3 exclusion in a single pass
        #  => for i in adjacency:
        #       for j in adjacency[i]: (j> i) => 1-2
        #         for k in adjacency[i]: if k> j => 1-3
        for i in adjacency:
            neigh_i = sorted(adjacency[i])  # ensure i< j < k if you like
            for j in neigh_i:
                if j > i:
                    # 1-2 => exclude
                    nonbonded_force.addException(i, j, 0.0, 1.0, 0.0, replace=True)
            # 1-3
            #   if i-j and i-k => j,k in adjacency[i]
            #   exclude i,k
            #   we do j<k for k>j
            #   also k != i => trivially, but i is fixed
            for idx_a in range(len(neigh_i)):
                j = neigh_i[idx_a]
                if j <= i:
                    continue
                for idx_b in range(idx_a + 1, len(neigh_i)):
                    k = neigh_i[idx_b]
                    if k > j:
                        # 1-3 => exclude
                        nonbonded_force.addException(i, k, 0.0, 1.0, 0.0, replace=True)

        nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        nonbonded_force.setCutoffDistance(1.1 * unit.nanometer)

        # (4) BondForce
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
                raise ValueError(f"Unknown bond param for: {label1}-{label2}")

            r0 = bp["r0"]  # nm
            k_val = bp["k"]  # kJ/(mol·nm^2)
            bond_force.addBond(
                i1,
                i2,
                r0 * unit.nanometer,
                k_val * unit.kilojoule_per_mole / (unit.nanometer**2),
            )

        # (5) AngleForce (angles in rad)
        angle_params = self.param_dict.get("ANGLE", {})
        angles = self._get_angles(topology)
        for i, j, k in angles:
            resA = atoms[i].residue.name
            resB = atoms[j].residue.name
            resC = atoms[k].residue.name
            key1 = (resA, resB, resC)
            if key1 in angle_params:
                ap = angle_params[key1]
                theta0_rad = ap["theta0"]  # rad
                k_rad = ap["k"]  # kJ/(mol·rad^2)
                angle_force.addAngle(
                    i,
                    j,
                    k,
                    theta0_rad * unit.radian,
                    k_rad * unit.kilojoule_per_mole / (unit.radian**2),
                )

        # (6) TorsionForce => single or multiple dict
        torsion_params = self.param_dict.get("TORSION", {})
        tors_list = self._get_torsions(topology)
        for i, j, k, l in tors_list:
            resA = atoms[i].residue.name
            resB = atoms[j].residue.name
            resC = atoms[k].residue.name
            resD = atoms[l].residue.name
            key1 = (resA, resB, resC, resD)
            if key1 in torsion_params:
                val = torsion_params[key1]
                if isinstance(val, dict):
                    # single dict
                    periodicity = val["periodicity"]
                    phase_rad = val["phase"]  # rad
                    k_tor = val["k"]  # kJ/mol
                    torsion_force.addTorsion(
                        i,
                        j,
                        k,
                        l,
                        periodicity,
                        phase_rad * unit.radian,
                        k_tor * unit.kilojoule_per_mole,
                    )
                elif isinstance(val, list):
                    # multiple
                    for tp in val:
                        periodicity = tp["periodicity"]
                        phase_rad = tp["phase"]
                        k_tor = tp["k"]
                        torsion_force.addTorsion(
                            i,
                            j,
                            k,
                            l,
                            periodicity,
                            phase_rad * unit.radian,
                            k_tor * unit.kilojoule_per_mole,
                        )
                else:
                    raise ValueError(
                        f"TORSION param for {key1} must be dict or list of dict"
                    )

        # (7) PBC
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

    # adjacency builder
    def _build_adjacency_list(self, topology: app.Topology) -> Dict[int, set]:
        adjacency = defaultdict(set)
        for bond in topology.bonds():
            a1, a2 = bond
            i1, i2 = a1.index, a2.index
            adjacency[i1].add(i2)
            adjacency[i2].add(i1)
        return adjacency

    def _get_angles(self, topology: app.Topology) -> List[Tuple[int, int, int]]:
        adjacency = self._build_adjacency_list(topology)
        angles = []
        for j in adjacency:
            neigh = sorted(adjacency[j])
            for idx_a in range(len(neigh)):
                i = neigh[idx_a]
                if i < j:
                    for idx_b in range(idx_a + 1, len(neigh)):
                        k = neigh[idx_b]
                        if k > j and k != i:
                            angles.append((i, j, k))
        return angles

    def _get_torsions(self, topology: app.Topology) -> List[Tuple[int, int, int, int]]:
        adjacency = self._build_adjacency_list(topology)
        # build angles first
        angles = []
        for j in adjacency:
            neigh = sorted(adjacency[j])
            for idx_a in range(len(neigh)):
                i = neigh[idx_a]
                if i < j:
                    for idx_b in range(idx_a + 1, len(neigh)):
                        k = neigh[idx_b]
                        if k > j and k != i:
                            angles.append((i, j, k))

        torsions = []
        for i, j, k in angles:
            neigh_k = adjacency[k]
            for l in neigh_k:
                if l > k and l not in (i, j):
                    torsions.append((i, j, k, l))
        return torsions
