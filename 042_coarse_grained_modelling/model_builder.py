import math
import random
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import openmm
import openmm.app as app
from openmm import Vec3

sys.setrecursionlimit(100_000)


class SmilesNode:
    def __init__(self, label: Optional[str]):
        self.label: Optional[str] = label
        self.children: List[SmilesNode] = []
        self.next_sibling: Optional[SmilesNode] = None
        self.atom_index: Optional[int] = None

    def __repr__(self):
        return self._repr_recursive(0)

    def _repr_recursive(self, level: int) -> str:
        indent = "  " * level
        lbl = self.label if self.label else "(Group)"
        idx = f"(idx={self.atom_index})" if self.atom_index is not None else ""
        s = f"{indent}{lbl}{idx}\n"
        for c in self.children:
            s += c._repr_recursive(level + 1)
        if self.next_sibling:
            s += self.next_sibling._repr_recursive(level)
        return s


_token_pattern = re.compile(
    r"""
    (?P<monomer>\[[^\]]*\]) |  # [EO], [OH]
    (?P<number>\d+)         |  # 3, 10 (rep count)
    (?P<lparen>\()          |  # (
    (?P<rparen>\))             # )
    """,
    re.VERBOSE,
)


def tokenize_smileslike(sequence: str) -> List[str]:
    tokens = []
    for m in _token_pattern.finditer(sequence):
        if m.lastgroup == "monomer":
            tokens.append(m.group("monomer"))
        elif m.lastgroup == "number":
            tokens.append(m.group("number"))
        elif m.lastgroup == "lparen":
            tokens.append("(")
        elif m.lastgroup == "rparen":
            tokens.append(")")
    return tokens


class TokenStream:
    def __init__(self, tokens: List[str]):
        self.tokens = tokens
        self.index = 0

    def peek(self) -> Optional[str]:
        if self.index < len(self.tokens):
            return self.tokens[self.index]
        return None

    def get(self) -> Optional[str]:
        if self.index < len(self.tokens):
            t = self.tokens[self.index]
            self.index += 1
            return t
        return None


def parse_smileslike(sequence: str) -> SmilesNode:
    tokens = tokenize_smileslike(sequence)
    stream = TokenStream(tokens)
    root = SmilesNode(label=None)
    _parse_group(stream, root)

    if stream.peek() is not None:
        raise ValueError("Unexpected token after parse.")
    return root


def _parse_group(stream: TokenStream, parent: SmilesNode):
    current_sibling: Optional[SmilesNode] = None
    while True:
        token = stream.peek()
        if token is None or token == ")":
            break

        if token == "(":
            stream.get()  # consume '('
            branch_root = SmilesNode(None)
            _parse_group(stream, branch_root)
            closing = stream.get()  # consume ')'
            if closing != ")":
                raise ValueError("Missing closing parenthesis.")

            rep_count = 1
            nxt = stream.peek()
            if nxt and nxt.isdigit():
                rep_count = int(stream.get())

            for _ in range(rep_count):
                clone = _clone_subtree(branch_root)
                if current_sibling is None:
                    parent.children.append(clone)
                else:
                    current_sibling.next_sibling = clone
                current_sibling = _get_last_sibling(clone)

        else:
            monomer = stream.get()  # 例: "[EO]"
            name = monomer.strip("[]")  # "EO"
            rep_count = 1
            nxt = stream.peek()
            if nxt and nxt.isdigit():
                rep_count = int(stream.get())

            for _ in range(rep_count):
                node = SmilesNode(name)
                if current_sibling is None:
                    parent.children.append(node)
                else:
                    current_sibling.next_sibling = node
                current_sibling = node


def _clone_subtree(node: SmilesNode) -> SmilesNode:
    new_node = SmilesNode(node.label)
    new_node.children = [_clone_subtree(c) for c in node.children]
    if node.next_sibling:
        new_node.next_sibling = _clone_subtree(node.next_sibling)
    return new_node


def _get_last_sibling(node: SmilesNode) -> SmilesNode:
    cur = node
    while cur.next_sibling:
        cur = cur.next_sibling
    return cur


@dataclass
class CGModelBuilder:
    smiles_sequence: str
    bond_length: float
    min_dist: float
    max_tries: int
    tolerance: float = 2.0
    filetype: str = "pdb"
    verbose: bool = False
    packmol_path: str = "packmol"

    def create_molecule(self) -> Tuple[app.Topology, List[openmm.Vec3]]:
        root = parse_smileslike(self.smiles_sequence)

        topology = app.Topology()
        chain = topology.addChain()
        atom_list = []
        bonds = []

        def traverse_create_atoms(node: SmilesNode, parent_idx: Optional[int]):
            if node.label is not None:
                residue = topology.addResidue(node.label, chain)
                elem = app.element.carbon
                new_atom = topology.addAtom(node.label, elem, residue)
                my_idx = len(atom_list)
                atom_list.append(new_atom)
                node.atom_index = my_idx

                if parent_idx is not None:
                    bonds.append((parent_idx, my_idx))
                current_parent = my_idx
            else:
                current_parent = parent_idx

            for child in node.children:
                traverse_create_atoms(child, current_parent)

            if node.next_sibling:
                traverse_create_atoms(node.next_sibling, current_parent)

        for c in root.children:
            traverse_create_atoms(c, None)

        for i1, i2 in bonds:
            topology.addBond(atom_list[i1], atom_list[i2])

        positions = [openmm.Vec3(0, 0, 0) for _ in range(len(atom_list))]

        def traverse_assign_positions(node: SmilesNode, parent_idx: Optional[int]):
            if node.label is not None:
                my_idx = node.atom_index
                if parent_idx is None:
                    positions[my_idx] = openmm.Vec3(0.0, 0.0, 0.0)
                else:
                    parent_pos = positions[parent_idx]
                    candidate = self._find_random_position(parent_pos, positions)
                    positions[my_idx] = candidate

                current_parent = my_idx
            else:
                current_parent = parent_idx

            for child in node.children:
                traverse_assign_positions(child, current_parent)

            if node.next_sibling:
                traverse_assign_positions(node.next_sibling, current_parent)

        for c in root.children:
            traverse_assign_positions(c, None)

        return topology, positions

    def _find_random_position(
        self, parent_pos: openmm.Vec3, positions: List[openmm.Vec3]
    ) -> openmm.Vec3:
        for _ in range(self.max_tries):
            direction = self._random_3d_unit_vector()  # unit vector
            candidate = parent_pos + direction * self.bond_length
            if not self._too_close(candidate, positions):
                return candidate
        raise RuntimeError("Could not find valid random position (max_tries exceeded).")

    def _random_3d_unit_vector(self) -> openmm.Vec3:
        u = random.random()  # 0..1
        v = random.random()  # 0..1
        theta = math.acos(1.0 - 2.0 * u)  # [0, π]
        phi = 2.0 * math.pi * v  # [0, 2π)

        sin_t = math.sin(theta)
        x = sin_t * math.cos(phi)
        y = sin_t * math.sin(phi)
        z = math.cos(theta)
        return openmm.Vec3(x, y, z)

    def _too_close(self, candidate: openmm.Vec3, positions: List[openmm.Vec3]) -> bool:
        for p in positions:
            dx = candidate.x - p.x
            dy = candidate.y - p.y
            dz = candidate.z - p.z
            dist_sq = dx * dx + dy * dy + dz * dz
            if dist_sq < (self.min_dist * self.min_dist):
                return True
        return False

    def create_packed_model(
        self,
        structures: List[Dict],
        box_size: Tuple[float, float, float],
        output_pdb: str = "mixture.pdb",
    ) -> Tuple[app.Topology, List[Vec3]]:
        packmol_input = self._build_packmol_input(structures, output_pdb)

        inp_file = "packmol_auto.inp"
        with open(inp_file, "w") as f:
            f.write(packmol_input)

        cmd = [self.packmol_path]
        with open("packmol.log", "w") as logf:
            proc = subprocess.Popen(
                cmd,
                stdin=open(inp_file, "r"),
                stdout=logf if not self.verbose else None,
                stderr=logf if not self.verbose else None,
            )
            proc.wait()

        if proc.returncode != 0:
            raise RuntimeError("PACKMOL failed. See packmol.log for details.")

        pdb = app.PDBFile(output_pdb)
        topology = pdb.topology
        positions_angs = pdb.positions  # unit: Å(A)
        positions_nm = []
        for pos in positions_angs:
            # Å -> nm
            x_nm = pos.x
            y_nm = pos.y
            z_nm = pos.z
            positions_nm.append(Vec3(x_nm, y_nm, z_nm))
        lx_nm = box_size[0]
        ly_nm = box_size[1]
        lz_nm = box_size[2]
        topology.setUnitCellDimensions((lx_nm, ly_nm, lz_nm))
        return topology, positions_nm

    def _build_packmol_input(self, structures: List[Dict], output_pdb: str) -> str:
        lines = []
        lines.append(f"tolerance {self.tolerance}")
        lines.append(f"filetype {self.filetype}")
        lines.append(f"output {output_pdb}")
        lines.append("")

        for st in structures:
            pdb_file = st["pdb"]
            num_mol = st["number"]
            box = st["box"]  # (x1,y1,z1, x2,y2,z2) in Å
            lines.append(f"structure {pdb_file}")
            lines.append(f"   number {num_mol}")
            lines.append(
                f"   inside box {box[0]} {box[1]} {box[2]} {box[3]} {box[4]} {box[5]}"
            )
            lines.append("end structure")
            lines.append("")

        return "\n".join(lines)

    def replicate_model(
        self,
        orig_topology: app.Topology,
        orig_positions: List[Vec3],
        nx: int = 2,
        ny: int = 2,
        nz: int = 2,
    ) -> Tuple[app.Topology, List[Vec3]]:
        orig_atoms = list(orig_topology.atoms())
        n_orig_atoms = len(orig_atoms)
        if n_orig_atoms != len(orig_positions):
            raise ValueError("Topology and positions do not match.")
        xs = [pos.x for pos in orig_positions]
        ys = [pos.y for pos in orig_positions]
        zs = [pos.z for pos in orig_positions]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        min_z, max_z = min(zs), max(zs)

        box_size_x = max_x - min_x
        box_size_y = max_y - min_y
        box_size_z = max_z - min_z

        new_topology = app.Topology()
        atom_mapping = {}

        orig_chains = list(orig_topology.chains())
        n_replicates = nx * ny * nz

        new_positions_count = n_replicates * n_orig_atoms
        new_positions = [None] * new_positions_count

        replicate_id = 0

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    shift_x = ix * box_size_x
                    shift_y = iy * box_size_y
                    shift_z = iz * box_size_z

                    for orig_chain in orig_chains:
                        new_chain = new_topology.addChain(id=orig_chain.id)

                        for orig_res in orig_chain.residues():
                            new_res = new_topology.addResidue(
                                orig_res.name, new_chain, id=orig_res.id
                            )

                            for orig_atom in orig_res.atoms():
                                elem = orig_atom.element
                                if elem is None:
                                    elem = app.element.get_by_symbol("C")
                                new_atom = new_topology.addAtom(
                                    orig_atom.name, elem, new_res
                                )

                                old_idx = orig_atom.index
                                new_pos_index = replicate_id * n_orig_atoms + old_idx

                                orig_pos = orig_positions[old_idx]
                                shifted_pos = Vec3(
                                    orig_pos.x + shift_x,
                                    orig_pos.y + shift_y,
                                    orig_pos.z + shift_z,
                                )
                                new_positions[new_pos_index] = shifted_pos
                                atom_mapping[(replicate_id, old_idx)] = new_atom

                    replicate_id += 1

        orig_bonds = list(orig_topology.bonds())
        for rep_id in range(n_replicates):
            for atom1, atom2 in orig_bonds:
                old1 = atom1.index
                old2 = atom2.index
                new_atom1 = atom_mapping[(rep_id, old1)]
                new_atom2 = atom_mapping[(rep_id, old2)]
                new_topology.addBond(new_atom1, new_atom2)
        replicate_box_x = nx * box_size_x
        replicate_box_y = ny * box_size_y
        replicate_box_z = nz * box_size_z
        new_topology.setUnitCellDimensions(
            (replicate_box_x, replicate_box_y, replicate_box_z)
        )

        return new_topology, new_positions
