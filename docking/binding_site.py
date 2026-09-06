# Binding site definition and grid calculations
import re
import numpy as np
from io import StringIO
from Bio.PDB import PDBParser
from docking.models import GridBox

def extract_reference_ligand(pdb_string: str, ligand_resname: str) -> tuple[np.ndarray, str]:
    """Extract coordinates and PDB lines for a specific co-crystallised reference ligand."""
    ligand_resname = (ligand_resname or "").strip().upper()
    if not ligand_resname:
        raise ValueError("No reference ligand code provided.")

    coords = []
    ligand_lines = []

    for line in pdb_string.splitlines():
        if line.startswith("HETATM") or line.startswith("ATOM"):
            rname = line[17:20].strip().upper()
            if rname == ligand_resname:
                ligand_lines.append(line)
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass

    if not coords:
        raise ValueError(f"Reference ligand '{ligand_resname}' was not found in PDB coordinates.")

    coords_arr = np.array(coords, dtype=np.float32)
    return coords_arr, "\n".join(ligand_lines)

def define_grid_from_ligand(pdb_string: str, ligand_resname: str, padding: float = 8.0) -> tuple[GridBox, np.ndarray, str]:
    """Define a finite 3D docking search box centered on the crystal reference ligand."""
    coords_arr, ligand_pdb_text = extract_reference_ligand(pdb_string, ligand_resname)
    center = coords_arr.mean(axis=0)

    ranges = coords_arr.max(axis=0) - coords_arr.min(axis=0)
    sizes = np.maximum(ranges + 2 * padding, 15.0)

    grid = GridBox(
        center_x=round(float(center[0]), 3),
        center_y=round(float(center[1]), 3),
        center_z=round(float(center[2]), 3),
        size_x=round(float(sizes[0]), 1),
        size_y=round(float(sizes[1]), 1),
        size_z=round(float(sizes[2]), 1)
    )
    return grid, coords_arr, ligand_pdb_text

def define_grid_from_residues(pdb_string: str, residue_specs: list[str], padding: float = 10.0) -> GridBox:
    """
    Define a finite 3D docking search box based on a validated list of pocket residues.
    Example residue_specs: ['A:790', 'A:858'] or ['790', '858'].
    Prohibits origin-centred fallbacks.
    """
    if not residue_specs:
        raise ValueError("No binding-site residues provided. Docking requires a validated pocket.")

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("target", StringIO(pdb_string))

    coords = []
    target_set = set(s.strip().upper() for s in residue_specs if s.strip())

    for model in struct:
        for chain in model:
            for res in chain:
                res_num = str(res.id[1])
                full_spec = f"{chain.id}:{res_num}".upper()
                if full_spec in target_set or res_num in target_set:
                    for atom in res.get_atoms():
                        coords.append(atom.get_coord())

    if not coords:
        raise ValueError(
            f"None of the specified residues {residue_specs} were found in the receptor coordinates."
        )

    coords_arr = np.array(coords, dtype=np.float32)
    center = coords_arr.mean(axis=0)
    ranges = coords_arr.max(axis=0) - coords_arr.min(axis=0)
    sizes = np.maximum(ranges + 2 * padding, 18.0)

    return GridBox(
        center_x=round(float(center[0]), 3),
        center_y=round(float(center[1]), 3),
        center_z=round(float(center[2]), 3),
        size_x=round(float(sizes[0]), 1),
        size_y=round(float(sizes[1]), 1),
        size_z=round(float(sizes[2]), 1)
    )
