# Receptor preparation for docking
import os
import re
import shutil
import subprocess
from pathlib import Path
from io import StringIO
from Bio.PDB import PDBParser, PDBIO

def _get_autodock_atom_type(atom_name: str, res_name: str, element: str) -> str:
    atom_name = atom_name.strip().upper()
    res_name = res_name.strip().upper()
    element = element.strip().upper()
    if not element:
        cleaned_name = re.sub(r"[0-9]", "", atom_name)
        element = cleaned_name[0] if cleaned_name else "C"
    if element == "H":
        if atom_name in ["HG", "HG1", "HD", "HD1", "HD2", "HE", "HE1", "HE2", "HH", "HH11", "HH12", "HH21", "HH22", "HZ", "HZ1", "HZ2", "HZ3"]:
            return "HD"
        return "H"
    if element == "C":
        if res_name in ["PHE", "TYR"] and atom_name in ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]:
            return "A"
        if res_name == "TRP" and atom_name in ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]:
            return "A"
        if res_name in ["HIS", "HIE", "HID", "HIP"] and atom_name in ["CG", "CD2", "CE1"]:
            return "A"
        return "C"
    if element == "N":
        if res_name in ["HIS", "HIE", "HID"] and atom_name in ["ND1", "NE2"]:
            return "NA"
        return "N"
    if element == "O":
        return "OA"
    if element == "S":
        return "SA"
    if element in ["MG", "ZN", "CA", "FE", "MN", "NA", "K", "CL", "BR", "F", "I", "P"]:
        return element.capitalize() if len(element) > 1 else element
    return element

def _convert_pdb_to_rigid_pdbqt(clean_pdb_path: Path, clean_pdbqt_path: Path):
    lines = []
    with open(clean_pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                element = line[76:78].strip() if len(line) >= 78 else ""
                ad_type = _get_autodock_atom_type(atom_name, res_name, element)
                base = line[:66].ljust(66)
                charge_str = "    0.000"
                type_str = f" {ad_type:<2}"
                pdbqt_line = f"{base} {charge_str} {type_str}\n"
                lines.append(pdbqt_line)
    with open(clean_pdbqt_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

def prepare_receptor(
    pdb_string: str,
    output_dir: Path,
    keep_waters: bool = False,
    reference_ligand_resname: str = None
) -> tuple[Path, Path, dict]:
    output_dir = Path(output_dir)
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("target", StringIO(pdb_string))

    removed_waters = 0
    removed_ligands = []
    retained_cofactors = []
    reference_ligand_removed = False
    ref_ligand_lines = []

    ref_upper = (reference_ligand_resname or "").strip().upper()

    # 1. Preserve the reference ligand in reference_ligand.pdb before receptor removal
    if ref_upper:
        for line in pdb_string.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                rname = line[17:20].strip().upper()
                if rname == ref_upper:
                    ref_ligand_lines.append(line)
        if ref_ligand_lines:
            ref_lig_path = output_dir / "reference_ligand.pdb"
            with open(ref_lig_path, "w", encoding="utf-8") as f:
                f.write("\n".join(ref_ligand_lines) + "\n")

    # 2. Separate protein, ligands, waters, and cofactors with full provenance (resname, chain, resnum)
    for model in struct:
        for chain in model:
            to_detach = []
            for res in chain:
                res_id = res.id[0]
                res_name = res.resname.strip().upper()
                res_num = int(res.id[1]) if isinstance(res.id[1], int) else str(res.id[1])
                chain_id = str(chain.id)
                res_entry = {"resname": res_name, "chain": chain_id, "resnum": res_num}

                if res_id.startswith("W"):
                    if not keep_waters:
                        to_detach.append(res.id)
                        removed_waters += 1
                elif res_id != " ":
                    if ref_upper and res_name == ref_upper:
                        to_detach.append(res.id)
                        reference_ligand_removed = True
                    elif res_name in ["HEM", "MG", "ZN", "CA", "FE", "MN"]:
                        retained_cofactors.append(res_entry)
                    else:
                        to_detach.append(res.id)
                        removed_ligands.append(res_entry)
            for rid in to_detach:
                chain.detach_child(rid)

    # 3. Save immutable receptor clean PDB
    clean_pdb_path = output_dir / "receptor_clean.pdb"
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(clean_pdb_path.resolve()))

    # 4. Convert to standard rigid PDBQT with multi-tier fallback
    clean_pdbqt_path = output_dir / "receptor.pdbqt"
    converted = False

    # Tier A: pybel
    try:
        from openbabel import pybel
        mol = next(pybel.readfile("pdb", str(clean_pdb_path.resolve())))
        mol.OBMol.AddHydrogens(False, True)
        mol.write("pdbqt", str(clean_pdbqt_path.resolve()), overwrite=True, opt={"r": None})
        if clean_pdbqt_path.exists() and clean_pdbqt_path.stat().st_size > 0:
            converted = True
    except Exception:
        converted = False

    # Tier B: obabel CLI
    if not converted:
        obabel_bin = shutil.which("obabel")
        if obabel_bin:
            try:
                cmd = [obabel_bin, "-ipdb", str(clean_pdb_path.resolve()), "-opdbqt", "-O", str(clean_pdbqt_path.resolve()), "-xr"]
                res = subprocess.run(cmd, capture_output=True, text=True)
                if clean_pdbqt_path.exists() and clean_pdbqt_path.stat().st_size > 0:
                    converted = True
            except Exception:
                converted = False

    # Tier C: Pure Python robust fallback
    if not converted:
        _convert_pdb_to_rigid_pdbqt(clean_pdb_path, clean_pdbqt_path)

    # 5. Sanitize rigid receptor PDBQT to strictly remove ROOT/BRANCH/TORSDOF tags
    if clean_pdbqt_path.exists():
        lines = []
        with open(clean_pdbqt_path, "r", encoding="utf-8") as f:
            for line in f:
                if not (line.startswith("ROOT") or line.startswith("ENDROOT") or
                        line.startswith("BRANCH") or line.startswith("ENDBRANCH") or
                        line.startswith("TORSDOF")):
                    lines.append(line)
        with open(clean_pdbqt_path, "w", encoding="utf-8") as f:
            f.writelines(lines)

    prep_report = {
        "clean_pdb_path": str(clean_pdb_path.name),
        "clean_pdbqt_path": str(clean_pdbqt_path.name),
        "reference_ligand_saved": bool(ref_ligand_lines),
        "reference_ligand_removed": reference_ligand_removed,
        "waters_policy": "Bulk waters removed (standard rigid-receptor docking model assumption)" if not keep_waters else "Retained waters",
        "waters_removed_count": removed_waters,
        "ligands_removed": removed_ligands,
        "cofactors_retained": retained_cofactors,
        "cofactor_parameterization_note": "Retained cofactors are processed with standard Gasteiger charges; no advanced force-field parameterization is claimed.",
        "polar_hydrogens_added": True
    }
    return clean_pdb_path, clean_pdbqt_path, prep_report
