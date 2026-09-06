# Receptor preparation for docking
import os
from pathlib import Path
from io import StringIO
from Bio.PDB import PDBParser, PDBIO
from openbabel import pybel

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

    # 4. Convert to standard rigid PDBQT
    clean_pdbqt_path = output_dir / "receptor.pdbqt"
    mol = next(pybel.readfile("pdb", str(clean_pdb_path.resolve())))
    mol.OBMol.AddHydrogens(False, True)
    mol.write("pdbqt", str(clean_pdbqt_path.resolve()), overwrite=True, opt={"r": None})

    # 5. Sanitize rigid receptor PDBQT to strictly remove ROOT/BRANCH/TORSDOF tags
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
