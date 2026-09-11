# Ligand preparation using RDKit and OpenBabel
import os
import shutil
import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

def prepare_ligand_from_smiles(smiles: str, ligand_name: str, output_dir: Path, ph: float = 7.4) -> tuple[Path, Path, dict]:
    output_dir = Path(output_dir)
    smiles = (smiles or "").strip()
    if not smiles:
        raise ValueError("Empty SMILES string provided.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: '{smiles}'")

    mol = Chem.AddHs(mol)
    embed_code = AllChem.EmbedMolecule(mol, randomSeed=42)
    if embed_code != 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass

    sdf_path = output_dir / f"{ligand_name}_3d.sdf"
    writer = Chem.SDWriter(str(sdf_path.resolve()))
    writer.write(mol)
    writer.close()

    pdbqt_path = output_dir / f"{ligand_name}.pdbqt"
    converted = False

    # Tier A: pybel readfile sdf
    try:
        from openbabel import pybel
        ob_mol = next(pybel.readfile("sdf", str(sdf_path.resolve())))
        ob_mol.write("pdbqt", str(pdbqt_path.resolve()), overwrite=True)
        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            converted = True
    except Exception:
        converted = False

    # Tier B: pybel readstring mol block
    if not converted:
        try:
            from openbabel import pybel
            mb = Chem.MolToMolBlock(mol)
            ob_mol = pybel.readstring("mol", mb)
            ob_mol.write("pdbqt", str(pdbqt_path.resolve()), overwrite=True)
            if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
                converted = True
        except Exception:
            converted = False

    # Tier C: obabel CLI
    if not converted:
        obabel_bin = shutil.which("obabel")
        if obabel_bin:
            try:
                cmd = [obabel_bin, "-isdf", str(sdf_path.resolve()), "-opdbqt", "-O", str(pdbqt_path.resolve())]
                subprocess.run(cmd, capture_output=True, text=True)
                if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
                    converted = True
            except Exception:
                converted = False

    # Tier D: pybel readstring pdb
    if not converted:
        try:
            from openbabel import pybel
            pb = Chem.MolToPDBBlock(mol)
            ob_mol = pybel.readstring("pdb", pb)
            ob_mol.write("pdbqt", str(pdbqt_path.resolve()), overwrite=True)
            if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
                converted = True
        except Exception:
            converted = False

    if not converted or not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
        raise RuntimeError(f"Failed to generate PDBQT for ligand '{ligand_name}'.")

    num_rotatable_bonds = AllChem.CalcNumRotatableBonds(mol)

    report = {
        "ligand_name": ligand_name,
        "smiles": smiles,
        "sdf_file": str(sdf_path.name),
        "pdbqt_file": str(pdbqt_path.name),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "total_atom_count": mol.GetNumAtoms(),
        "rotatable_bonds": num_rotatable_bonds,
        "ph_assumed": ph
    }
    return sdf_path, pdbqt_path, report
