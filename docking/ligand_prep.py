# Ligand preparation using RDKit and OpenBabel
import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel

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
    ob_mol = next(pybel.readfile("sdf", str(sdf_path.resolve())))
    ob_mol.write("pdbqt", str(pdbqt_path.resolve()), overwrite=True)

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
