# Protein-Ligand Docking Scientific Methods & Architecture

## 1. Overview & Scientific Purpose
This module implements rigid-receptor, flexible-ligand, non-covalent protein-ligand docking in accordance with the project docking brief and literature citations (Ferreira et al., 2015; Whalen et al., Lippincott Pharmacology).

Docking is treated strictly as an **in-silico hypothesis-generation and pose-prioritisation workflow**. Docking scores rank candidate conformations within a uniform protocol; they do not constitute measured binding affinity, clinical efficacy, or drug dosing approval.

## 2. Pinned Software Stack & Tools
- **Docking Engine**: AutoDock Vina v1.2.7 (pinned local binary under tools/vina/vina.exe)
- **Ligand 3D Conformer & Protonation**: RDKit & Meeko
- **Receptor & Ligand PDBQT Conversion**: OpenBabel / Meeko
- **Coordinate System**: Biological Assembly 1 (RCSB PDB) via existing verified pipeline

## 3. Strict Scientific Gates
- **Gate A (Redocking Validation)**: The co-crystallised small molecule from a validated target (e.g. 1M17 with ligand AQ4 Erlotinib) is removed, prepared, and re-docked. An RMSD threshold of <= 2.0 A against crystal coordinates is required.
- **Gate B (Physical Inspection)**: Poses are evaluated for steric clashes, pocket occupancy, and pose clustering across independent seeds.
- **Gate C (Matched Healthy vs Mutant Comparison)**: Requires identical ligand state, identical binding site grid, and identical exhaustiveness. Delta score is reported as a model-specific ranking difference.

## 4. Run Manifest & Reproducibility
Every run writes an immutable manifest to docking_runs/<run_id>/manifest.json detailing all run parameters, seeds, scores, and coordinates.
