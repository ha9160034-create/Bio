# Docking Code Repair Instructions

## Objective

Repair the protein-ligand docking implementation in \`C:\Users\hassan\Desktop\MY_VScode\` so that its results are reproducible, scientifically limited, and safe to interpret. Do not redesign the existing protein-comparison workflow. Do not make clinical, dosage, safety, \`Kd\`, \`Ki\`, \`IC50\`, \`EC50\`, efficacy, or treatment claims.

The project is **not accepted** until every P0 item and the acceptance tests below pass.

## P0 — required before any scientific result is shown

### 1. Replace the invalid redocking RMSD calculation

**Problem:** \`docking/validation.py\` compares coordinate arrays by order and silently truncates to the shorter array. This is not a valid ligand RMSD calculation. Atom order can differ between the crystal ligand and PDBQT output; a mismatched molecule could pass.

**Required implementation:**

1. Remove \`min_len\` truncation completely.
2. Create a representation for ligand atoms containing at minimum: atom name, element, coordinates, and a stable index from the source molecule.
3. For redocking, require that the reference ligand and the docked ligand represent the same chemical entity. Verify equal heavy-atom count and a documented atom mapping. If this cannot be established, set validation to \`validation unavailable\` and stop; never calculate a partial RMSD.
4. Superpose the receptor structures, or establish that both ligand coordinate sets are already in the same receptor coordinate frame. Then calculate RMSD over all mapped heavy atoms.
5. Save atom mapping method, atom count, RMSD, reference-ligand coordinates, and predicted-pose coordinates in \`manifest.json\`.
6. The redocking gate passes only when **the best valid predicted pose has heavy-atom RMSD <= 2.0 Å**. The result status must be \`validation passed\` or \`validation failed\`, not a green success badge for an arbitrary selected pose.

**Important:** the values \`rmsd_lb\` and \`rmsd_ub\` printed by Vina describe similarity between Vina poses. They are not RMSD to the experimental ligand.

### 2. Prohibit unvalidated default grids

**Problem:** \`docking/ui.py\` catches a reference-ligand error and docks around \`(0, 0, 0)\`. The no-reference path also uses a default origin-centred grid.

**Required implementation:**

- Remove both origin-centred fallback grids.
- Allow exactly two site-definition paths:
  1. a successfully extracted co-crystal/reference ligand; or
  2. a residue-defined pocket where the user supplies residues and the code computes/displays a finite grid.
- If neither is valid, stop before receptor/ligand preparation and show: \`Docking was not run: define a validated binding site.\`
- Record \`site_method\`, reference ligand code or residues, grid centre, grid dimensions, and padding in the manifest.

### 3. Execute all seeds and report reproducibility

**Problem:** the model defines \`[42, 101, 2024]\`, but UI code runs Vina only with seed \`42\`.

**Required implementation:**

1. Execute Vina once per requested seed.
2. Store one output PDBQT and one log file per seed.
3. Parse all poses with their source seed and rank.
4. Cluster poses across seeds by documented heavy-atom RMSD after atom mapping. Do not set \`cluster_id = rank\`.
5. Report each seed’s top score, the median score for each cross-seed cluster, and the number of seeds supporting the dominant cluster.
6. Flag the run if no pose cluster is reproduced by at least two independent seeds.
7. A single-seed run may be saved only as \`exploratory\`; it must never show \`validation passed\`.

### 4. Make manifests complete and mandatory

**Problem:** the existing run folders contain no \`manifest.json\`; the small manifest in the UI is insufficient.

**Required implementation:**

- Write \`manifest.initial.json\` before any external executable runs.
- Update it atomically to \`manifest.json\` after the run ends, whether the run succeeds or fails.
- Never delete prior manifests or overwrite another run directory.
- A run directory must contain at least:

\`\`\`text
manifest.json
receptor_source.pdb                 # immutable input coordinate file
receptor_clean.pdb
receptor.pdbqt
ligand_input.*
ligand_prepared.*
vina_seed_<seed>.log                # one per seed
out_seed_<seed>.pdbqt               # one per completed seed
reference_ligand.pdb                # when redocking
\`\`\`

The manifest must include:

\`\`\`json
{
  "schema_version": 1,
  "run_id": "...",
  "status": "completed | failed | validation_failed | exploratory",
  "created_at_utc": "...",
  "target": {
    "pdb_id": "...",
    "coordinate_source": "rcsb_biological_assembly | deposited | upload",
    "assembly_id": 1,
    "chains": ["..."],
    "mutation_label": "..."
  },
  "binding_site": {"method": "reference_ligand | residues"},
  "ligand": {"name": "...", "input_smiles": "...", "protonation_note": "..."},
  "receptor_preparation": {"waters": "...", "ligands_removed": [], "cofactors_retained": []},
  "engine": {"name": "AutoDock Vina", "version": "...", "seeds": [], "exhaustiveness": 8},
  "grid": {"center": ["x", "y", "z"], "size": ["x", "y", "z"]},
  "inputs": {"sha256": {}},
  "results": {"per_seed": [], "clusters": [], "redocking": {}},
  "warnings": []
}
\`\`\`

### 5. Correct result language and gate display

- Change the column title from \`Binding Energy\` / \`طاقة الارتباط\` to \`Docking score (Vina, kcal/mol)\` / \`درجة الإرساء (Vina، kcal/mol)\`.
- Display the supplied scientific disclaimer beside every result table.
- Do not use a red/green “effective” or “ineffective” outcome.
- Use only these statuses: \`validation passed\`, \`validation failed\`, \`exploratory—no reference validation\`, \`run failed\`, and \`validation unavailable\`.
- Do not use the old \`2.5 Å\` pass threshold. Use the defined \`<= 2.0 Å\` gate only after correct mapping/RMSD calculation.

## P1 — required for a reliable application

### 6. Use project-root-relative paths

**Problem:** \`VINA_PATH = os.path.abspath(os.path.join('tools', 'vina', 'vina.exe'))\` and \`BASE_RUNS_DIR = 'docking_runs'\` depend on the shell’s current working directory.

**Required implementation:**

\`\`\`python
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
VINA_PATH = PROJECT_ROOT / "tools" / "vina" / "vina.exe"
BASE_RUNS_DIR = PROJECT_ROOT / "docking_runs"
\`\`\`

Use \`Path\` consistently and create paths only below \`BASE_RUNS_DIR\`. Validate that \`VINA_PATH\` is an existing file before the UI enables execution.

### 7. Eliminate the circular import from \`b.py\`

**Problem:** \`docking/ui.py\` imports \`fetch_deposited_pdb\` from \`b.py\` inside the UI callback. When Streamlit executes \`b.py\` as \`__main__\`, this can import the application a second time.

**Required implementation:**

1. Move the RCSB coordinate-download functions to a small neutral module, for example \`docking/rcsb.py\` or a shared \`data_fetching.py\`.
2. Import that neutral module from both \`b.py\` and \`docking/ui.py\`.
3. Do not import the Streamlit application module from any docking module.
4. Test \`python -c "import docking.ui"\` and a Streamlit application start after this change.

### 8. Preserve coordinate provenance and biological assembly

- For an RCSB-based docking run, default to the selected RCSB Biological Assembly if available.
- Save the exact source PDB text before processing and record the PDB ID, assembly ID, URL, coordinate-source mode, and selected chains.
- Never remove chains through sequence de-duplication.
- For the \`1M17\` demo, obtain coordinates through the same shared coordinate-selection policy. Do not silently force deposited coordinates when the application’s selected mode is biological assembly.

### 9. Make receptor preparation explicit rather than silent

\`receptor_prep.py\` currently removes all non-protein residues except a small hard-coded allow-list.

- Keep reference-ligand removal explicit for redocking and save it to \`reference_ligand.pdb\` before removal.
- Do not silently delete cofactors, metal ions, prosthetic groups, or structural waters. Present an explicit preparation summary and require a recorded choice for each such component.
- Retaining a component is valid only when the selected PDBQT preparation path supports it and charges/atom types are valid.
- Save the complete preparation report in the manifest.

### 10. Harden engine-output parsing

- If an output model lacks \`REMARK VINA RESULT\`, raise a parse error; never create a score of \`0.0\`.
- Capture command arguments, exit code, stdout, stderr, and Vina \`--version\` in the manifest.
- Keep use of \`subprocess.run([...])\` as an argument list; do not replace it with \`shell=True\` or a concatenated shell string.
- Validate numeric user settings (\`exhaustiveness\`, grid dimensions, seed values, number of modes) before execution.

### 11. Provide a requirements file and controlled environment

Create \`requirements.txt\` (or an equivalent pinned environment file) with tested versions for at least:

\`\`\`text
biopython
numpy
pandas
streamlit
py3Dmol
rdkit
openbabel-wheel
\`\`\`

Document the tested Python version, platform, Vina version, installation command, and command to start the app. Do not claim a successful end-to-end run until the same environment imports all required modules.

## P2 — required before healthy-versus-mutant conclusions

### 12. Enforce matched-comparison conditions

Add a separate \`healthy_vs_mutant\` workflow. It must refuse to compare scores unless these match:

- same target domain and comparable functional state;
- same ligand identity, stereochemistry, protonation/tautomer state, and preparation toolchain;
- equivalent validated pocket and comparable grid policy;
- same engine version, scoring function, search settings, and seed list;
- same receptor preparation policy regarding waters, cofactors, metals, and ligands.

Report a score difference only as a **protocol-specific model ranking difference** with per-seed values and cluster support. It is not a measured affinity change or a predicted clinical response.

## Required automated tests

Add tests that do not need internet access or the full Vina executable:

1. Invalid/missing reference ligand stops before docking and creates a failed manifest.
2. No-reference site requires valid residue-defined input; it never uses origin-centred fallback coordinates.
3. RMSD refuses unequal atom counts and refuses absent atom mapping.
4. Correctly mapped identical ligand coordinates produce RMSD \`0.0\`.
5. A known translated coordinate set produces the expected RMSD after the documented alignment rule.
6. Three requested seeds cause three engine-adapter calls and three saved logs/outputs.
7. A missing Vina score remark raises a parser error.
8. Run directory and Vina binary resolution work when the process current directory is not the project root.
9. Manifest is written on both success and engine failure.
10. Importing \`docking.ui\` does not import \`b.py\`.

## Manual acceptance procedure

1. Create a clean environment from the new requirements file.
2. Start the Streamlit app from a directory other than the project root.
3. Select a protein-ligand co-crystal with a verified reference ligand and biological-assembly provenance.
4. Run redocking with three seeds.
5. Confirm that every seed log/output, reference ligand, source receptor, and complete manifest exists.
6. Confirm that the mapped heavy-atom RMSD calculation produces a reproducible value and that only \`<= 2.0 Å\` passes Gate A.
7. Inspect the docked complex in 3-D and verify that the ligand is in the intended pocket without obvious clashes.
8. Run one intentionally invalid-site case and confirm that Vina is never invoked.
9. Only after all of the above pass, run one matched healthy-versus-mutant experiment.

## Final rule

Do not scale to a ligand library or the project’s 500 protein pairs until the redocking pilot passes, manifests are complete, the RMSD mapping is valid, and the independent-seed protocol is implemented.

