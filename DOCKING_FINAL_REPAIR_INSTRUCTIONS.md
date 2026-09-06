# Final Docking Repair Instructions

## Purpose

Apply these remaining repairs to C:\Users\hassan\Desktop\MY_VScode. Do not claim clinical benefit, safety, dose, affinity, Kd, Ki, IC50, or treatment efficacy.

## P0: RMSD mapping must inspect every valid mapping

File: docking/validation.py

The current redocking calculation gets every reference and pose mapping but only uses the first pose mapping. This can produce a wrong RMSD for symmetric ligands.

Required changes:

1. Iterate through every pair of reference and pose mappings.
2. Calculate heavy-atom RMSD for each complete pair.
3. Keep the smallest valid RMSD.
4. Refuse incomplete mappings or unequal heavy-atom counts.
5. Save mapping method, reference and pose heavy-atom counts, matched count, number of mappings evaluated, selected mapping pairs, and final RMSD in the manifest.
6. Do not fit or rotate the ligand before RMSD. The docked and co-crystal ligand already share the receptor coordinate frame in redocking.

Use this structure:

    for m_ref in match_ref_all:
        for m_pose in match_pose_all:
            calculate RMSD for this complete mapping pair

If full identity cannot be established, return validation unavailable, never partial RMSD.

## P0: Gate A needs reproducible cluster support

File: docking/ui.py

Validation passed is permitted only when all are true:

- co-crystal reference ligand mode was used;
- three independent seeds were executed;
- mapped heavy-atom RMSD is at most 2.0 Angstrom;
- the validated pose belongs to a cluster supported by at least two distinct seeds.

Required statuses:

- RMSD above 2.0 Angstrom: validation failed.
- Mapping unavailable: validation unavailable.
- One seed: exploratory—no reference validation.
- Valid RMSD but less than two-seed cluster support: validation failed, with an explicit warning.

Check the cluster containing the actual validated pose; do not assume the first cluster is that cluster.

## P0: Create a new complete validation run

The existing test_repaired_pipeline_20260906_181005_98aafc run is not a valid final record:

- redocking RMSD is 8.547 Angstrom;
- its status is completed, not validation failed;
- receptor_source.pdb and reference_ligand.pdb are absent;
- created_at_utc is empty.

Do not edit or relabel this historical run. Treat it as a failed exploratory development run.

After repair, create a new run through the Streamlit UI with:

    manifest.initial.json
    manifest.json
    receptor_source.pdb
    reference_ligand.pdb
    receptor_clean.pdb
    receptor.pdbqt
    ligand SDF and PDBQT
    out_seed_42.pdbqt
    out_seed_101.pdbqt
    out_seed_2024.pdbqt
    vina_seed_42.log
    vina_seed_101.log
    vina_seed_2024.log

Final statuses may only be:

    validation passed
    validation failed
    validation unavailable
    exploratory—no reference validation
    run failed

Do not emit the vague status completed.

## P1: Correct coordinate provenance

File: docking/ui.py

The 1M17 demo can fall back to deposited coordinates after biological-assembly download fails, but still records biological assembly. Record the source actually used.

Required policy:

1. Attempt biological assembly 1.
2. If it exists, set coordinate_source to rcsb_biological_assembly and assembly_id to 1.
3. If it fails and deposited coordinates are used, set coordinate_source to deposited and assembly_id to null.
4. Record source URL, PDB ID, assembly ID, selected chains, and non-empty UTC creation time.
5. Uploaded structures must use coordinate_source upload and record the original file name.

## P1: Complete every manifest

File: docking/ui.py and docking/reporting.py

Required changes:

1. Add mapping_details to results.redocking.
2. Save validated pose seed, rank, cluster ID, and output file.
3. Derive top_score as the minimum score across every pose from every seed, not poses[0].
4. Preserve target, ligand, grid, binding site, engine, and created_at_utc when preparation or engine failures finalize a manifest.
5. Save per-seed command, Vina version, exit code, log file, output file, pose count, and best score.
6. Save per-cluster seed support, median score, top score, and representative pose source.

## P1: Record receptor-preparation decisions

File: docking/receptor_prep.py

1. Preserve the reference ligand in reference_ligand.pdb before receptor removal.
2. Record each removed ligand and retained cofactor with residue name, chain, and residue number.
3. Display preparation choices in the UI.
4. Mark removed waters as a model assumption.
5. Do not claim retained cofactors are parameterized unless their PDBQT types and charges are verified.

## P2: Implement matched healthy-versus-mutant docking

The application currently has no docking comparison workflow for healthy versus mutant structures. Sequence/SASA comparison is not docking comparison.

Create docking/comparison.py. Refuse comparison unless both runs have:

- same target domain and comparable functional state;
- exact same prepared ligand state and ligand hash;
- equivalent validated functional pocket and grid policy;
- same receptor-preparation policy;
- same Vina version, search settings, seed list, scoring protocol, and number of modes.

For an accepted comparison:

1. Create independent healthy and mutant run directories/manifests.
2. Dock the same prepared ligand with three seeds per receptor.
3. Require a reproducible cluster in each receptor.
4. Report cluster-level score distributions and medians.
5. Calculate delta_score as mutant cluster median minus healthy cluster median.
6. Label it: Protocol-specific docking-score ranking difference. It is not a measured affinity difference, drug-response prediction, or clinical conclusion.

## P3: Tests and environment

The project now has requirements.txt and test_docking_repairs.py, but the inspected Python runtime could not import OpenBabel. Create the declared project environment and run:

    python -m pip install -r requirements.txt
    python -m unittest -v test_docking_repairs.py

Add tests for:

1. all symmetry-equivalent mapping pairs;
2. mapping failure and unequal atom counts;
3. RMSD threshold with one-seed cluster cannot pass;
4. RMSD threshold with two-seed cluster can pass;
5. deposited fallback provenance;
6. global top score across seeds;
7. non-empty created_at_utc and stored mapping details;
8. failed manifests retaining initial metadata;
9. no completed validation status;
10. healthy/mutant setting mismatch rejection.

## Acceptance checklist

- [ ] All tests pass in the declared environment.
- [ ] A fresh UI redocking run contains every required artifact.
- [ ] RMSD mapping evaluates all valid mapping pairs.
- [ ] Gate A requires RMSD at most 2.0 Angstrom and at least two-seed pose reproduction.
- [ ] The final manifest has correct provenance and full validation details.
- [ ] A healthy/mutant docking result is available only through the matched-comparison workflow.
- [ ] The scientific limitation disclaimer remains visible beside results.

