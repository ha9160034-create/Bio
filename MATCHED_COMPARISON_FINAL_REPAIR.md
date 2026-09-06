# Remaining Repair: Healthy-versus-Mutant Docking Comparison

## Scope

This is the final remaining repair for the docking implementation. It applies only to:

    C:\Users\hassan\Desktop\MY_VScode\docking\comparison.py

The automated test suite has passed. Do not change the single-protein docking workflow, RMSD validation, or RCSB biological-assembly logic.

## Problem

The matched healthy-versus-mutant comparison currently does this:

1. Runs docking for healthy and mutant receptors.
2. Detects when either dominant pose cluster is not reproduced by at least two independent seeds.
3. Adds only a warning.
4. Still calculates and displays delta_score.

This is scientifically invalid. A delta score must not exist when either receptor has no reproducible dominant pose cluster.

A second issue is that the code uses 0.0 when a cluster is missing. This can create a fabricated delta score.

## Required behavior

After clustering the healthy and mutant poses:

1. Identify the dominant cluster for each receptor.
2. Require each dominant cluster to exist.
3. Require each dominant cluster to have seed_count >= 2.
4. If either requirement fails:
   - do not calculate delta_score;
   - do not return comparison results to the UI;
   - write a final failed manifest for both run directories;
   - include the exact failure reason and cluster summaries in both manifests;
   - raise ValueError with a clear user-facing message.

Use this logic before any delta-score calculation:

    h_top_cluster = h_clusters[0] if h_clusters else None
    m_top_cluster = m_clusters[0] if m_clusters else None

    if h_top_cluster is None or m_top_cluster is None:
        raise ValueError(
            "Matched comparison rejected: no dominant pose cluster was produced for one or both receptors."
        )

    if h_top_cluster.seed_count < 2 or m_top_cluster.seed_count < 2:
        raise ValueError(
            "Matched comparison rejected: the dominant pose was not reproduced by at least two independent seeds "
            "for both healthy and mutant receptors."
        )

Only after those checks are true may the code calculate:

    delta_score = mutant_cluster_median_score - healthy_cluster_median_score

Never use 0.0 as a substitute for missing cluster scores.

## Required manifest behavior

Before receptor preparation or Vina execution, create an initial manifest for both runs:

    healthy_run/manifest.initial.json
    mutant_run/manifest.initial.json

If preparation, docking, clustering, or comparison validation fails:

1. Finalize a manifest in both run directories.
2. Use status:

    run failed

3. Preserve:
   - run ID;
   - UTC timestamp;
   - target metadata;
   - ligand metadata;
   - receptor-preparation policy;
   - grid;
   - engine settings;
   - per-seed results already completed;
   - cluster summaries;
   - explicit failure reason.

Do not leave a partially completed comparison directory without a manifest.

## Result behavior after successful comparison

Only show the matched-comparison UI result when both clusters are reproducible.

Use this label without alteration:

> Protocol-specific docking-score ranking difference. It is not a measured affinity difference, drug-response prediction, or clinical conclusion.

The result must display:

- healthy dominant cluster median score;
- mutant dominant cluster median score;
- delta_score;
- supporting seed list and seed count for both clusters;
- the existing scientific disclaimer.

## Tests to add

Add at least these two automated tests to test_docking_repairs.py:

1. A comparison is rejected when the healthy or mutant dominant cluster has seed_count = 1.
2. A comparison is rejected when either receptor has no cluster.

The tests must verify that delta_score is never calculated or returned in either failure case.

## Acceptance criteria

The repair is complete only when:

- a non-reproducible healthy cluster rejects the comparison;
- a non-reproducible mutant cluster rejects the comparison;
- a missing cluster rejects the comparison;
- no result uses 0.0 as a missing cluster score;
- both run folders receive final failed manifests;
- a reproducible cluster for both receptors still permits delta_score calculation;
- the full test suite still ends with:

    OK

