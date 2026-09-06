# Matched Healthy-versus-Mutant Docking Comparison Module
import datetime
from pathlib import Path
from docking.models import GridBox, PoseCluster
from docking.storage import create_run_directory, write_initial_manifest
from docking.receptor_prep import prepare_receptor
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.binding_site import define_grid_from_residues
from docking.engine import run_vina_multi_seeds, get_vina_version
from docking.validation import cluster_poses_across_seeds
from docking.reporting import build_and_save_manifest, compute_file_sha256, DISCLAIMER_TEXT

COMPARISON_LABEL = (
    "Protocol-specific docking-score ranking difference. "
    "It is not a measured affinity difference, drug-response prediction, or clinical conclusion."
)

def validate_matched_conditions(
    healthy_domain: str,
    mutant_domain: str,
    healthy_ligand_hash: str,
    mutant_ligand_hash: str,
    healthy_grid_size: list,
    mutant_grid_size: list,
    healthy_waters_policy: str,
    mutant_waters_policy: str,
    healthy_seeds: list,
    mutant_seeds: list,
    healthy_exhaustiveness: int,
    mutant_exhaustiveness: int
) -> tuple[bool, str]:
    """Strictly verify matched comparison requirements before allowing delta calculations."""
    if healthy_domain != mutant_domain:
        return False, f"Domain mismatch: healthy domain '{healthy_domain}' != mutant domain '{mutant_domain}'."
    if healthy_ligand_hash != mutant_ligand_hash:
        return False, "Ligand mismatch: prepared ligand SHA-256 hashes must be identical."
    if healthy_grid_size != mutant_grid_size:
        return False, f"Grid box size mismatch: {healthy_grid_size} != {mutant_grid_size}."
    if healthy_waters_policy != mutant_waters_policy:
        return False, "Receptor preparation mismatch: waters policy differs."
    if healthy_seeds != mutant_seeds or len(healthy_seeds) < 3:
        return False, "Seed protocol mismatch: must run identical 3 independent seeds."
    if healthy_exhaustiveness != mutant_exhaustiveness:
        return False, f"Search exhaustiveness mismatch: {healthy_exhaustiveness} != {mutant_exhaustiveness}."
    return True, "Matched comparison conditions satisfied."

def verify_comparison_clusters(h_clusters: list[PoseCluster], m_clusters: list[PoseCluster]) -> float:
    """
    Verify that both healthy and mutant receptors produced a dominant pose cluster
    supported by at least two independent seeds.
    Refuses calculation and never substitutes 0.0 for missing clusters.
    Returns delta_score = mutant_median - healthy_median.
    """
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

    return round(m_top_cluster.median_score - h_top_cluster.median_score, 3)

def run_matched_docking_comparison(
    healthy_pdb_text: str,
    mutant_pdb_text: str,
    healthy_target_meta: dict,
    mutant_target_meta: dict,
    ligand_smiles: str,
    ligand_name: str,
    pocket_residues: list[str],
    target_domain: str = "Kinase Domain",
    padding: float = 8.0,
    seeds: list[int] = [42, 101, 2024],
    exhaustiveness: int = 8
) -> dict:
    """
    Execute controlled matched docking comparison between healthy and mutant structures.
    Refuses comparison unless all protocol settings match strictly and both dominant
    clusters have at least 2 independent supporting seeds.
    """
    if not healthy_pdb_text or not mutant_pdb_text:
        raise ValueError("Missing coordinates for healthy or mutant protein.")
    if not pocket_residues:
        raise ValueError("Docking was not run: define a validated binding site.")
    if len(seeds) < 3:
        raise ValueError("Matched comparison requires at least 3 independent seeds.")

    vina_ver = get_vina_version()
    h_created_at = healthy_target_meta.get("created_at_utc") or datetime.datetime.now(datetime.timezone.utc).isoformat()
    m_created_at = mutant_target_meta.get("created_at_utc") or datetime.datetime.now(datetime.timezone.utc).isoformat()

    # 1. Setup independent run directories
    h_run_id, h_run_dir = create_run_directory("dock_matched_healthy")
    m_run_id, m_run_dir = create_run_directory("dock_matched_mutant")

    # Save immutable source PDB files
    with open(h_run_dir / "receptor_source.pdb", "w", encoding="utf-8") as f:
        f.write(healthy_pdb_text)
    with open(m_run_dir / "receptor_source.pdb", "w", encoding="utf-8") as f:
        f.write(mutant_pdb_text)

    # 2. Define equivalent pocket grids
    h_grid = define_grid_from_residues(healthy_pdb_text, pocket_residues, padding=padding)
    m_grid = define_grid_from_residues(mutant_pdb_text, pocket_residues, padding=padding)

    grid_meta_h = {"center": [h_grid.center_x, h_grid.center_y, h_grid.center_z], "size": [h_grid.size_x, h_grid.size_y, h_grid.size_z]}
    grid_meta_m = {"center": [m_grid.center_x, m_grid.center_y, m_grid.center_z], "size": [m_grid.size_x, m_grid.size_y, m_grid.size_z]}
    engine_meta = {"name": "AutoDock Vina", "version": vina_ver, "seeds": seeds, "exhaustiveness": exhaustiveness, "num_modes": 9, "energy_range": 3.0}
    binding_site_meta = {"method": "residues", "residues": pocket_residues, "padding": padding}

    # 3. Create initial manifest for both runs BEFORE preparation or engine runs
    h_initial = {
        "schema_version": 1,
        "run_id": h_run_id,
        "status": "running",
        "created_at_utc": h_created_at,
        "target": healthy_target_meta,
        "binding_site": binding_site_meta,
        "ligand": {"name": ligand_name, "input_smiles": ligand_smiles},
        "grid": grid_meta_h,
        "engine": engine_meta
    }
    write_initial_manifest(h_run_dir, h_initial)

    m_initial = {
        "schema_version": 1,
        "run_id": m_run_id,
        "status": "running",
        "created_at_utc": m_created_at,
        "target": mutant_target_meta,
        "binding_site": binding_site_meta,
        "ligand": {"name": ligand_name, "input_smiles": ligand_smiles},
        "grid": grid_meta_m,
        "engine": engine_meta
    }
    write_initial_manifest(m_run_dir, m_initial)

    # Helper to finalize both manifests on any failure
    def _finalize_both_failed(err_msg: str, h_meta=None, m_meta=None, h_cl=None, m_cl=None):
        build_and_save_manifest(
            run_dir=h_run_dir,
            run_id=h_run_id,
            status="run failed",
            target_meta=healthy_target_meta,
            binding_site_meta=binding_site_meta,
            ligand_meta={"name": ligand_name, "input_smiles": ligand_smiles},
            prep_meta={},
            engine_meta=engine_meta,
            grid_meta=grid_meta_h,
            results_meta={"clusters": h_cl or [], "per_seed": h_meta or [], "failure_reason": err_msg},
            warnings=[err_msg]
        )
        build_and_save_manifest(
            run_dir=m_run_dir,
            run_id=m_run_id,
            status="run failed",
            target_meta=mutant_target_meta,
            binding_site_meta=binding_site_meta,
            ligand_meta={"name": ligand_name, "input_smiles": ligand_smiles},
            prep_meta={},
            engine_meta=engine_meta,
            grid_meta=grid_meta_m,
            results_meta={"clusters": m_cl or [], "per_seed": m_meta or [], "failure_reason": err_msg},
            warnings=[err_msg]
        )

    # 4. Receptor preparation with identical policy
    try:
        h_clean_pdb, h_rec_pdbqt, h_prep_rep = prepare_receptor(healthy_pdb_text, h_run_dir, keep_waters=False)
        m_clean_pdb, m_rec_pdbqt, m_prep_rep = prepare_receptor(mutant_pdb_text, m_run_dir, keep_waters=False)
    except Exception as e:
        _finalize_both_failed(f"Receptor preparation failure: {e}")
        raise

    # 5. Prepare ligand identically
    try:
        h_sdf, h_lig_pdbqt, h_lig_rep = prepare_ligand_from_smiles(ligand_smiles, ligand_name, h_run_dir, ph=7.4)
        m_sdf, m_lig_pdbqt, m_lig_rep = prepare_ligand_from_smiles(ligand_smiles, ligand_name, m_run_dir, ph=7.4)
    except Exception as e:
        _finalize_both_failed(f"Ligand preparation failure: {e}")
        raise

    h_lig_hash = compute_file_sha256(h_lig_pdbqt)
    m_lig_hash = compute_file_sha256(m_lig_pdbqt)

    # 6. Verify matched conditions
    valid, reason = validate_matched_conditions(
        healthy_domain=target_domain,
        mutant_domain=target_domain,
        healthy_ligand_hash=h_lig_hash,
        mutant_ligand_hash=m_lig_hash,
        healthy_grid_size=[h_grid.size_x, h_grid.size_y, h_grid.size_z],
        mutant_grid_size=[m_grid.size_x, m_grid.size_y, m_grid.size_z],
        healthy_waters_policy=h_prep_rep["waters_policy"],
        mutant_waters_policy=m_prep_rep["waters_policy"],
        healthy_seeds=seeds,
        mutant_seeds=seeds,
        healthy_exhaustiveness=exhaustiveness,
        mutant_exhaustiveness=exhaustiveness
    )
    if not valid:
        err_msg = f"Matched comparison rejected: {reason}"
        _finalize_both_failed(err_msg)
        raise ValueError(err_msg)

    # 7. Execute multi-seed Vina for Healthy
    try:
        h_poses, h_meta = run_vina_multi_seeds(
            receptor_pdbqt=h_rec_pdbqt,
            ligand_pdbqt=h_lig_pdbqt,
            grid=h_grid,
            output_dir=h_run_dir,
            seeds=seeds,
            exhaustiveness=exhaustiveness
        )
        h_clusters = cluster_poses_across_seeds(h_poses, rmsd_threshold=2.0)
    except Exception as e:
        _finalize_both_failed(f"Healthy docking engine failure: {e}")
        raise

    # 8. Execute multi-seed Vina for Mutant
    try:
        m_poses, m_meta = run_vina_multi_seeds(
            receptor_pdbqt=m_rec_pdbqt,
            ligand_pdbqt=m_lig_pdbqt,
            grid=m_grid,
            output_dir=m_run_dir,
            seeds=seeds,
            exhaustiveness=exhaustiveness
        )
        m_clusters = cluster_poses_across_seeds(m_poses, rmsd_threshold=2.0)
    except Exception as e:
        _finalize_both_failed(f"Mutant docking engine failure: {e}", h_meta=h_meta, h_cl=h_clusters)
        raise

    # 9. Strict reproducibility check before calculating delta_score
    try:
        delta_score = verify_comparison_clusters(h_clusters, m_clusters)
    except ValueError as val_err:
        _finalize_both_failed(str(val_err), h_meta=h_meta, m_meta=m_meta, h_cl=h_clusters, m_cl=m_clusters)
        raise

    h_top_cluster = h_clusters[0]
    m_top_cluster = m_clusters[0]

    # 10. Save successful manifests for both runs
    build_and_save_manifest(
        run_dir=h_run_dir,
        run_id=h_run_id,
        status="exploratory—no reference validation",
        target_meta=healthy_target_meta,
        binding_site_meta=binding_site_meta,
        ligand_meta=h_lig_rep,
        prep_meta=h_prep_rep,
        engine_meta=engine_meta,
        grid_meta=grid_meta_h,
        results_meta={
            "top_score": min(p.score for p in h_poses),
            "clusters_count": len(h_clusters),
            "clusters": [
                {
                    "cluster_id": cl.cluster_id,
                    "supporting_seeds": cl.supporting_seeds,
                    "seed_count": cl.seed_count,
                    "top_score": cl.top_score,
                    "median_score": cl.median_score,
                    "representative_pose_source": cl.representative_pose_source
                }
                for cl in h_clusters
            ],
            "per_seed": h_meta
        },
        warnings=[]
    )
    build_and_save_manifest(
        run_dir=m_run_dir,
        run_id=m_run_id,
        status="exploratory—no reference validation",
        target_meta=mutant_target_meta,
        binding_site_meta=binding_site_meta,
        ligand_meta=m_lig_rep,
        prep_meta=m_prep_rep,
        engine_meta=engine_meta,
        grid_meta=grid_meta_m,
        results_meta={
            "top_score": min(p.score for p in m_poses),
            "clusters_count": len(m_clusters),
            "clusters": [
                {
                    "cluster_id": cl.cluster_id,
                    "supporting_seeds": cl.supporting_seeds,
                    "seed_count": cl.seed_count,
                    "top_score": cl.top_score,
                    "median_score": cl.median_score,
                    "representative_pose_source": cl.representative_pose_source
                }
                for cl in m_clusters
            ],
            "per_seed": m_meta
        },
        warnings=[]
    )

    return {
        "comparison_label": COMPARISON_LABEL,
        "delta_score": delta_score,
        "healthy": {
            "run_id": h_run_id,
            "target_id": healthy_target_meta.get("pdb_id"),
            "dominant_cluster_top_score": h_top_cluster.top_score,
            "dominant_cluster_median": h_top_cluster.median_score,
            "dominant_cluster_seeds": h_top_cluster.supporting_seeds,
            "seed_count": h_top_cluster.seed_count,
            "clean_pdb": h_clean_pdb,
            "poses": h_poses,
            "clusters": h_clusters
        },
        "mutant": {
            "run_id": m_run_id,
            "target_id": mutant_target_meta.get("pdb_id"),
            "dominant_cluster_top_score": m_top_cluster.top_score,
            "dominant_cluster_median": m_top_cluster.median_score,
            "dominant_cluster_seeds": m_top_cluster.supporting_seeds,
            "seed_count": m_top_cluster.seed_count,
            "clean_pdb": m_clean_pdb,
            "poses": m_poses,
            "clusters": m_clusters
        },
        "warnings": [],
        "disclaimer": DISCLAIMER_TEXT
    }
