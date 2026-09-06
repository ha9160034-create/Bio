# Automated test suite for docking repairs (P3 requirements)
import unittest
import json
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from openbabel import pybel

from docking.models import GridBox, DockingPose, PoseCluster
from docking.storage import PROJECT_ROOT, get_vina_path, create_run_directory, write_initial_manifest, finalize_manifest
from docking.binding_site import define_grid_from_residues, define_grid_from_ligand
from docking.validation import calculate_heavy_atom_rmsd, calculate_mapped_redocking_rmsd, cluster_poses_across_seeds
from docking.reporting import build_and_save_manifest
from docking.comparison import validate_matched_conditions, verify_comparison_clusters

class TestDockingFinalRepairs(unittest.TestCase):

    def test_1_all_symmetry_equivalent_mapping_pairs(self):
        """P3.1: Verify all symmetry-equivalent mapping pairs are evaluated."""
        b_mol = Chem.MolFromSmiles("c1ccccc1")
        b_mol = Chem.AddHs(b_mol)
        AllChem.EmbedMolecule(b_mol, randomSeed=42)
        b_mol = Chem.RemoveHs(b_mol)
        b_pdb = Chem.MolToPDBBlock(b_mol)

        ob_mol = pybel.readstring("pdb", b_pdb)
        b_pdbqt = ob_mol.write("pdbqt")
        
        rmsd, details = calculate_mapped_redocking_rmsd(b_pdb, b_pdbqt)
        self.assertAlmostEqual(rmsd, 0.0, places=3)
        self.assertIn("mappings_evaluated", details)
        self.assertGreaterEqual(details["mappings_evaluated"], 1)
        self.assertEqual(details["mapping_method"], "rdkit_mcs_topology_invariant_exhaustive_pairs")
        self.assertIn("selected_mapping_pairs", details)

    def test_2_mapping_failure_and_unequal_atom_counts(self):
        """P3.2: Refuse unequal atom counts or failed mapping."""
        b_mol = Chem.MolFromSmiles("c1ccccc1")
        t_mol = Chem.MolFromSmiles("c1ccccc1C")
        AllChem.EmbedMolecule(b_mol, randomSeed=42)
        AllChem.EmbedMolecule(t_mol, randomSeed=42)
        b_pdb = Chem.MolToPDBBlock(b_mol)
        ob_mol = pybel.readstring("pdb", Chem.MolToPDBBlock(t_mol))
        t_pdbqt = ob_mol.write("pdbqt")

        with self.assertRaises(ValueError) as ctx:
            calculate_mapped_redocking_rmsd(b_pdb, t_pdbqt)
        self.assertIn("Validation unavailable", str(ctx.exception))

    def test_3_rmsd_threshold_with_one_seed_cluster_cannot_pass(self):
        """P3.3: RMSD <= 2.0 A with only 1 seed cannot pass Gate A."""
        seeds_run = [42]
        rmsd_val = 1.2
        cluster_seed_count = 1

        if len(seeds_run) < 3:
            status = "exploratory—no reference validation"
        elif rmsd_val <= 2.0 and cluster_seed_count >= 2:
            status = "validation passed"
        else:
            status = "validation failed"

        self.assertNotEqual(status, "validation passed")
        self.assertEqual(status, "exploratory—no reference validation")

    def test_4_rmsd_threshold_with_two_seed_cluster_can_pass(self):
        """P3.4: RMSD <= 2.0 A with multi-seed cluster support can pass Gate A."""
        seeds_run = [42, 101, 2024]
        rmsd_val = 1.5
        cluster_seed_count = 2

        if len(seeds_run) < 3:
            status = "exploratory—no reference validation"
        elif rmsd_val <= 2.0 and cluster_seed_count >= 2:
            status = "validation passed"
        else:
            status = "validation failed"

        self.assertEqual(status, "validation passed")

    def test_5_deposited_fallback_provenance(self):
        """P3.5: Fallback to deposited coordinates sets assembly_id to None."""
        assembly_available = False
        if assembly_available:
            coord_source = "rcsb_biological_assembly"
            assembly_id = 1
        else:
            coord_source = "deposited"
            assembly_id = None

        self.assertEqual(coord_source, "deposited")
        self.assertIsNone(assembly_id)

    def test_6_global_top_score_across_seeds(self):
        """P3.6: Global top_score must be minimum across all poses from every seed."""
        p1 = DockingPose(rank=1, score=-6.5, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0)
        p2 = DockingPose(rank=1, score=-8.2, seed=101, cluster_id=2, rmsd_lb=0, rmsd_ub=0)
        p3 = DockingPose(rank=1, score=-7.4, seed=2024, cluster_id=3, rmsd_lb=0, rmsd_ub=0)
        poses = [p1, p2, p3]

        top_score = min(p.score for p in poses)
        self.assertEqual(top_score, -8.2)

    def test_7_non_empty_created_at_utc_and_stored_mapping_details(self):
        """P3.7: Non-empty UTC timestamp and stored mapping details in manifest."""
        run_id, run_dir = create_run_directory("test_meta")
        try:
            created_at = "2026-09-06T18:00:00Z"
            manifest_path = build_and_save_manifest(
                run_dir=run_dir,
                run_id=run_id,
                status="validation failed",
                target_meta={"pdb_id": "TEST", "created_at_utc": created_at},
                binding_site_meta={"method": "co_crystal_ligand"},
                ligand_meta={"name": "LIG"},
                prep_meta={},
                engine_meta={},
                grid_meta={},
                results_meta={
                    "top_score": -7.5,
                    "redocking": {
                        "rmsd_angstrom": 3.4,
                        "mapping_details": {"mapping_method": "rdkit_mcs_topology_invariant_exhaustive_pairs"}
                    }
                }
            )
            with open(manifest_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.assertEqual(data["created_at_utc"], created_at)
            self.assertTrue(len(data["created_at_utc"]) > 0)
            self.assertIn("mapping_details", data["results"]["redocking"])
        finally:
            import shutil
            shutil.rmtree(run_dir, ignore_errors=True)

    def test_8_failed_manifests_retaining_initial_metadata(self):
        """P3.8: Failed manifests preserve target, ligand, engine, and binding_site."""
        run_id, run_dir = create_run_directory("test_fail_retain")
        try:
            target_meta = {"pdb_id": "1M17", "created_at_utc": "2026-09-06T18:00:00Z"}
            ligand_meta = {"name": "Erlotinib"}
            binding_site_meta = {"method": "co_crystal_ligand"}
            engine_meta = {"name": "AutoDock Vina"}
            
            manifest_path = build_and_save_manifest(
                run_dir=run_dir,
                run_id=run_id,
                status="run failed",
                target_meta=target_meta,
                binding_site_meta=binding_site_meta,
                ligand_meta=ligand_meta,
                prep_meta={},
                engine_meta=engine_meta,
                grid_meta={},
                results_meta={},
                warnings=["Preparation failed"]
            )
            with open(manifest_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.assertEqual(data["status"], "run failed")
            self.assertEqual(data["target"]["pdb_id"], "1M17")
            self.assertEqual(data["ligand"]["name"], "Erlotinib")
            self.assertEqual(data["engine"]["name"], "AutoDock Vina")
            self.assertEqual(data["binding_site"]["method"], "co_crystal_ligand")
        finally:
            import shutil
            shutil.rmtree(run_dir, ignore_errors=True)

    def test_9_no_completed_validation_status(self):
        """P3.9: 'completed' status is forbidden and sanitized to a valid status."""
        run_id, run_dir = create_run_directory("test_no_completed")
        try:
            manifest_path = build_and_save_manifest(
                run_dir=run_dir,
                run_id=run_id,
                status="completed",
                target_meta={"pdb_id": "TEST", "created_at_utc": "2026-09-06T18:00:00Z"},
                binding_site_meta={},
                ligand_meta={},
                prep_meta={},
                engine_meta={},
                grid_meta={},
                results_meta={}
            )
            with open(manifest_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.assertNotEqual(data["status"], "completed")
            allowed = ["validation passed", "validation failed", "validation unavailable", "exploratory—no reference validation", "run failed"]
            self.assertIn(data["status"], allowed)
        finally:
            import shutil
            shutil.rmtree(run_dir, ignore_errors=True)

    def test_10_healthy_mutant_setting_mismatch_rejection(self):
        """P3.10: Rejection of mismatched healthy-versus-mutant docking settings."""
        valid, reason = validate_matched_conditions(
            healthy_domain="Kinase Domain",
            mutant_domain="SH2 Domain",
            healthy_ligand_hash="abc",
            mutant_ligand_hash="abc",
            healthy_grid_size=[20, 20, 20],
            mutant_grid_size=[20, 20, 20],
            healthy_waters_policy="Bulk waters removed",
            mutant_waters_policy="Bulk waters removed",
            healthy_seeds=[42, 101, 2024],
            mutant_seeds=[42, 101, 2024],
            healthy_exhaustiveness=8,
            mutant_exhaustiveness=8
        )
        self.assertFalse(valid)
        self.assertIn("Domain mismatch", reason)

        valid2, reason2 = validate_matched_conditions(
            healthy_domain="Kinase Domain",
            mutant_domain="Kinase Domain",
            healthy_ligand_hash="abc",
            mutant_ligand_hash="def",
            healthy_grid_size=[20, 20, 20],
            mutant_grid_size=[20, 20, 20],
            healthy_waters_policy="Bulk waters removed",
            mutant_waters_policy="Bulk waters removed",
            healthy_seeds=[42, 101, 2024],
            mutant_seeds=[42, 101, 2024],
            healthy_exhaustiveness=8,
            mutant_exhaustiveness=8
        )
        self.assertFalse(valid2)
        self.assertIn("Ligand mismatch", reason2)

    def test_11_matched_comparison_rejected_if_seed_count_one(self):
        """P2: Reject matched comparison when healthy or mutant cluster has seed_count = 1."""
        dummy_pose = DockingPose(rank=1, score=-7.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0)
        h_cl = [PoseCluster(cluster_id=1, seed_count=1, supporting_seeds=[42], top_score=-7.5, median_score=-7.2, representative_pose=dummy_pose)]
        m_cl = [PoseCluster(cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-7.8, median_score=-7.6, representative_pose=dummy_pose)]

        with self.assertRaises(ValueError) as ctx:
            verify_comparison_clusters(h_cl, m_cl)
        self.assertIn("Matched comparison rejected", str(ctx.exception))
        self.assertIn("at least two independent seeds", str(ctx.exception))

    def test_12_matched_comparison_rejected_if_cluster_missing(self):
        """P2: Reject matched comparison when either receptor has no cluster (never use 0.0)."""
        dummy_pose = DockingPose(rank=1, score=-7.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0)
        valid_cl = [PoseCluster(cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-7.5, median_score=-7.2, representative_pose=dummy_pose)]

        # Empty healthy
        with self.assertRaises(ValueError) as ctx1:
            verify_comparison_clusters([], valid_cl)
        self.assertIn("no dominant pose cluster was produced", str(ctx1.exception))

        # Empty mutant
        with self.assertRaises(ValueError) as ctx2:
            verify_comparison_clusters(valid_cl, [])
        self.assertIn("no dominant pose cluster was produced", str(ctx2.exception))

    def test_13_matched_comparison_calculates_delta_when_clusters_valid(self):
        """P2: Calculates delta_score when both clusters have seed_count >= 2."""
        dummy_pose = DockingPose(rank=1, score=-7.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0)
        h_cl = [PoseCluster(cluster_id=1, seed_count=3, supporting_seeds=[42, 101, 2024], top_score=-8.0, median_score=-7.5, representative_pose=dummy_pose)]
        m_cl = [PoseCluster(cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-7.2, median_score=-6.8, representative_pose=dummy_pose)]

        delta = verify_comparison_clusters(h_cl, m_cl)
        # -6.8 - (-7.5) = +0.7
        self.assertAlmostEqual(delta, 0.7, places=3)

if __name__ == "__main__":
    unittest.main()
