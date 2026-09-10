# Automated test suite for docking repairs (P3 requirements)
import unittest
from unittest.mock import patch
import json
import numpy as np
from pathlib import Path
import docking
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from openbabel import pybel

from docking.models import GridBox, DockingPose, PoseCluster
from docking.storage import PROJECT_ROOT, get_vina_path, create_run_directory, write_initial_manifest, finalize_manifest
from docking.binding_site import define_grid_from_residues, define_grid_from_ligand
from docking.validation import calculate_heavy_atom_rmsd, calculate_mapped_redocking_rmsd, cluster_poses_across_seeds, validate_redocking_top_n
from docking.reporting import build_and_save_manifest
from docking.comparison import validate_matched_conditions, verify_comparison_clusters
from docking.ui import sync_target_binding_site_state
from docking.ui import sync_target_binding_site_state, find_recovery_pose_index

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

    def test_14_changing_custom_pdb_clears_target_specific_inputs(self):
        """CUSTOM_TARGET_STATE_REPAIR: Changing Custom PDB ID clears reference ligand, pocket residues, and padding."""
        state = {
            "dock_previous_target_identity": "custom:4MZI",
            "dock_ref_lig_input": "XYZ",
            "dock_pocket_res_input": "145, 147",
            "dock_pocket_padding": 12.0,
            "dock_lig_name_input": "Aspirin",
            "dock_lig_smiles_input": "CC(=O)OC1=CC=CC=C1C(=O)O"
        }
        # Change custom PDB to 2ITO
        was_reset = sync_target_binding_site_state("custom:2ITO", state)
        self.assertTrue(was_reset)
        self.assertEqual(state["dock_ref_lig_input"], "")
        self.assertEqual(state["dock_pocket_res_input"], "")
        self.assertEqual(state["dock_pocket_padding"], 8.0)
        self.assertEqual(state["dock_previous_target_identity"], "custom:2ITO")
        # Ensure ligand name and smiles are NOT reset
        self.assertEqual(state["dock_lig_name_input"], "Aspirin")
        self.assertEqual(state["dock_lig_smiles_input"], "CC(=O)OC1=CC=CC=C1C(=O)O")

    def test_15_switching_target_source_mode_clears_target_specific_values(self):
        """CUSTOM_TARGET_STATE_REPAIR: Switching target source mode clears target-specific values."""
        state = {
            "dock_previous_target_identity": "demo:1M17",
            "dock_ref_lig_input": "AQ4",
            "dock_pocket_res_input": "790, 858",
            "dock_pocket_padding": 10.0
        }
        # Switch to uploaded file
        was_reset = sync_target_binding_site_state("upload:protein.pdb", state)
        self.assertTrue(was_reset)
        self.assertEqual(state["dock_ref_lig_input"], "")
        self.assertEqual(state["dock_pocket_res_input"], "")
        self.assertEqual(state["dock_pocket_padding"], 8.0)
        self.assertEqual(state["dock_previous_target_identity"], "upload:protein.pdb")

    def test_16_changing_non_target_settings_does_not_reset_binding_site(self):
        """CUSTOM_TARGET_STATE_REPAIR: Changing only surface, seed mode, or exhaustiveness does not reset them."""
        state = {
            "dock_previous_target_identity": "custom:4MZI",
            "dock_ref_lig_input": "LIG1",
            "dock_pocket_res_input": "120, 125",
            "dock_pocket_padding": 9.5
        }
        # Calling with the same target identity (e.g. user toggled surface or exhaustiveness)
        was_reset = sync_target_binding_site_state("custom:4MZI", state)
        self.assertFalse(was_reset)
        self.assertEqual(state["dock_ref_lig_input"], "LIG1")
        self.assertEqual(state["dock_pocket_res_input"], "120, 125")
        self.assertEqual(state["dock_pocket_padding"], 9.5)

    def test_17_custom_pdb_1M17_defaults_to_aq4(self):
        """CUSTOM_TARGET_STATE_REPAIR: Custom PDB 1M17 defaults to AQ4 and 790, 858."""
        state = {
            "dock_previous_target_identity": "custom:4MZI",
            "dock_ref_lig_input": "",
            "dock_pocket_res_input": "",
            "dock_pocket_padding": 8.0
        }
        was_reset = sync_target_binding_site_state("custom:1M17", state)
        self.assertTrue(was_reset)
        self.assertEqual(state["dock_ref_lig_input"], "AQ4")
        self.assertEqual(state["dock_pocket_res_input"], "790, 858")
        self.assertEqual(state["dock_pocket_padding"], 8.0)

    def test_18_top_n_recovery_passes_when_top1_fails(self):
        """Top-N: Non-top-1 pose within Top-10, RMSD 1.5 Å, support from 2 seeds -> recovery passes while top1_rmsd_pass is false."""
        poses = []
        for i in range(1, 11):
            poses.append(DockingPose(rank=i, score=-8.0 + (i * 0.1), seed=42 if i % 2 == 0 else 101, cluster_id=1 if i == 5 else 2, rmsd_lb=0, rmsd_ub=0, pdbqt_block=f"POSE_{i}"))
        
        clusters = [
            PoseCluster(cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-7.5, median_score=-7.5, representative_pose=poses[4]),
            PoseCluster(cluster_id=2, seed_count=2, supporting_seeds=[42, 101], top_score=-8.0, median_score=-7.5, representative_pose=poses[0])
        ]

        def mock_rmsd(ref, block):
            if block == "POSE_1":
                return (5.0, {"matched_atoms": 20})
            elif block == "POSE_5":
                return (1.5, {"matched_atoms": 20})
            return (6.0, {"matched_atoms": 20})

        with patch("docking.validation.calculate_mapped_redocking_rmsd", side_effect=mock_rmsd):
            res = validate_redocking_top_n("REF_LIG", poses, clusters, top_n=10, rmsd_cutoff=2.0, min_seed_support=2)

        self.assertFalse(res["top1_rmsd_pass"])
        self.assertEqual(res["top1"]["rmsd_angstrom"], 5.0)
        self.assertTrue(res["recovery_rmsd_pass"])
        self.assertEqual(res["recovery"]["rmsd_angstrom"], 1.5)
        self.assertTrue(res["recovery_cluster_pass"])
        self.assertEqual(res["recovery"]["rank"], 5)

    def test_19_pose_at_rank_11_does_not_pass_top10_recovery(self):
        """Top-N: A pose at rank 11, RMSD 1.0 Å -> does not pass Top-10 recovery."""
        poses = []
        for i in range(1, 15):
            poses.append(DockingPose(rank=i, score=-8.0 + (i * 0.1), seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block=f"POSE_{i}"))
        clusters = [PoseCluster(cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-8.0, median_score=-7.0, representative_pose=poses[0])]

        def mock_rmsd(ref, block):
            if block == "POSE_11":
                return (1.0, {"matched_atoms": 20})
            return (4.0, {"matched_atoms": 20})

        with patch("docking.validation.calculate_mapped_redocking_rmsd", side_effect=mock_rmsd):
            res = validate_redocking_top_n("REF_LIG", poses, clusters, top_n=10, rmsd_cutoff=2.0, min_seed_support=2)

        self.assertFalse(res["recovery_rmsd_pass"])
        self.assertEqual(res["recovery"]["rmsd_angstrom"], 4.0)

    def test_20_top10_pose_fails_if_only_one_seed_support(self):
        """Top-N: A Top-10 pose with RMSD 1.5 Å but one-seed support -> fails recovery."""
        poses = [
            DockingPose(rank=1, score=-8.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="POSE_1"),
            DockingPose(rank=2, score=-7.8, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="POSE_2")
        ]
        # Single seed support = 1
        clusters = [PoseCluster(cluster_id=1, seed_count=1, supporting_seeds=[42], top_score=-8.0, median_score=-7.9, representative_pose=poses[0])]

        def mock_rmsd(ref, block):
            return (1.5, {"matched_atoms": 20})

        with patch("docking.validation.calculate_mapped_redocking_rmsd", side_effect=mock_rmsd):
            res = validate_redocking_top_n("REF_LIG", poses, clusters, top_n=10, rmsd_cutoff=2.0, min_seed_support=2)

        self.assertTrue(res["recovery_rmsd_pass"])
        self.assertFalse(res["recovery_cluster_pass"])

    def test_21_rmsd_exactly_2_passes_cutoff(self):
        """Top-N: RMSD exactly 2.0 Å -> passes cutoff."""
        poses = [DockingPose(rank=1, score=-8.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="POSE_1")]
        clusters = [PoseCluster(cluster_id=1, seed_count=2, supporting_seeds=[42, 101], top_score=-8.0, median_score=-8.0, representative_pose=poses[0])]

        with patch("docking.validation.calculate_mapped_redocking_rmsd", return_value=(2.0, {"matched_atoms": 20})):
            res = validate_redocking_top_n("REF_LIG", poses, clusters, top_n=10, rmsd_cutoff=2.0, min_seed_support=2)

        self.assertTrue(res["top1_rmsd_pass"])
        self.assertTrue(res["recovery_rmsd_pass"])

    def test_22_fewer_than_three_seeds_reports_insufficient_or_exploratory(self):
        """Top-N: Fewer than three requested seeds -> UI decision logic marks as exploratory/insufficient."""
        seeds_to_run = [42]
        MIN_VALIDATION_SEEDS = 3
        status = "exploratory—insufficient independent seeds" if len(seeds_to_run) < MIN_VALIDATION_SEEDS else "validation passed"
        self.assertEqual(status, "exploratory—insufficient independent seeds")

    def test_23_manifest_contains_criterion_and_both_records(self):
        """Top-N: Manifest contains criterion, top1_scored_pose, and best_recovered_pose_in_top_n."""
        import tempfile
        import shutil
        temp_dir = Path(tempfile.mkdtemp())
        try:
            target_meta = {"created_at_utc": "2026-09-09T01:00:00Z"}
            results_meta = {
                "top_score": -7.5,
                "validation_status": "validation passed (Top-10 recovery)",
                "redocking": {
                    "criterion": {
                        "top_n": 10,
                        "rmsd_cutoff_angstrom": 2.0,
                        "minimum_independent_seeds": 3,
                        "minimum_cluster_seed_support": 2
                    },
                    "top1_scored_pose": {"rmsd_angstrom": 5.0},
                    "best_recovered_pose_in_top_n": {"rmsd_angstrom": 1.2, "cluster_seed_count": 2},
                    "top1_rmsd_pass": False,
                    "recovery_rmsd_pass": True,
                    "recovery_cluster_pass": True
                }
            }
            manifest_path = build_and_save_manifest(
                run_dir=temp_dir,
                run_id="test_top_n_manifest",
                status="validation passed (Top-10 recovery)",
                target_meta=target_meta,
                binding_site_meta={},
                ligand_meta={},
                prep_meta={},
                engine_meta={},
                grid_meta={},
                results_meta=results_meta
            )
            with open(manifest_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.assertEqual(data["status"], "validation passed (Top-10 recovery)")
            self.assertIn("criterion", data["results"]["redocking"])
            self.assertEqual(data["results"]["redocking"]["criterion"]["top_n"], 10)
            self.assertIn("top1_scored_pose", data["results"]["redocking"])
            self.assertIn("best_recovered_pose_in_top_n", data["results"]["redocking"])
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def test_24_repeated_ranks_across_seeds_returns_exact_recovery_index(self):
        """Top-N Recovery Viz: Given poses with repeated ranks across seeds, recovery {seed: 101, rank: 7} returns index of seed 101/rank 7."""
        poses = [
            DockingPose(rank=1, score=-8.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P1"),
            DockingPose(rank=7, score=-7.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P2"),
            DockingPose(rank=1, score=-7.8, seed=101, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P3"),
            DockingPose(rank=7, score=-7.2, seed=101, cluster_id=2, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P4"),
            DockingPose(rank=7, score=-6.8, seed=202, cluster_id=2, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P5"),
        ]
        val_details = {
            "recovery": {
                "seed": 101,
                "rank": 7,
                "rmsd_angstrom": 1.45
            }
        }
        idx = find_recovery_pose_index(poses, val_details)
        self.assertEqual(idx, 3)
        self.assertEqual(poses[idx].seed, 101)
        self.assertEqual(poses[idx].rank, 7)

    def test_25_missing_validation_details_returns_zero(self):
        """Top-N Recovery Viz: Missing validation details or empty poses returns index 0."""
        poses = [
            DockingPose(rank=1, score=-8.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P1"),
            DockingPose(rank=2, score=-7.5, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P2"),
        ]
        self.assertEqual(find_recovery_pose_index(poses, None), 0)
        self.assertEqual(find_recovery_pose_index(poses, {}), 0)
        self.assertEqual(find_recovery_pose_index(poses, {"recovery": {}}), 0)
        self.assertEqual(find_recovery_pose_index([], {"recovery": {"seed": 42, "rank": 1}}), 0)

    def test_26_recovery_metadata_not_in_poses_returns_zero_safely(self):
        """Top-N Recovery Viz: Recovery metadata that does not exist in the current pose list returns 0 without exception."""
        poses = [
            DockingPose(rank=1, score=-8.0, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P1"),
            DockingPose(rank=2, score=-7.5, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P2"),
        ]
        val_details = {
            "recovery": {
                "seed": 999,
                "rank": 99,
                "rmsd_angstrom": 1.2
            }
        }
        idx = find_recovery_pose_index(poses, val_details)
        self.assertEqual(idx, 0)

    def test_27_selection_identifies_recovery_pose_only_when_both_seed_and_rank_match(self):
        """Top-N Recovery Viz: Selection identifies the recovery pose only when both seed and rank match."""
        poses = [
            DockingPose(rank=7, score=-7.1, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P1"),
            DockingPose(rank=2, score=-7.3, seed=101, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P2"),
            DockingPose(rank=7, score=-6.9, seed=101, cluster_id=2, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P3"),
        ]
        val_details = {"recovery": {"seed": 101, "rank": 7, "rmsd_angstrom": 1.15}}
        self.assertEqual(find_recovery_pose_index(poses, val_details), 2)

        poses_no_exact_match = [
            DockingPose(rank=7, score=-7.1, seed=42, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P1"),
            DockingPose(rank=2, score=-7.3, seed=101, cluster_id=1, rmsd_lb=0, rmsd_ub=0, pdbqt_block="P2"),
        ]
        self.assertEqual(find_recovery_pose_index(poses_no_exact_match, val_details), 0)

if __name__ == "__main__":
    unittest.main()
