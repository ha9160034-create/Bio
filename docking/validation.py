# Docking validation, atom mapping, and RMSD calculations
import numpy as np
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdFMCS
from docking.models import DockingPose, PoseCluster

def parse_pdbqt_coords(pdbqt_block: str) -> np.ndarray:
    """Extract coordinates of heavy atoms from a PDBQT block."""
    coords = []
    for line in pdbqt_block.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_type = line[77:].strip() if len(line) >= 78 else ""
            if not atom_type.startswith("H"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass
    return np.array(coords, dtype=np.float32)

def calculate_heavy_atom_rmsd(ref_coords: np.ndarray, pred_coords: np.ndarray) -> float:
    """Calculate heavy-atom RMSD requiring strictly identical atom counts."""
    if len(ref_coords) != len(pred_coords) or len(ref_coords) == 0:
        raise ValueError(
            f"Cannot calculate RMSD: unequal atom count (ref={len(ref_coords)}, pred={len(pred_coords)})."
        )
    diff = ref_coords - pred_coords
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

def calculate_mapped_redocking_rmsd(ref_ligand_pdb_text: str, pose_pdbqt_block: str) -> tuple[float, dict]:
    """
    Calculate heavy-atom RMSD between reference crystal ligand and a docked pose
    using Maximum Common Substructure (MCS) topology-invariant atom mapping.
    Evaluates every pair of reference and pose mappings to handle chemical symmetry correctly.
    Refuses incomplete mappings or unequal heavy-atom counts.
    Do not fit or rotate the ligand before RMSD; the docked and crystal ligand share the receptor coordinate frame.
    """
    if not ref_ligand_pdb_text.strip() or not pose_pdbqt_block.strip():
        raise ValueError("Missing reference ligand or pose coordinates.")

    # 1. Parse reference ligand into RDKit molecule
    mol_ref = Chem.MolFromPDBBlock(ref_ligand_pdb_text, sanitize=False)
    if mol_ref is None:
        ob_ref = pybel.readstring("pdb", ref_ligand_pdb_text)
        mol_ref = Chem.MolFromMolBlock(ob_ref.write("mol"), sanitize=False)

    if mol_ref is None:
        raise ValueError("Could not parse reference ligand coordinates into a chemical molecule.")

    # 2. Convert pose PDBQT block into RDKit molecule via OpenBabel to preserve coordinates
    ob_pose = pybel.readstring("pdbqt", pose_pdbqt_block)
    mol_pose = Chem.MolFromPDBBlock(ob_pose.write("pdb"), sanitize=False)
    if mol_pose is None:
        mol_pose = Chem.MolFromMolBlock(ob_pose.write("mol"), sanitize=False)

    if mol_pose is None:
        raise ValueError("Could not parse docked pose coordinates into a chemical molecule.")

    ref_heavy_count = mol_ref.GetNumHeavyAtoms()
    pose_heavy_count = mol_pose.GetNumHeavyAtoms()

    if ref_heavy_count != pose_heavy_count or ref_heavy_count == 0:
        raise ValueError(
            f"Validation unavailable: reference ligand has {ref_heavy_count} heavy atoms, "
            f"docked ligand has {pose_heavy_count} heavy atoms. Molecules must match exactly."
        )

    # 3. Establish documented atom mapping via MCS
    mcs = rdFMCS.FindMCS(
        [mol_ref, mol_pose],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        matchValences=False,
        ringMatchesRingOnly=False,
        completeRingsOnly=False
    )

    if mcs.numAtoms != ref_heavy_count:
        raise ValueError(
            f"Validation unavailable: atom mapping could only match {mcs.numAtoms}/{ref_heavy_count} heavy atoms. "
            "Chemical identity mismatch between reference and docked ligand."
        )

    query_mol = Chem.MolFromSmarts(mcs.smartsString)
    match_ref_all = mol_ref.GetSubstructMatches(query_mol, uniquify=False)
    match_pose_all = mol_pose.GetSubstructMatches(query_mol, uniquify=False)

    if not match_ref_all or not match_pose_all:
        raise ValueError("Validation unavailable: could not map atoms from MCS query.")

    conf_ref = mol_ref.GetConformer()
    conf_pose = mol_pose.GetConformer()

    # 4. Iterate through every pair of reference and pose mappings to find minimum RMSD
    best_rmsd = 999.0
    best_mapping = None
    mappings_evaluated = 0

    for m_ref in match_ref_all:
        for m_pose in match_pose_all:
            mappings_evaluated += 1
            diffs = []
            mapping_pairs = []
            for idx_q in range(len(m_ref)):
                idx_r = m_ref[idx_q]
                idx_p = m_pose[idx_q]
                p_r = conf_ref.GetAtomPosition(idx_r)
                p_p = conf_pose.GetAtomPosition(idx_p)
                diffs.append([p_r.x - p_p.x, p_r.y - p_p.y, p_r.z - p_p.z])
                mapping_pairs.append({"ref_atom_idx": int(idx_r), "pose_atom_idx": int(idx_p)})

            diffs_arr = np.array(diffs, dtype=np.float32)
            rmsd = float(np.sqrt(np.mean(np.sum(diffs_arr**2, axis=1))))
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_mapping = mapping_pairs

    details = {
        "mapping_method": "rdkit_mcs_topology_invariant_exhaustive_pairs",
        "reference_heavy_atom_count": ref_heavy_count,
        "pose_heavy_atom_count": pose_heavy_count,
        "matched_atom_count": mcs.numAtoms,
        "mappings_evaluated": mappings_evaluated,
        "selected_mapping_pairs": best_mapping,
        "final_rmsd": round(best_rmsd, 3)
    }
    return round(best_rmsd, 3), details

def cluster_poses_across_seeds(poses: list[DockingPose], rmsd_threshold: float = 2.0) -> list[PoseCluster]:
    """
    Cluster poses across independent seeds using heavy-atom RMSD.
    Returns clustered groups with seed counts, median scores, and representative poses.
    """
    if not poses:
        return []

    # Sort poses by binding score (ascending, best/lowest first)
    sorted_poses = sorted(poses, key=lambda p: p.score)
    clusters: list[list[DockingPose]] = []

    for pose in sorted_poses:
        p_coords = parse_pdbqt_coords(pose.pdbqt_block)
        assigned = False
        for cl in clusters:
            rep_coords = parse_pdbqt_coords(cl[0].pdbqt_block)
            if len(p_coords) == len(rep_coords) and len(p_coords) > 0:
                rmsd = float(np.sqrt(np.mean(np.sum((p_coords - rep_coords)**2, axis=1))))
                if rmsd <= rmsd_threshold:
                    cl.append(pose)
                    assigned = True
                    break
        if not assigned:
            clusters.append([pose])

    result_clusters = []
    for idx, cl in enumerate(clusters, start=1):
        supporting_seeds = sorted(list(set(p.seed for p in cl)))
        scores = [p.score for p in cl]
        median_score = float(np.median(scores))
        top_score = float(min(scores))
        for p in cl:
            p.cluster_id = idx
        result_clusters.append(PoseCluster(
            cluster_id=idx,
            seed_count=len(supporting_seeds),
            supporting_seeds=supporting_seeds,
            top_score=round(top_score, 2),
            median_score=round(median_score, 2),
            representative_pose=cl[0],
            representative_pose_source={"seed": cl[0].seed, "rank": cl[0].rank}
        ))

    return result_clusters
