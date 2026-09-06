# Models for Protein-Ligand Docking Pipeline
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any

@dataclass
class GridBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float = 20.0
    size_y: float = 20.0
    size_z: float = 20.0

@dataclass
class DockingPose:
    rank: int
    score: float            # kcal/mol
    seed: int
    cluster_id: int
    rmsd_lb: float
    rmsd_ub: float
    rmsd_to_reference: Optional[float] = None
    pdbqt_block: str = ""

@dataclass
class PoseCluster:
    cluster_id: int
    seed_count: int
    supporting_seeds: List[int]
    top_score: float
    median_score: float
    representative_pose: DockingPose
    representative_pose_source: Dict[str, int] = field(default_factory=dict)

@dataclass
class DockingResult:
    run_id: str
    status: str             # "validation passed" | "validation failed" | "exploratory" | "run failed" | "validation unavailable"
    created_at_utc: str
    target_pdb: str
    ligand_name: str
    grid: GridBox
    poses: List[DockingPose]
    clusters: List[PoseCluster]
    redocking_rmsd: Optional[float] = None
    validation_status: str = "exploratory"
    prepared_receptor_path: str = ""
    prepared_ligand_path: str = ""
    manifest_path: str = ""
    warnings: List[str] = field(default_factory=list)
