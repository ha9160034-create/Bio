# AutoDock Vina Engine Adapter
import os
import subprocess
from pathlib import Path
from docking.models import GridBox, DockingPose
from docking.storage import get_vina_path

def get_vina_version() -> str:
    vina_exec = get_vina_path()
    if not vina_exec.exists():
        return "not_installed"
    try:
        res = subprocess.run([str(vina_exec), "--version"], capture_output=True, text=True, check=True)
        return res.stdout.strip() or res.stderr.strip()
    except Exception as e:
        return f"unknown ({e})"

def run_vina_single_seed(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    grid: GridBox,
    output_dir: Path,
    seed: int = 42,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0
) -> tuple[list[DockingPose], dict]:
    vina_exec = get_vina_path()
    if not vina_exec.exists():
        raise FileNotFoundError(f"Vina binary not found at: {vina_exec}")

    # Validate numeric parameters
    if exhaustiveness < 1:
        raise ValueError("Exhaustiveness must be >= 1")
    if num_modes < 1:
        raise ValueError("num_modes must be >= 1")
    if grid.size_x <= 0 or grid.size_y <= 0 or grid.size_z <= 0:
        raise ValueError("Grid box dimensions must be strictly positive.")

    out_pdbqt = output_dir / f"out_seed_{seed}.pdbqt"
    log_file = output_dir / f"vina_seed_{seed}.log"

    cmd = [
        str(vina_exec),
        "--receptor", str(receptor_pdbqt.resolve()),
        "--ligand", str(ligand_pdbqt.resolve()),
        "--center_x", str(grid.center_x),
        "--center_y", str(grid.center_y),
        "--center_z", str(grid.center_z),
        "--size_x", str(grid.size_x),
        "--size_y", str(grid.size_y),
        "--size_z", str(grid.size_z),
        "--out", str(out_pdbqt.resolve()),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--seed", str(seed)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
    with open(log_file, "w", encoding="utf-8") as f:
        f.write(result.stdout + "\n" + result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"Vina failed with exit code {result.returncode}: {result.stderr}")

    if not out_pdbqt.exists():
        raise FileNotFoundError(f"Expected Vina output file not found: {out_pdbqt}")

    with open(out_pdbqt, "r", encoding="utf-8") as f:
        content = f.read()

    models = content.split("MODEL ")
    if len(models) <= 1:
        raise ValueError(f"No models found in Vina output for seed {seed}.")

    poses: list[DockingPose] = []
    for rank_idx, m in enumerate(models[1:], start=1):
        lines = m.splitlines()
        found_score = False
        score = 0.0
        rmsd_lb = 0.0
        rmsd_ub = 0.0
        for line in lines:
            if "REMARK VINA RESULT:" in line:
                parts = line.split()
                if len(parts) >= 6:
                    score = float(parts[3])
                    rmsd_lb = float(parts[4])
                    rmsd_ub = float(parts[5])
                    found_score = True
                    break

        if not found_score:
            raise ValueError(
                f"Corrupted model in Vina output: missing REMARK VINA RESULT tag for model {rank_idx} (seed {seed})."
            )

        pose_block = "MODEL " + m
        poses.append(DockingPose(
            rank=rank_idx,
            score=score,
            seed=seed,
            cluster_id=0,
            rmsd_lb=rmsd_lb,
            rmsd_ub=rmsd_ub,
            pdbqt_block=pose_block
        ))

    run_meta = {
        "seed": seed,
        "command": cmd,
        "vina_version": get_vina_version(),
        "exit_code": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "out_pdbqt": str(out_pdbqt.name),
        "log_file": str(log_file.name),
        "pose_count": len(poses),
        "best_score": min((p.score for p in poses), default=None)
    }

    return poses, run_meta

def run_vina_multi_seeds(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    grid: GridBox,
    output_dir: Path,
    seeds: list[int] = [42, 101, 2024],
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0
) -> tuple[list[DockingPose], list[dict]]:
    all_poses: list[DockingPose] = []
    runs_meta: list[dict] = []

    for s in seeds:
        seed_poses, meta = run_vina_single_seed(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            grid=grid,
            output_dir=output_dir,
            seed=s,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            energy_range=energy_range
        )
        all_poses.extend(seed_poses)
        runs_meta.append(meta)

    return all_poses, runs_meta
