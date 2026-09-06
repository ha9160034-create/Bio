# Storage and Manifest Management
import os
import json
import uuid
import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
VINA_PATH = PROJECT_ROOT / "tools" / "vina" / "vina.exe"
BASE_RUNS_DIR = PROJECT_ROOT / "docking_runs"

def get_vina_path() -> Path:
    return VINA_PATH

def create_run_directory(prefix="run") -> tuple[str, Path]:
    BASE_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
    uid = uuid.uuid4().hex[:6]
    run_id = f"{prefix}_{ts}_{uid}"
    run_dir = BASE_RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_id, run_dir

def write_initial_manifest(run_dir: Path, data: dict) -> Path:
    initial_path = run_dir / "manifest.initial.json"
    with open(initial_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    return initial_path

def finalize_manifest(run_dir: Path, data: dict) -> Path:
    manifest_path = run_dir / "manifest.json"
    temp_path = run_dir / "manifest.json.tmp"
    with open(temp_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    temp_path.replace(manifest_path)
    return manifest_path

def load_manifest(manifest_path: Path) -> dict:
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)
