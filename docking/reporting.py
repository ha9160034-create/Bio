# Reporting and Manifest Generation
import hashlib
from pathlib import Path
from dataclasses import asdict
from docking.models import DockingResult
from docking.storage import finalize_manifest

DISCLAIMER_TEXT = (
    "تنبيه علمي: هذه النتائج تمثل فرضية ارتباط حاسوبية (In-Silico Docking Hypothesis) "
    "تهدف لترتيب الوضعيات والمركبات ضمن هذا البروتوكول المحدد. "
    "لا تُعد هذه الدرجات قياساً معملياً لطاقة الارتباط الحرة أو دليلاً على الفاعلية السريرية أو الجرعة العلاجية."
)

def compute_file_sha256(filepath: Path) -> str:
    if not filepath.exists() or not filepath.is_file():
        return ""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def build_and_save_manifest(
    run_dir: Path,
    run_id: str,
    status: str,
    target_meta: dict,
    binding_site_meta: dict,
    ligand_meta: dict,
    prep_meta: dict,
    engine_meta: dict,
    grid_meta: dict,
    results_meta: dict,
    warnings: list[str] = None
) -> Path:
    run_dir = Path(run_dir)

    # Valid status values per repair guidelines
    valid_statuses = [
        "validation passed",
        "validation passed (Top-10 recovery)",
        "validation failed",
        "validation unavailable",
        "exploratory—no reference validation",
        "exploratory—insufficient independent seeds",
        "run failed"
    ]
    if status not in valid_statuses:
        status = "run failed" if "fail" in status.lower() else "exploratory—no reference validation"

    created_at_utc = target_meta.get("created_at_utc") or ""

    manifest_data = {
        "schema_version": 1,
        "run_id": run_id,
        "status": status,
        "created_at_utc": created_at_utc,
        "target": target_meta,
        "binding_site": binding_site_meta,
        "ligand": ligand_meta,
        "receptor_preparation": prep_meta,
        "engine": engine_meta,
        "grid": grid_meta,
        "inputs": {
            "sha256": {
                p.name: compute_file_sha256(p)
                for p in sorted(run_dir.iterdir())
                if p.is_file() and not p.name.startswith("manifest")
            }
        },
        "results": results_meta,
        "warnings": warnings or [],
        "scientific_disclaimer": DISCLAIMER_TEXT
    }
    return finalize_manifest(run_dir, manifest_data)
