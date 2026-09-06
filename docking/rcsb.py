"""
RCSB PDB and Biological Assembly coordinate downloading service.
Independent module to avoid circular imports between b.py and docking modules.
"""
import gzip
import requests

def fetch_deposited_pdb_text(pdb_id: str, timeout: int = 60):
    """Download deposited coordinates (asymmetric unit) from RCSB."""
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.text
    except Exception:
        return None

def fetch_biological_assembly_text(pdb_id: str, assembly_id: int = 1, timeout: int = 60):
    """Download official RCSB biological assembly (.pdb1.gz)."""
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb{assembly_id}.gz"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return gzip.decompress(response.content).decode("utf-8")
    except Exception:
        return None
