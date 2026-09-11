# Docking package
import sys
import os
import types
import re
from pathlib import Path

def _configure_openbabel_paths():
    """Detect OpenBabel plugins and data directories dynamically for Linux, Windows, and Cloud."""
    plugin_dirs = [
        "/usr/lib/x86_64-linux-gnu/openbabel/3.1.1",
        "/usr/lib/x86_64-linux-gnu/openbabel/3.1.0",
        "/usr/lib/openbabel/3.1.1",
        "/usr/lib/openbabel/3.1.0",
        "/usr/local/lib/openbabel/3.1.1",
    ]
    data_dirs = [
        "/usr/share/openbabel/3.1.1",
        "/usr/share/openbabel/3.1.0",
        "/usr/share/openbabel",
    ]
    try:
        import openbabel
        ob_root = Path(openbabel.__file__).parent
        for p in ob_root.rglob("*"):
            if p.is_file():
                if p.suffix in [".so", ".obf", ".dylib"] or "format" in p.name:
                    plugin_dirs.insert(0, str(p.parent))
                elif p.suffix == ".txt" and ("aromatic" in p.name or "element" in p.name or "types" in p.name):
                    data_dirs.insert(0, str(p.parent))
    except Exception:
        pass

    for cand in plugin_dirs:
        if os.path.isdir(cand):
            os.environ["BABEL_LIBDIR"] = cand
            break

    for cand in data_dirs:
        if os.path.isdir(cand):
            os.environ["BABEL_DATADIR"] = cand
            break

_configure_openbabel_paths()

# Self-healing patch for upstream OpenBabel pybel.py unpacking bug on Linux / Python 3.14+
def _apply_openbabel_patch():
    if "openbabel.pybel" in sys.modules:
        return
    try:
        import openbabel
        pybel_path = Path(openbabel.__file__).parent / "pybel.py"
        if not pybel_path.exists():
            return
        text = pybel_path.read_text(encoding="utf-8", errors="ignore")
        pat = r"broken\s*=\s*\[\(x,\s*y\.strip\(\)\)\s*for\s*x,\s*y\s*in\s*broken\]"
        fixed = "broken = [(p[0], p[1].strip()) if len(p) >= 2 else (p[0], p[0]) for p in broken]"
        if re.search(pat, text):
            patched_text = re.sub(pat, fixed, text, count=1)
            try:
                pybel_path.write_text(patched_text, encoding="utf-8")
            except Exception:
                pass
            mod = types.ModuleType("openbabel.pybel")
            mod.__file__ = str(pybel_path)
            mod.__package__ = "openbabel"
            exec(patched_text, mod.__dict__)
            sys.modules["openbabel.pybel"] = mod
            openbabel.pybel = mod
    except Exception:
        pass

_apply_openbabel_patch()
