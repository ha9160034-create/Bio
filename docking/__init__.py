# Docking package
import sys
import types
import re
from pathlib import Path

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

