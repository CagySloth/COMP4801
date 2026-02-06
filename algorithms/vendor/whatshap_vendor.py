from __future__ import annotations

import importlib
import sys
from pathlib import Path


def import_whatshap_vendor():
    """
    Ensure `import whatshap` resolves to COMP4801/vendor/whatshap_core
    (not a pip-installed whatshap), then return (whatshap, core, readselect).
    """
    repo_root = Path(__file__).resolve().parents[2]  # .../COMP4801
    vendor_root = (repo_root / "vendor" / "whatshap_core").resolve()
    vendor_root_str = str(vendor_root)

    # Make vendor highest priority on sys.path
    if sys.path[0] != vendor_root_str:
        if vendor_root_str in sys.path:
            sys.path.remove(vendor_root_str)
        sys.path.insert(0, vendor_root_str)

    # If whatshap is already imported from elsewhere, purge it.
    m = sys.modules.get("whatshap")
    if m is not None:
        mod_file = str(getattr(m, "__file__", "") or "")
        if mod_file and not str(Path(mod_file).resolve()).startswith(vendor_root_str):
            for k in list(sys.modules.keys()):
                if k == "whatshap" or k.startswith("whatshap."):
                    del sys.modules[k]

    wh = importlib.import_module("whatshap")
    from whatshap import core, readselect  # type: ignore

    # Hard check: must be vendor
    wh_file = str(Path(wh.__file__).resolve())
    if not wh_file.startswith(vendor_root_str):
        raise RuntimeError(
            "Loaded non-vendored `whatshap`.\n"
            f"Expected under: {vendor_root_str}\n"
            f"Got: {wh_file}\n\n"
            "Fix: ensure vendor path is first on sys.path and/or uninstall pip whatshap."
        )

    return wh, core, readselect
