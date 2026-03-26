# vendor/whatshap_core/whatshap/__init__.py

import sys

# On Linux, other compiled extensions (readselect/priorityqueue) may need C++ symbols
# from core.so. Make core load with RTLD_GLOBAL so its symbols are visible.
try:
    import ctypes
    RTLD_GLOBAL = ctypes.RTLD_GLOBAL
except Exception:
    RTLD_GLOBAL = 0

if hasattr(sys, "getdlopenflags") and RTLD_GLOBAL:
    _old_flags = sys.getdlopenflags()
    sys.setdlopenflags(_old_flags | RTLD_GLOBAL)
    from . import core
    sys.setdlopenflags(_old_flags)
else:
    from . import core

from . import readselect
from . import priorityqueue
