"""
Legacy matrix simulator compatibility shim.

Real implementation moved to: dataset.legacy.simulate

This file exists so older commands still work:
  python -m dataset.simulate
  simulate (console script entrypoint)
"""
from dataset.legacy.simulate import *  # re-export legacy API
from dataset.legacy.simulate import main  # explicit for entrypoints

if __name__ == "__main__":
    main()