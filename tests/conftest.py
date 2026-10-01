"""Shared pytest setup for the KITE test suite.

Two things are settled here rather than left to the caller's environment.

**Which `kite` is under test.** An editable install of this package may point
at a *different* checkout, in which case `import kite` inside this repository
silently loads that one and the suite passes while testing code nobody is
editing. The repository's own `python/` is put first on `sys.path` so the tests
always exercise the tree they live in.

**The legacy large tests.** `tests/large001_magdos2d` and `tests/large002_magdos3d`
read `dos.dat` at module scope, so pytest cannot import them without a prior
KITEx run; collection fails before any other test runs. They are excluded from
collection here and remain runnable by hand.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(REPO / "python"))

collect_ignore_glob = ["large*/*"]
