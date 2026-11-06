"""`import kite` must not make the interpreter complain.

A docstring once carried `\\in` without an `r` prefix, so every import printed
`SyntaxWarning: invalid escape sequence '\\i'`. It is only a docstring, but the
user sees it on every session, and CPython has said this becomes a `SyntaxError`
-- under `-W error::SyntaxWarning` it already is one, which is what this file
turns into a failure.

The check is on the whole package rather than on the one line that was wrong:
any module reachable from `import kite` is covered, and a new bad escape in a
file nobody is looking at fails here instead of surfacing in a user's terminal.

It runs in a subprocess because warnings are raised at compile time and a module
already imported into this pytest process would be served from the cache.
"""

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PACKAGE = REPO / "python"


def _import_under(*flags):
    return subprocess.run(
        [sys.executable, *flags, "-c", "import kite"],
        cwd=PACKAGE, capture_output=True, text=True)


def test_importing_kite_raises_no_syntax_warning():
    done = _import_under("-W", "error::SyntaxWarning")
    assert done.returncode == 0, (
        "importing kite under -W error::SyntaxWarning failed. An invalid escape "
        "in a docstring or string literal is the usual cause; prefix that string "
        f"with r.\n\n{done.stderr}")


def test_importing_kite_is_silent():
    """Nothing on stderr at all -- a warning the user sees is a defect."""
    done = _import_under()
    assert done.returncode == 0, done.stderr
    assert not done.stderr.strip(), (
        f"importing kite wrote to stderr:\n{done.stderr}")
