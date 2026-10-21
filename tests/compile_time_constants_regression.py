"""Python's copy of `MEMORY` must agree with the one KITEx was built with.

`MEMORY` is a compile-time constant in `Src/Generic.hpp`, and it is
`#ifndef`-guarded, so a build can override it with `-DMEMORY=...`. The Python
exporter cannot see that flag: it carries its own `KITE_MEMORY` and rejects a
`num_moments` that is not a multiple of it.

That check is the only real protection. The four C++ guards -- `Gamma2D.cpp`,
`Gamma3D.cpp`, `Custom/SimulationRankTwo.cpp`, `Custom/SimulationRankTwoLocal.cpp`
-- test `% 2`, i.e. evenness, while the loops they guard stride by `MEMORY`. That
gap is F-01: `num_moments=520` is even, is not a multiple of 16, passes every
guard and overruns the heap.

So if the two copies of the constant ever diverge, the export check silently
stops matching the loop it exists to protect. This asserts they agree, and that
the Python guard is still wired to the constant rather than to a literal.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
GENERIC = REPO / "Src" / "Generic.hpp"
API = REPO / "python" / "kite" / "__init__.py"


def _cpp_default(name):
    """The value `Generic.hpp` compiles with when no -D flag overrides it."""
    for line in GENERIC.read_text().splitlines():
        m = re.match(rf"\s*#define\s+{name}\s+(\d+)\s*$", line)
        if m:
            return int(m.group(1))
    return None


def test_python_and_cpp_agree_on_the_block_size():
    cpp = _cpp_default("MEMORY")
    assert cpp is not None, "Src/Generic.hpp no longer defines MEMORY"
    py = re.search(r"^KITE_MEMORY\s*=\s*(\d+)", API.read_text(), re.M)
    assert py is not None, "python/kite no longer defines KITE_MEMORY"
    assert int(py.group(1)) == cpp, (
        f"KITE_MEMORY is {py.group(1)} but Src/Generic.hpp compiles MEMORY as "
        f"{cpp}. The export check and the C++ loop stride no longer agree, which "
        "is the F-01 failure mode: a num_moments that passes export overruns the "
        "block loop")


def test_memory_is_still_overridable_at_build_time():
    """If the guard ever goes, the hazard changes shape and §3.5 is stale."""
    text = GENERIC.read_text()
    assert re.search(r"#ifndef\s+MEMORY", text), (
        "MEMORY is no longer #ifndef-guarded; ARCHITECTURE §3.5 describes the "
        "build-flag hazard that guard creates")


def test_the_export_check_is_tied_to_the_constant_not_a_literal():
    text = API.read_text()
    assert re.search(r"num_moments\s*%\s*KITE_MEMORY", text), (
        "the block-stepping export check no longer reads KITE_MEMORY; if it "
        "hardcodes a number it cannot be kept in step with Src/Generic.hpp")
