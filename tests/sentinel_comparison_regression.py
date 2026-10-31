"""A KITE-tools override must compare against the sentinel its field was initialised to.

Each `--Flag -x value` option is held in a `variables.<Field>`, initialised in
`parse_input.cpp` to a value no user would pick, and each reader's
`override_parameters()` decides "the user supplied this" by comparing against
that same value. The two constants are written in different files, several
hundred lines apart, and nothing ties them together.

`conductivity_2order.cpp` compared `variables.CondOpt2_ratio` against `8.8888`
where the field is initialised to `-8.8888`, so the test was true whenever `-R`
was absent and the sentinel itself was assigned as the frequency ratio. Every
`--CondOpt2` run without `-R` evaluated the response at `w2 = -8.8888*w1`
instead of the intended `w2 = w1` (BUG-011).

The mismatch is mechanical, so the check is mechanical: it reads both files and
compares every pair. That is a stronger guarantee than a test of the one option
that was wrong, because the next one to drift will not be this one.

The second half of this file covers why that value went unnoticed. The banner
is the only place a sentinel ever becomes visible, and in this reader the
"(default)" marker was printed on every line regardless -- neither
`fetch_parameters()` nor `override_parameters()` cleared a `default_` flag, so
`-8.8888` was displayed exactly as a deliberate default would be.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TOOLS_SRC = REPO / "tools" / "Src"
PARSE_INPUT = TOOLS_SRC / "Tools" / "parse_input.cpp"
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

#: `Field = value;` in parse_input's parse_* functions.
_INIT = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(-?[\d.]+)\s*;", re.M)
#: `variables.Field != value` in a reader's override_parameters().
_COMPARE = re.compile(r"variables\.([A-Za-z0-9_]+)\s*!=\s*(-?[\d.]+)")


def _initialisers():
    out = {}
    for m in _INIT.finditer(PARSE_INPUT.read_text()):
        out.setdefault(m.group(1), m.group(2))
    return out


def _comparisons():
    init = _initialisers()
    for path in sorted(TOOLS_SRC.rglob("*.cpp")):
        for m in _COMPARE.finditer(path.read_text()):
            if m.group(1) in init:
                yield path.relative_to(REPO).as_posix(), m.group(1), m.group(2)


def test_there_are_sentinel_comparisons_to_check():
    """A regex that silently matches nothing would pass every assertion below."""
    found = list(_comparisons())
    assert len(found) > 40, (
        f"only {len(found)} sentinel comparisons found; the pattern no longer "
        "matches how parse_input or the readers are written, so this file is "
        "no longer checking anything")


@pytest.mark.parametrize("path,field,compared", sorted(_comparisons()))
def test_the_comparison_matches_the_initialiser(path, field, compared):
    expected = _initialisers()[field]
    assert float(compared) == float(expected), (
        f"{path} tests `variables.{field} != {compared}` but parse_input.cpp "
        f"initialises {field} to {expected}. The condition is then true when "
        "the user supplied nothing, and the sentinel is used as if it were a "
        "value they chose")


# --- and the behaviour that mismatch produced ---------------------------------

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

a, acc = 0.24595, 0.142
lat = latt.Lattice(a1=a*np.array([1.0, 0.0]), a2=a*np.array([0.5, 0.5*np.sqrt(3)]))
lat.add_sublattices(('A', [0, -acc/2], 0.0), ('B', [0, acc/2], 0.0))
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([1, -1], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0))

cfg = kite.Configuration(divisions=[2, 2], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=False,
                         precision=1, spectrum_range=[-4, 4])
calc = kite.Calculation(cfg)
calc.conductivity_optical_nonlinear(num_points=100, num_moments=32, num_random=1,
                                    num_disorder=1, direction='xxx', temperature=0.01)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''

_RATIO = re.compile(r"^\s*Ratio:\s*(\S+)", re.MULTILINE)


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    if not (KITEX.exists() and KITE_TOOLS.exists()):
        pytest.skip("needs a built KITEx and KITE-tools")
    d = tmp_path_factory.mktemp("sentinel")
    subprocess.run([sys.executable, "-c", _BUILDER, "n.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "n.h5"], cwd=d, check=True, capture_output=True)
    return d


def _ratio(d, *extra):
    out = subprocess.run([str(KITE_TOOLS), "n.h5", *extra], cwd=d,
                         check=True, capture_output=True, text=True).stdout
    found = _RATIO.findall(out)
    assert len(found) == 1, f"expected one Ratio line, got {found}"
    return float(found[0])


def test_omitting_the_ratio_flag_keeps_the_built_in_default(run):
    """set_default_parameters sets ratio = 1.0; nothing should overwrite it."""
    assert _ratio(run) == pytest.approx(1.0), (
        "--CondOpt2 without -R no longer uses the built-in ratio of 1")


@pytest.mark.parametrize("value", [2.0, -1.0])
def test_the_ratio_flag_is_still_honoured(run, value):
    """-R -1 selects the photocurrent branch, so it must survive the fix."""
    assert _ratio(run, "--CondOpt2", "-R", str(value)) == pytest.approx(value)


# --- and the marker that made it look deliberate ------------------------------

_BANNER_LINE = re.compile(r"^\s{3}([A-Za-z][^:]*):\s*(\S+)(\s*\(default\))?\s*$",
                          re.MULTILINE)


def _banner(d, *extra):
    out = subprocess.run([str(KITE_TOOLS), "n.h5", *extra], cwd=d,
                         check=True, capture_output=True, text=True).stdout
    block = out[out.index("second-order optical conductivity"):]
    return {name.strip(): bool(flag) for name, _, flag in _BANNER_LINE.findall(block)}


def test_a_value_read_from_the_file_is_not_marked_default(run):
    """The script asked for these; the built-in defaults are 0.001 and 128."""
    marked = _banner(run)
    for name in ("Temperature", "Number of frequencies"):
        assert marked[name] is False, (
            f"'{name}' came from the .h5 file but is marked (default)")


def test_a_value_the_user_did_not_set_is_still_marked_default(run):
    """The marker has to stay meaningful in the other direction too."""
    marked = _banner(run)
    assert marked["Ratio"] is True, (
        "Ratio was not overridden, so it must be marked (default) -- otherwise "
        "the marker says nothing")
    assert marked["Broadening"] is True


@pytest.mark.parametrize("flag,value,name", [
    ("-T", "0.05", "Temperature"),
    ("-R", "2", "Ratio"),
    ("-E", "201", "Num integration points"),
    ("-F", "0.3", "Fermi energy"),
])
def test_an_overridden_value_is_not_marked_default(run, flag, value, name):
    """Every override_parameters branch must clear its own flag."""
    marked = _banner(run, "--CondOpt2", flag, value)
    assert marked[name] is False, (
        f"'{name}' was set with {flag} {value} but is still marked (default)")
