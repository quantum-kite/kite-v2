"""A KITE-tools flag must reach the output, not just the banner.

Every reader is constructed the same way: `set_default_parameters()`, then
`fetch_parameters()` reading the .h5, then `override_parameters()` reading the
command line. The command line is applied last, so it wins -- that is the
documented post-processing workflow, and the ordering is what guarantees it.

This file asserts the guarantee twice. Structurally, every reader must call
`fetch_parameters()` before `override_parameters()`; a reader that reversed them
would silently hand the stored value back whatever the user typed. Behaviourally,
a flag passed to a working entry point must change its output file.

The investigation behind this (BUG-003) was opened on an observation that
`--CondOpt2` wrote byte-identical output with and without `-T` and `-O`. That
turned out not to be a provenance problem: the second-order path returns NaN for
its stored default temperature, and one NaN table is byte-identical to another.
With a temperature that does not overflow the Fermi factor, every documented
`--CondOpt2` flag changes the result. The NaN itself belongs to T-062.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TOOLS_SRC = REPO / "tools" / "Src"
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

_CALL = re.compile(
    r"\b(set_default_parameters|fetch_parameters|override_parameters)\s*\(\s*\)")


def _readers():
    for path in sorted(TOOLS_SRC.rglob("*.cpp")):
        seq = _CALL.findall(path.read_text())
        if "override_parameters" in seq and "fetch_parameters" in seq:
            yield path.relative_to(REPO).as_posix(), seq


def test_every_reader_was_found():
    """A pattern that stopped matching would make the next test vacuous."""
    found = list(_readers())
    assert len(found) >= 13, (
        f"only {len(found)} readers matched; the construction sequence is no "
        "longer written the way this file looks for it")


@pytest.mark.parametrize("path,seq", list(_readers()))
def test_the_command_line_is_applied_after_the_stored_values(path, seq):
    """Reversing these two would make every documented flag a no-op."""
    first = seq.index("fetch_parameters")
    override = seq.index("override_parameters")
    assert first < override, (
        f"{path} calls override_parameters() before fetch_parameters(), so the "
        "stored value overwrites whatever the user passed on the command line")


# --- and the behaviour that ordering is supposed to produce -------------------

_BUILDER = r'''
import sys, numpy as np, kite
from kite import custom
from kite import lattice as latt

a, acc = 0.24595, 0.142
lat = latt.Lattice(a1=a*np.array([1.0, 0.0]), a2=a*np.array([0.5, 0.5*np.sqrt(3)]))
lat.add_sublattices(('A', [0, -acc/2], 0.0), ('B', [0, acc/2], 0.0))
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([1, -1], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0))

cfg = kite.Configuration(divisions=[2, 2], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-4, 4])
calc = kite.Calculation(cfg)
calc.dos(num_points=500, num_moments=64, num_random=1, num_disorder=1)
calc.conductivity_dc(num_points=200, num_moments=64, num_random=1,
                     num_disorder=1, direction='xx', temperature=0.01)
A = custom.Vertex(64, [[1.0, "vx"]])
B = custom.Vertex(64, [[1.0, "vx"]])
calc.custom_two(stream_=[A, B], num_random_=1, num_disorder_=1,
                num_points_=200, temperature_=0.01)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    if not (KITEX.exists() and KITE_TOOLS.exists()):
        pytest.skip("needs a built KITEx and KITE-tools")
    d = tmp_path_factory.mktemp("provenance")
    subprocess.run([sys.executable, "-c", _BUILDER, "p.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "p.h5"], cwd=d, check=True, capture_output=True)
    return d


def _output(d, flag, *extra):
    name = f"out_{abs(hash((flag, extra)))}.dat"
    subprocess.run([str(KITE_TOOLS), "p.h5", flag, *extra, "-N", name],
                   cwd=d, check=True, capture_output=True, text=True)
    return (d / name).read_text()


@pytest.mark.parametrize("flag,extra", [
    ("--DOS", ("-E", "123")),
    ("--DOS", ("-M", "32")),
    ("--CondDC", ("-T", "0.5")),
    ("--CondDC", ("-E", "301")),
    ("--CustomTwo", ("-T", "0.5")),
    ("--CustomTwo", ("-E", "301")),
])
def test_the_flag_changes_the_output(run, flag, extra):
    """Not just the banner: the numbers the user takes away have to move."""
    base = _output(run, flag)
    changed = _output(run, flag, *extra)
    assert changed != base, (
        f"{flag} {' '.join(extra)} produced output identical to the run without "
        "it, so the stored value won over the command line")
