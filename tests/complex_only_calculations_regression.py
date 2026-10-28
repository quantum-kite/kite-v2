"""`ldos_map` and `spectral_map` must refuse a real Hamiltonian, not no-op.

Both drivers seed a random-PHASE vector, so KITEx compiles their bodies only for
a complex scalar type: `SimulationLDoS.cpp` and `SimulationSpectral.cpp` wrap
everything in `if constexpr (is_tt<std::complex, T>::value)`. Built with
`is_complex=False` the calculation prints "Calculating LDoS." and exits zero
having written no Map dataset at all -- the same silent-success failure as
BUG-009, where the user cannot tell a missing result from a real one.

`config_system` therefore refuses at export, before a run is ever started.

Two halves, and the second is the one that rots: the refusal must fire, and it
must NOT fire for a complex Hamiltonian, or it would block the only supported
way of asking for these maps.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "Src" / "Simulation"

#: driver file -> the Calculation method whose request it serves
COMPLEX_ONLY = {
    "SimulationLDoS.cpp": "ldos_map",
    "SimulationSpectral.cpp": "spectral_map",
}

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt
lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
lat.add_sublattices(['A', [0, 0], 0.0], ['B', [0.5, 0.5], -0.0])
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([-1, 0], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0), ([-1, -1], 'A', 'B', -1.0))
cfg = kite.Configuration(divisions=[1, 1], length=[16, 16],
                         boundaries=['periodic', 'periodic'],
                         is_complex=(sys.argv[2] == "complex"),
                         precision=1, spectrum_range=[-5, 5])
calc = kite.Calculation(cfg)
getattr(calc, sys.argv[3])(energy_=0.2, sigma_=0.15, vectors_=8)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
print("EXPORTED")
'''


def _export(tmp_path, scalar, method):
    return subprocess.run(
        [sys.executable, "-c", _BUILDER, str(tmp_path / "x.h5"), scalar, method],
        cwd=tmp_path, capture_output=True, text=True)


@pytest.mark.parametrize("driver,method", sorted(COMPLEX_ONLY.items()))
def test_the_driver_is_still_complex_only(driver, method):
    """If the guard goes, the refusal is wrong and must go with it."""
    text = (SRC / driver).read_text()
    assert re.search(r"if constexpr \(is_tt<std::complex, T>::value\)", text), (
        f"{driver} is no longer compiled for complex scalars only, so "
        f"config_system's refusal of {method} on a real Hamiltonian is now "
        "rejecting something that would work")


@pytest.mark.parametrize("method", sorted(COMPLEX_ONLY.values()))
def test_a_real_hamiltonian_is_refused(tmp_path, method):
    result = _export(tmp_path, "real", method)
    assert result.returncode != 0, (
        f"{method} exported with is_complex=False; KITEx would report success "
        "and write no map")
    assert "complex Hamiltonian" in result.stdout + result.stderr, (
        f"{method} was refused without saying a complex Hamiltonian is needed")
    assert "is_complex=True" in result.stdout + result.stderr, (
        "the refusal does not tell the user what to change")


@pytest.mark.parametrize("method", sorted(COMPLEX_ONLY.values()))
def test_a_complex_hamiltonian_is_allowed(tmp_path, method):
    result = _export(tmp_path, "complex", method)
    assert result.returncode == 0 and "EXPORTED" in result.stdout, (
        f"{method} no longer exports with is_complex=True:\n"
        + result.stdout + result.stderr)
