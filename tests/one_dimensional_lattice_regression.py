"""A 1D lattice must be refused at export, with a workaround that works.

KITEx specializes `KPM_Vector` for D=2 and D=3 only. D=1 resolves to the primary
template, whose methods are empty stubs, so a 1D system used to export cleanly,
run to completion, print "Done." and return identically zero moments -- including
mu_0, which is 1 by normalisation alone. Nothing at any stage said 1D was
unsupported.

`config_system` now stops instead, and the message suggests building the chain as
a 2D lattice with no hopping in the second direction. The chains decouple, so one
run computes all of them and the stochastic trace averages over every chain.

The third test is the important one. A message that gives advice is a promise, and
this checks the promise: the suggested construction is run and compared against
the exact 1D tight-binding density of states, rho(E) = 1/(pi*sqrt(4t^2 - E^2)).

That comparison also has a trap worth keeping in the file, because it caught me:
the chains must be long and `num_moments` modest. A finite chain has a discrete
spectrum, and if the KPM resolution is finer than the level spacing you resolve
individual levels rather than the band. At length 512 with 2048 moments the result
is wildly oscillatory and looks like a bug; at length 4096 with 256 moments it
reproduces the analytic curve to a few percent.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
TOOLS = REPO / "build" / "KITE-tools"

_ONE_D = '''
import kite
from kite import lattice as latt
lattice = latt.Lattice(a1=[1.0])
lattice.add_sublattices(['A', [0.0], 0.0])
lattice.add_hoppings([[-1], 'A', 'A', -1.0])
configuration = kite.Configuration(
    divisions=[1], length=[64], boundaries=["periodic"],
    is_complex=False, precision=1, spectrum_range=[-2.2, 2.2])
calculation = kite.Calculation(configuration)
calculation.dos(num_points=200, num_moments=64, num_random=1)
kite.config_system(lattice, configuration, calculation, filename="d1.h5")
'''

_CHAINS = '''
import kite
from kite import lattice as latt
lattice = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lattice.add_sublattices(['A', [0.0, 0.0], 0.0])
lattice.add_hoppings([[-1, 0], 'A', 'A', -1.0])          # x only: chains decouple
configuration = kite.Configuration(
    divisions=[1, 1], length=[4096, 8], boundaries=["periodic"] * 2,
    is_complex=False, precision=1, spectrum_range=[-2.2, 2.2])
calculation = kite.Calculation(configuration)
calculation.dos(num_points=2000, num_moments=256, num_random=1, num_disorder=1)
kite.config_system(lattice, configuration, calculation, filename="chains.h5")
'''


def _run_export(script, tmp_path):
    return subprocess.run([sys.executable, "-c", script], cwd=str(tmp_path),
                          capture_output=True, text=True)


def test_a_one_dimensional_lattice_is_refused(tmp_path):
    result = _run_export(_ONE_D, tmp_path)
    assert result.returncode != 0, (
        "a 1D lattice exported without complaint; it would run and return zeros")
    assert not (tmp_path / "d1.h5").exists(), "a 1D lattice still wrote an HDF5"


@pytest.mark.parametrize("phrase", [
    "not implemented",
    "leaves the second direction uncoupled",
    "Section_4_G_molecules",
    "resolution must stay coarser than the level spacing",
])
def test_the_refusal_explains_what_to_do_instead(tmp_path, phrase):
    # The message is wrapped, so compare on collapsed whitespace.
    result = _run_export(_ONE_D, tmp_path)
    # SystemExit's message goes to stderr; KITE's banner goes to stdout.
    output = " ".join((result.stdout + result.stderr).split())
    assert phrase in output, (
        f"the 1D refusal no longer says {phrase!r}; without it the user is told "
        "no but not told how")


@pytest.mark.skipif(not (KITEX.exists() and TOOLS.exists()),
                    reason="needs a built KITEx and KITE-tools")
def test_the_suggested_workaround_reproduces_the_exact_1d_dos(tmp_path):
    """The advice in the message is a promise. This is the promise."""
    assert _run_export(_CHAINS, tmp_path).returncode == 0, "the suggested lattice failed to export"
    for binary in (KITEX, TOOLS):
        subprocess.run([str(binary), "chains.h5"], cwd=str(tmp_path),
                       capture_output=True, check=True)

    data = np.loadtxt(tmp_path / "dos.dat")
    energy, rho = data[:, 0], data[:, 1]
    assert np.trapezoid(rho, energy) == pytest.approx(1.0, abs=0.02), (
        "the density of states is not normalised; if it is identically zero the "
        "chains are not being propagated at all")

    window = np.abs(energy) < 1.8          # stay off the van Hove edges
    exact = 1 / (np.pi * np.sqrt(4 - energy[window] ** 2))
    error = np.median(np.abs(rho[window] - exact) / exact)
    assert error < 0.10, (
        f"the decoupled-chain construction gives a median {error:.1%} error "
        "against the exact 1D tight-binding DOS, so the workaround the refusal "
        "recommends no longer produces 1D physics")
