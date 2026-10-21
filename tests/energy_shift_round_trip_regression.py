"""A dropped energy shift is invisible unless spectrum_range is asymmetric.

KITE rescales the Hamiltonian to the Chebyshev interval as (H - shift)/scale.
Every entry point that takes an energy must apply the SAME affine map to it.
Omitting the shift evaluates at E + shift instead of E -- and is completely
invisible whenever spectrum_range is symmetric about zero, because shift is then
zero. That is how it reached production in `SimulationLDoS.cpp` once before (the
comment at `:90-99` records it) and how BUG-010 reached `custom_singleshot_two`.

The check here is convention-independent: it does not need to know the right
answer, only that **two valid rescalings of one model must agree**. Both
spectrum_ranges are chosen with the SAME WIDTH, so `energy_scale` is identical
and the rescaled broadenings and resolution are identical too; only `shift`
differs. Any disagreement is therefore the shift alone.

It compares the *position* of a spectral feature rather than its amplitude.
`seed_h` does not seed the random-vector stream, so two runs of the same input
file differ by a few percent and an amplitude assertion would either flake or be
too loose to catch anything. An energy offset moves the feature; noise does not.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

#: Square lattice with a NON-ZERO on-site energy, so the spectrum is not
#: symmetric about zero and the shift is not trivially zero for both ranges.
_BUILD = '''
import sys, numpy as np, kite
from kite import lattice as latt
from kite import custom
lo, hi, out = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
lattice = latt.Lattice(a1=[1, 0], a2=[0, 1])
lattice.add_sublattices(['A', [0, 0], 1.0])
lattice.add_hoppings([[-1, 0], 'A', 'A', -1.0], [[0, -1], 'A', 'A', -1.0])
configuration = kite.Configuration(
    divisions=[1, 1], length=[64, 64], boundaries=["periodic"] * 2,
    is_complex=True, precision=1, spectrum_range=[lo, hi])
calculation = kite.Calculation(configuration)
A = custom.Vertex(64, [[1.0, "vx"]])
B = custom.Vertex(64, [[1.0, "vx"]])
calculation.custom_singleshot_two(
    stream_=[A, B], num_random_=8, num_disorder_=1, gamma_=0.1, sigma_=0.1,
    energies_=np.linspace(0.2, 1.8, 9))
kite.config_system(lattice, configuration, calculation, filename=out)
'''

ENERGIES = np.linspace(0.2, 1.8, 9)


def _sweep(tmp_path, lo, hi, name):
    import h5py

    subprocess.run([sys.executable, "-c", _BUILD, str(lo), str(hi), name],
                   cwd=str(tmp_path), capture_output=True, check=True)
    subprocess.run([str(KITEX), name], cwd=str(tmp_path),
                   capture_output=True, check=True)
    with h5py.File(tmp_path / name, "r") as handle:
        return np.real(np.asarray(
            handle["/Calculation/CustomSingleTwo/Result"])).ravel()


def test_two_spectrum_ranges_of_equal_width_put_the_feature_at_the_same_energy(tmp_path):
    """BUG-010. Same width, so only the shift differs; the peak must not move."""
    symmetric = _sweep(tmp_path, -6, 6, "sym.h5")     # scale 6, shift 0
    asymmetric = _sweep(tmp_path, -5, 7, "asym.h5")   # scale 6, shift 1

    peak_sym = ENERGIES[np.argmax(np.abs(symmetric))]
    peak_asym = ENERGIES[np.argmax(np.abs(asymmetric))]
    step = ENERGIES[1] - ENERGIES[0]

    assert abs(peak_sym - peak_asym) <= step, (
        f"the feature sits at E={peak_sym:.2f} with spectrum_range=[-6,6] but at "
        f"E={peak_asym:.2f} with [-5,7]. Both ranges have the same width, so "
        f"energy_scale and every rescaled broadening are identical and only the "
        f"shift differs -- a displacement of {peak_asym - peak_sym:+.2f} against a "
        f"shift of 1.0 means the shift is being applied to one and not the other")


def test_every_driver_that_reads_the_scale_also_reads_the_shift(tmp_path):
    """The static half, stated as a rule rather than a pattern match.

    A first attempt here searched for an energy-named variable divided by
    `energy_scale`, and missed BUG-010 entirely -- the broken line read
    `const value_type en = vv / energy_scale;` and `vv` is not energy-named. So
    the rule is structural instead: a driver that reads `/EnergyScale` to rescale
    a user-supplied *energy* must also read `/EnergyShift`.

    The exemptions are explicit and each is a quantity that is a difference, for
    which no shift applies: `Hamiltonian.cpp` rescales an additive on-site
    correction, `GlobalSimulation.cpp` only reads the value, and
    `SimulationSingleShot.cpp` carries it through for KITE-tools to apply on
    output. `SimulationLCM` and `SimulationSLCM` are genuine omissions of the
    same kind as BUG-010, kept here because they are experimental, have no Python
    entry point, and would otherwise be forgotten.
    """
    import re

    exempt_difference_only = {
        "Hamiltonian.cpp", "GlobalSimulation.cpp", "SimulationSingleShot.cpp",
    }
    known_omissions = {"SimulationLCM.cpp", "SimulationSLCM.cpp"}

    missing = []
    for source in list((REPO / "Src").rglob("*.cpp")):
        text = "\n".join(re.sub(r"//.*$", "", line)
                          for line in source.read_text(errors="ignore").splitlines())
        if '"/EnergyScale"' not in text:
            continue
        if source.name in exempt_difference_only:
            continue
        if '"/EnergyShift"' not in text:
            missing.append(source.name)

    unexpected = sorted(set(missing) - known_omissions)
    assert not unexpected, (
        f"{unexpected} read /EnergyScale to rescale an energy but never read "
        "/EnergyShift. That evaluates at E + shift and is invisible whenever "
        "spectrum_range is symmetric about zero -- BUG-010's exact failure mode")

    # And the converse: if one of the known omissions is ever fixed, say so, so
    # that this list does not quietly outlive the problem.
    fixed = sorted(n for n in known_omissions
                   if '"/EnergyShift"' in (REPO / "Src" / "Simulation" / n).read_text())
    assert not fixed, (
        f"{fixed} now read /EnergyShift; remove them from known_omissions")
