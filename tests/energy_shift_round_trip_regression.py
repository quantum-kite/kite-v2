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


def test_kitex_does_not_rescale_user_supplied_energies(tmp_path):
    """Python owns the rescaling. KITEx reads the shift in exactly one place.

    After T-073 every user-supplied energy and width is rescaled in
    `config_system` before export, so no driver needs `/EnergyShift` to
    interpret its input. One file still reads it, for a reason that is physics
    rather than unit conversion:

      Both wave-packet drivers form `exp(-i * energy_shift * t)`, the global phase
      the shift contributes to time evolution:
      `SimulationLocalizedWavePacket.cpp:464` and the Bessel-coefficient block of
      `SimulationGaussianWavePacket.cpp`. That cannot move to Python -- it is a
      phase applied per timestep inside the propagator, not a rescaling of input.

      The same file also rescales its `energy_window` at `:66,71`. That one was
      deliberately left in C++: the "no filter" sentinel is `(0, 0)` in physical
      units, so pre-rescaling would turn it into `(-shift/scale, -shift/scale)`
      and silently switch filtering on for any asymmetric `spectrum_range`.
      Moving it means changing the sentinel too.

    Anything else reading `/EnergyShift` is a driver re-doing work Python has
    already done -- the double-rescale half of T-073, as opposed to BUG-010's
    dropped-shift half.
    """
    allowed = {"SimulationLocalizedWavePacket.cpp",
               "SimulationGaussianWavePacket.cpp"}
    readers = {f.name for f in (REPO / "Src").rglob("*.cpp")
               if '"/EnergyShift"' in f.read_text(errors="ignore")}

    unexpected = sorted(readers - allowed)
    assert not unexpected, (
        f"{unexpected} read /EnergyShift. Since T-073 the shift is applied in "
        "python/kite/__init__.py before export, so a driver that reads it is "
        "either rescaling twice or disagreeing with the exporter")

    for driver in allowed:
        assert driver in readers, (
            f"{driver} no longer reads /EnergyShift. Its propagator then computes "
            "exp(-i*Htilde*ttilde) instead of exp(-i*H*t): unobservable in any "
            "expectation value, but the state itself is wrong")


def test_the_exporter_rescales_every_energy_it_writes():
    """The other side of the same invariant, checked at the source.

    Every dataset `config_system` writes that holds a user-supplied *absolute*
    energy must carry `- config.energy_shift`. A width carries only the scale.
    """
    import re

    text = (REPO / "python" / "kite" / "__init__.py").read_text()
    pattern = re.compile(r"create_dataset\('(Energy|Energies)'[^\n]*", re.M)
    absolute = list(pattern.finditer(text))
    assert absolute, "config_system no longer writes an Energy dataset"
    for m in absolute:
        assert "energy_shift" in m.group(0), (
            "an Energy/Energies dataset is written without subtracting "
            f"config.energy_shift: {m.group(0)[:110]!r}")
