"""custom_two and conductivity_dc must report the same number for the same pair.

The suite already had `custom_two_same_vertex_sign_regression.py`, which asserts a
correlation above 0.7 between the two. Correlation is invariant under any positive
rescaling and under nothing else, so it cannot see a magnitude error and it cannot
see a sign error that is uniform across the curve. Both were present: custom_two
carried a literal factor 2 in its normalisation and transposed the moment matrix on
read, and the pair sat at -1.9954 times conductivity_dc while that test passed.

The model here is spinless on purpose. Haldane has no spin degree of freedom, so no
degeneracy, no s_z normalisation and no units convention can account for a
difference: the two reconstructions must agree. Both are requested from the same
Calculation, so one KITEx run produces both from one lattice at pinned seeds, and
the only thing that differs between the two numbers is the reconstruction.
"""

import importlib.util
import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITE_TOOLS.exists()), reason="needs a built KITEx/KITE-tools")

MOMENTS = 64
SEED = 20261101


def _haldane():
    """The shipped model, imported rather than retyped."""
    spec = importlib.util.spec_from_file_location(
        "hal", REPO / "examples" / "dos_dccond_haldane.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.haldane()


@pytest.fixture(scope="module")
def both_reconstructions(tmp_path_factory):
    kite = pytest.importorskip("kite")
    from kite import custom

    work = tmp_path_factory.mktemp("custom_two_vs_dc")
    h5 = work / "hal.h5"

    configuration = kite.Configuration(
        divisions=[2, 2], length=[128, 128],
        boundaries=["periodic", "periodic"], is_complex=True,
        precision=0, spectrum_range=[-10, 10],
        seed_h=SEED, seed_v=SEED)
    calculation = kite.Calculation(configuration)
    calculation.conductivity_dc(
        num_points=200, num_moments=MOMENTS, num_random=4,
        num_disorder=1, direction="xy", temperature=0.05)
    calculation.custom_two(
        stream_=[custom.Vertex(MOMENTS, [[1.0, "vx"]]),
                 custom.Vertex(MOMENTS, [[1.0, "vy"]])],
        num_random_=4, num_disorder_=1, num_points_=200, temperature_=0.05)
    kite.config_system(_haldane(), configuration, calculation, filename=str(h5))

    env = {"OMP_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"}
    subprocess.run([str(KITEX), str(h5)], cwd=work, check=True,
                   capture_output=True, env=env)
    subprocess.run([str(KITE_TOOLS), str(h5), "-t", "1"], cwd=work, check=True,
                   capture_output=True, env=env)

    dc = np.loadtxt(work / "condDC.dat")
    custom_two = np.loadtxt(work / "custom_two.dat")
    return dc, custom_two


def _ratio_on_common_grid(dc, custom_two):
    """The two readers use different Fermi grids; compare at the same energies."""
    energies = custom_two[:, 0]
    reference = np.interp(energies, dc[:, 0], dc[:, 1])
    carries_signal = np.abs(reference) > 0.25 * np.abs(reference).max()
    assert carries_signal.sum() > 20, "nothing above the noise floor to compare"
    return custom_two[carries_signal, 1] / reference[carries_signal]


def test_the_two_reconstructions_agree(both_reconstructions):
    """Equality, not correlation: the factor and the sign both have to be right."""
    ratio = np.median(_ratio_on_common_grid(*both_reconstructions))
    assert ratio == pytest.approx(1.0, abs=0.05), (
        f"custom_two/conductivity_dc = {ratio:+.4f}; a spinless model leaves no "
        "convention that could explain a difference")


def test_the_comparison_can_see_a_factor(both_reconstructions):
    """The assertion above would fail on the defect it was written for."""
    dc, custom_two = both_reconstructions
    corrupted = custom_two.copy()
    corrupted[:, 1] *= -2.0
    ratio = np.median(_ratio_on_common_grid(dc, corrupted))
    assert ratio != pytest.approx(1.0, abs=0.05), (
        "a -2x error has to be visible, or this file proves nothing")
