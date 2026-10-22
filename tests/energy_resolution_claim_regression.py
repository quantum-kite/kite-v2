"""The documented Jackson resolution must match the kernel's actual width.

`docs/documentation/optimization.md` tells the reader how finely an order-$M$
expansion resolves energy. The claim is load-bearing: it is how a user chooses
`num_moments`, and the page's own Rule of Thumb warns that too large an $M$
resolves individual levels and produces finite-size noise.

It used to state a single width, `pi*de/M`, which reads as uniform across the
band. It is not uniform -- the Jackson-broadened delta narrows towards the band
edges as `sqrt(1 - Etilde^2)`, so `pi*de/M` is the band-centre *maximum*.

This test rebuilds the Jackson-regularized delta from its definition and fits
its standard deviation, rather than trusting either the page or the code. If the
two ever disagree the page is wrong, and the page is what users act on.
"""

from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
PAGE = REPO / "docs" / "documentation" / "optimization.md"


def _jackson(n, moments):
    """The Jackson kernel coefficient g_n, from its closed form."""
    from math import cos, pi, sin, tan

    return ((moments - n + 1) * cos(pi * n / (moments + 1))
            + sin(pi * n / (moments + 1)) / tan(pi / (moments + 1))) / (moments + 1)


def _measured_width(moments, centre):
    """Standard deviation of the Jackson-broadened delta at a rescaled energy."""
    x = np.linspace(-0.999, 0.999, 200001)
    delta = np.zeros_like(x)
    for n in range(moments):
        delta += ((2 - (n == 0)) * _jackson(n, moments)
                  * np.cos(n * np.arccos(centre)) * np.cos(n * np.arccos(x)))
    delta /= np.pi * np.sqrt(1 - x ** 2)
    delta = np.clip(delta, 0, None)
    delta /= np.trapezoid(delta, x)
    mean = np.trapezoid(x * delta, x)
    return np.sqrt(np.trapezoid((x - mean) ** 2 * delta, x))


@pytest.mark.parametrize("moments", [128, 512])
@pytest.mark.parametrize("centre", [0.0, 0.5, 0.9])
def test_the_resolution_carries_the_square_root_factor(moments, centre):
    """sigma(Etilde) = (pi/M) sqrt(1 - Etilde^2), in rescaled units."""
    predicted = np.pi * np.sqrt(1 - centre ** 2) / moments
    measured = _measured_width(moments, centre)
    assert measured == pytest.approx(predicted, rel=0.03), (
        f"at M={moments}, Etilde={centre} the Jackson delta has width "
        f"{measured:.6f} against the documented {predicted:.6f}")


def test_the_band_centre_is_the_widest_point():
    """The direction of the dependence, which is what the page now asserts."""
    centre = _measured_width(256, 0.0)
    edge = _measured_width(256, 0.9)
    assert edge < centre, (
        f"the expansion is not coarsest at the band centre ({centre:.6f} vs "
        f"{edge:.6f} at Etilde=0.9); optimization.md says it is")
    assert centre / edge == pytest.approx(1 / np.sqrt(1 - 0.9 ** 2), rel=0.05), (
        "the ratio between band-centre and band-edge width is not 1/sqrt(1-E^2)")


@pytest.mark.parametrize("phrase", [
    "not uniform across the band",
    r"\sqrt{1-\tilde\lambda^{2}}",
    "coarsest at the band centre",
    "standard deviation, not a full width at half maximum",
])
def test_the_page_states_the_energy_dependence(phrase):
    """AC-DOCS-02: the energy dependence and the edge behaviour must be stated."""
    assert phrase in PAGE.read_text(), (
        f"optimization.md no longer says {phrase!r}; a reader choosing "
        "num_moments from it would take the band-centre width as uniform")


def test_the_page_scopes_the_energy_dependence_to_the_kernel_route():
    """The Jackson width varies across the band; the CPGF broadening does not.

    KITE offers both. `-K jackson` damps a truncated expansion, and because a
    kernel of uniform width in `theta = arccos(x)` maps to `sqrt(1-x^2)` in
    energy, its resolution is energy-dependent. `-K green` is not a kernel: it
    is an analytic resolvent at finite `eta`, whose broadening is `eta`
    everywhere. Measured at M=256: Jackson gives sigma = 0.0122, 0.0106, 0.0054
    at Etilde = 0, 0.5, 0.9 -- a factor 2.3 across the band -- while the CPGF
    width is `eta` by construction.

    Stating the first without the second would leave a reader thinking every
    KITE broadening narrows at the edges.
    """
    text = PAGE.read_text()
    assert "belongs to the kernel route, not to KITE's broadening in" in text, (
        "optimization.md no longer scopes the energy dependence to the kernel "
        "route; the Green's-function route broadens uniformly")
    assert "uniform across the band" in text, (
        "the page no longer says the Green's-function broadening is uniform")
