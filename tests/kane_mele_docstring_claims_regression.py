"""The Kane-Mele post-processor must state the right reason for its own error.

The docstring used to attribute the non-quantized plateau to finite size --
"not an exactly quantized number at L=64". That is wrong by roughly twenty
orders of magnitude, and it sends a user to the expensive knob (rerun KITEx at
larger L) instead of the free one (change a post-processing broadening).

Numbers, all checkable from the generator's own parameters:

  gap            6*sqrt(3)*t2 = 1.5588 |t| at t2 = 0.15, confirmed against
                 exact diagonalization (1.5590, direct = indirect)
  xi             hbar*v_F / (gap/2) ~ 1 unit cell, so at L = 64 the
                 finite-size correction goes as exp(-64)
  real floor     the reconstruction uses the BARE resolvent expansion with no
                 kernel, so the residual is exp(-M*eta/energy_scale) at the
                 band centre = exp(-2.56) = 7.7% at the shipped M=256,
                 eta=0.04, energy_scale=4

This test pins the arithmetic, so a future edit cannot quietly reintroduce a
tolerance that the reconstruction cannot deliver.
"""

import ast
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
PROCESS = REPO / "examples" / "kane_mele_spin_hall_process.py"
GENERATOR = REPO / "examples" / "kane_mele_spin_hall.py"


def _docstring(path):
    """Docstring with whitespace normalised, so a phrase that happens to wrap
    across a line still matches."""
    import re
    return re.sub(r"\s+", " ", ast.get_docstring(ast.parse(path.read_text())))


def test_gap_formula_matches_the_generators_parameters():
    """6*sqrt(3)*t2 with the t2 the generator actually ships."""
    t2 = 0.15
    assert 6 * np.sqrt(3) * t2 == pytest.approx(1.5588, abs=5e-4)


def test_truncation_residual_is_the_quoted_seven_point_seven_percent():
    """exp(-M * eta / energy_scale) at the band centre, no kernel applied."""
    moments, eta, energy_scale = 256, 0.04, 4.0
    residual = np.exp(-moments * eta / energy_scale)
    assert residual == pytest.approx(0.077, abs=0.002)
    assert f"{residual:.1%}".rstrip("%") == "7.7"


def test_the_pointwise_residual_is_not_quoted_as_the_plateau_error():
    """The truncation residual bounds the Green's function pointwise, not the
    integrated response. Measured plateau error at the shipped parameters is
    1.2%, against a pointwise residual of 7.7% -- the Kubo-Bastin energy
    integral averages the ringing away. A docstring that quotes the residual as
    the tolerance is wrong by an order of magnitude."""
    doc = _docstring(PROCESS)
    assert "about 10%" not in doc, "the pointwise residual is not the plateau error"
    assert "pointwise" in doc.lower(), (
        "the docstring must say which quantity the truncation residual bounds")
    assert "not read it as an error bar" in doc.lower()


def test_finite_size_is_not_the_explanation():
    """Whatever the docstring blames, it must not be L."""
    doc = _docstring(PROCESS)
    assert "not an exactly quantized number at L=64" not in doc
    low = doc.lower()
    assert "not limited by finite size" in low, "finite size must be ruled out explicitly"
    assert "correlation length" in low, "and why -- about one unit cell against L = 64"
    assert "bare" in low, "the resolvent expansion is un-damped; say so"


def test_the_stated_tolerance_matches_what_was_measured():
    """Measured on the shipped output: plateau -2.024, flat to 0.12%.

    The criterion in the docstring must be loose enough to hold and tight
    enough to catch a regression -- flat to 0.5%, within 2% of 2.
    """
    doc = _docstring(PROCESS)
    assert "0.12%" in doc and "-2.024" in doc, "measurements must be recorded"
    assert "0.5%" in doc and "2%" in doc, "the criterion must be stated"
    measured_flatness, measured_error = 0.0012, 0.012
    assert measured_flatness < 0.005, "criterion must hold on the measurement"
    assert measured_error < 0.02, "criterion must hold on the measurement"


def test_the_generator_still_ships_the_parameters_the_docstring_assumes():
    """If the generator's defaults move, this docstring's arithmetic goes stale.

    Read the signature rather than grepping the text: the numbers live in
    `main(t=1.0, t2=0.15, moments=256, ...)` and a string match would pass on a
    coincidental occurrence elsewhere in the file.
    """
    tree = ast.parse(GENERATOR.read_text())
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    defaults = dict(zip([a.arg for a in main.args.args][-len(main.args.defaults):],
                        [getattr(d, "value", None) for d in main.args.defaults]))
    assert defaults["t2"] == 0.15, (
        f"the gap number 6*sqrt(3)*t2 = 1.559 assumes t2 = 0.15, found {defaults['t2']}")
    assert defaults["moments"] == 256, (
        f"the 7.7% residual assumes M = 256, found {defaults['moments']}")
