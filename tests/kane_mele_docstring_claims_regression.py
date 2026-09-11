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
    return ast.get_docstring(ast.parse(path.read_text()))


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


def test_raising_eta_to_the_recommended_value_reaches_sub_percent():
    """The docstring claims 0.08 gives 0.6% and stays well inside the gap."""
    residual = np.exp(-256 * 0.08 / 4.0)
    assert residual == pytest.approx(0.006, abs=0.001)
    assert 0.08 < (6 * np.sqrt(3) * 0.15) / 10, "eta must stay well inside the gap"


def test_finite_size_is_not_the_explanation():
    """Whatever the docstring blames, it must not be L."""
    doc = _docstring(PROCESS)
    assert "not an exactly quantized number at L=64" not in doc
    assert "NOT" in doc and "finite size" in doc
    assert "bare" in doc.lower() or "BARE" in doc


def test_the_stated_tolerance_is_achievable():
    """A docstring must not promise tighter than the reconstruction allows."""
    doc = _docstring(PROCESS)
    assert "about 10%" in doc, "the shipped-parameter tolerance must be stated"
    residual = np.exp(-256 * 0.04 / 4.0)
    assert residual < 0.10, "10% must be achievable at the shipped parameters"


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
