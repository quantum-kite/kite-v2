"""Static sweep: every example's declared spectrum_range must contain its spectrum.

No KITEx run needed. Each example hands over its lattice and its declared range
(by intercepting `config_system`), and the range is checked against the band
extrema computed from the lattice itself.

**What this test is careful about, because getting it wrong produces noise that
buries the one real defect.** Two distinct bounds are in play:

  tight  = clean band extrema + on-site disorder spread. Rigorous, and close to
           the truth: Anderson disorder really does shift every site.
  loose  = tight + the structural-defect term, max_node(|dU| + sum|dt|). Also
           rigorous, by Weyl, but a worst case that a single impurity does not
           come near -- a rank-1 perturbation moves at most one eigenvalue, and
           it moves toward the sign of the perturbation, not both ways.

Only `tight` is a failure. Falling between the two is reported, not failed:
checked by exact diagonalization, `ldos_spin_operator_validation.py`'s
`limp = -3.0` impurity puts the true band edge at -4.2834 against a clean
-4.0792 and a Weyl worst case of -7.0792. Its declared +-6 is correct, and a
test using `loose` as the criterion would call it broken.

Margin is likewise reported, not enforced. `dos_graphene.py` declares
`3.01*t` -- 0.3% headroom, deliberate and documented in the example itself.
"""

import contextlib
import importlib.util
import io
from pathlib import Path

import numpy as np
import pytest

import kite

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = sorted((REPO / "examples").glob("*.py"))


class _Captured(Exception):
    """Raised to stop an example once it has handed over its configuration."""


def _capture(path):
    """Run an example far enough to see its config_system call, writing nothing."""
    captured = {}
    real = kite.config_system

    def spy(lattice, configuration, calculation, **kwargs):
        captured.update(lattice=lattice, configuration=configuration, kwargs=kwargs)
        raise _Captured()

    kite.config_system = spy
    try:
        spec = importlib.util.spec_from_file_location("swept_" + path.stem, path)
        module = importlib.util.module_from_spec(spec)
        with contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            try:
                spec.loader.exec_module(module)
            except _Captured:
                pass
            except BaseException:
                return None
            if not captured and hasattr(module, "main"):
                try:
                    module.main()
                except _Captured:
                    pass
                except BaseException:
                    return None
    finally:
        kite.config_system = real
    return captured or None


def _cases():
    for path in EXAMPLES:
        captured = _capture(path)
        if not captured:
            continue
        configuration = captured["configuration"]
        if not configuration.energy_scale:
            continue                       # automatic scaling is computed, not declared
        yield pytest.param(captured, id=path.name)


CASES = list(_cases())


@pytest.mark.parametrize("captured", CASES)
def test_declared_range_contains_the_spectrum(captured):
    lattice = captured["lattice"]
    configuration = captured["configuration"]
    kwargs = captured["kwargs"]

    declared_lo = configuration.energy_shift - configuration.energy_scale
    declared_hi = configuration.energy_shift + configuration.energy_scale

    tight_lo, tight_hi = kite.estimate_bounds(
        lattice, disorder=kwargs.get("disorder"),
        disorder_structural=None,                    # the loose term, deliberately excluded
        magnetic_field=bool(kwargs.get("modification")), margin=0.0)

    assert declared_lo <= tight_lo and declared_hi >= tight_hi, (
        f"declared [{declared_lo:.3f}, {declared_hi:.3f}] does not contain the "
        f"spectrum [{tight_lo:.3f}, {tight_hi:.3f}]. The rescaled spectrum would "
        f"leave [-1, 1], where T_n grows exponentially and the run still exits 0.")


def test_the_sweep_actually_covers_the_examples():
    """A sweep that silently captures nothing would pass forever."""
    assert len(CASES) >= 20, f"only {len(CASES)} examples captured; the sweep is not working"


def test_the_sweep_catches_a_range_that_is_too_narrow():
    """The F-63 configuration: a third of the true bandwidth."""
    from kite import lattice as latt
    a_cc = 1 / np.sqrt(3)
    lat = latt.Lattice(a1=[1, 0], a2=[0.5, np.sqrt(3) / 2])
    lat.add_sublattices(("A", [0, -a_cc / 2], 0.0), ("B", [0, a_cc / 2], 0.0))
    lat.add_hoppings(([0, 0], "A", "B", -2.8), ([-1, 0], "A", "B", -2.8),
                     ([-1, 1], "A", "B", -2.8))
    tight_lo, tight_hi = kite.estimate_bounds(lat, margin=0.0)
    declared_lo, declared_hi = -2.8, 2.8              # what ldos_graphene.py used to declare
    assert not (declared_lo <= tight_lo and declared_hi >= tight_hi), (
        "the sweep must reject a range covering a third of the bandwidth")
