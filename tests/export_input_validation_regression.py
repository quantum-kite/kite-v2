"""Two export-time guards: the spinor length, and an unpopulated Disorder.

Both fail in Python, where the message can name the mistake, rather than
downstream where it cannot.

**Spinor.** The seeded wave packet carries one amplitude per orbital of the unit
cell. KITEx reads a fixed stride, so a spinor of the wrong length is not caught
there -- it silently seeds the wrong orbitals, which is a wrong answer with no
signal.

**Disorder.** A `Disorder` object constructed but never given an
`add_disorder()` call leaves `_orbital` as an empty list, and
`disorder._orbital > -1` then raised

    TypeError: '>' not supported between instances of 'list' and 'int'

naming neither the object nor the omission. This was already fail-closed -- it
raised rather than exporting something wrong -- so the fix is the message, not
the behaviour, and the test asserts the message rather than the exception type.
"""

import numpy as np
import pytest

import kite
from kite import lattice as latt


def _square(num_sublattices=1):
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(*[(chr(65 + i), [0, 0], 0.0) for i in range(num_sublattices)])
    for i in range(num_sublattices - 1):
        lat.add_hoppings(([0, 0], chr(65 + i), chr(65 + i + 1), 0.1))
    lat.add_hoppings(([1, 0], "A", "A", -1.0))
    return lat


def _config(is_complex=False):
    return kite.Configuration(
        divisions=[1, 1], length=[16, 16], boundaries=["periodic", "periodic"],
        is_complex=is_complex, precision=1, spectrum_range=[-5, 5])


def _wave_packet(calculation, spinor_length):
    calculation.gaussian_wave_packet(
        num_points=10, num_moments=16, timestep=0.5, num_disorder=1,
        k_vector=np.array([[0.0, 0.0]]),
        spinor=np.array([[1.0 + 0j] * spinor_length]),
        width=2.0, mean_value=[8, 8])


@pytest.mark.parametrize("num_orbitals", [1, 2, 4])
def test_matching_spinor_is_accepted(num_orbitals, tmp_path):
    cfg = _config(is_complex=True)
    calc = kite.Calculation(cfg)
    _wave_packet(calc, num_orbitals)
    kite.config_system(_square(num_orbitals), cfg, calc,
                       filename=str(tmp_path / "ok.h5"))


@pytest.mark.parametrize("num_orbitals,spinor_length", [(2, 1), (2, 3), (4, 2), (1, 2)])
def test_mismatched_spinor_is_refused(num_orbitals, spinor_length, tmp_path):
    cfg = _config(is_complex=True)
    calc = kite.Calculation(cfg)
    _wave_packet(calc, spinor_length)
    with pytest.raises(SystemExit) as e:
        kite.config_system(_square(num_orbitals), cfg, calc,
                           filename=str(tmp_path / "bad.h5"))
    message = str(e.value)
    assert str(num_orbitals) in message and str(spinor_length) in message, (
        "the message must name both the spinor length and the orbital count")


def test_unpopulated_disorder_names_the_omission(tmp_path):
    lat = _square()
    cfg = _config()
    calc = kite.Calculation(cfg)
    calc.dos(num_points=100, num_moments=32, num_random=1)
    with pytest.raises(SystemExit) as e:
        kite.config_system(lat, cfg, calc, disorder=kite.Disorder(lat),
                           filename=str(tmp_path / "empty.h5"))
    message = str(e.value)
    assert "add_disorder" in message, "the message must say what to call"
    assert "not supported between instances" not in message, "no raw TypeError"


def test_populated_disorder_still_exports(tmp_path):
    lat = _square()
    cfg = _config()
    calc = kite.Calculation(cfg)
    calc.dos(num_points=100, num_moments=32, num_random=1)
    disorder = kite.Disorder(lat)
    disorder.add_disorder("A", "Uniform", 0.0, 0.1)
    kite.config_system(lat, cfg, calc, disorder=disorder,
                       filename=str(tmp_path / "ok.h5"))


def test_no_disorder_argument_still_exports(tmp_path):
    """The guard must not fire on a clean system."""
    lat = _square()
    cfg = _config()
    calc = kite.Calculation(cfg)
    calc.dos(num_points=100, num_moments=32, num_random=1)
    kite.config_system(lat, cfg, calc, filename=str(tmp_path / "clean.h5"))
