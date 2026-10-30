"""The moment count must be rejected in Python, not overrun in C++.

Two distinct constraints, imposed by how the C++ accumulation loops advance:

  pair-stepping   n += 2        dos, ldos, arpes, custom_one, custom_one_local
  block-stepping  n += MEMORY   conductivity_dc, conductivity_optical,
                                conductivity_optical_nonlinear, custom_two,
                                custom_two_local

They are not interchangeable, and that is the whole point of F-01: 1000 is
even, so an evenness check passes it, but 1000 mod 16 = 8 and the loop overruns
by 8008 elements. Whether that faults depends on where the allocation lands --
at num_moments=1000 it crashes, at 520 it completes and returns numbers, which
is the more dangerous outcome.

Deliberately unconstrained here, and asserted so:

  * singleshot_conductivity_dc steps by one;
  * ldos_map, spectral_map, custom_singleshot_two, lcm and st_lcm carry no
    num_moments at all -- resolution comes from a width parameter.

A check that over-rejects would be its own defect.
"""

import numpy as np
import pytest

import kite


def _config():
    return kite.Configuration(
        divisions=[1, 1], length=[16, 16], boundaries=["periodic", "periodic"],
        is_complex=False, precision=1, spectrum_range=[-4, 4])


@pytest.mark.parametrize("num_moments", [2, 64, 512, 6144])
def test_even_counts_accepted_for_pair_stepping(num_moments):
    kite.Calculation(_config()).dos(
        num_points=100, num_moments=num_moments, num_random=1)


@pytest.mark.parametrize("num_moments", [1, 63, 511, 6143])
def test_odd_counts_refused_for_pair_stepping(num_moments):
    with pytest.raises(SystemExit) as e:
        kite.Calculation(_config()).dos(
            num_points=100, num_moments=num_moments, num_random=1)
    assert "EVEN" in str(e.value)
    assert str(num_moments) in str(e.value), "the message must name the value given"


@pytest.mark.parametrize("num_moments", [16, 512, 1024, 6144])
def test_multiples_of_memory_accepted_for_block_stepping(num_moments):
    kite.Calculation(_config()).conductivity_dc(
        direction="xy", num_points=100, num_moments=num_moments, num_random=1)


@pytest.mark.parametrize("num_moments", [1000, 520, 8, 100])
def test_non_multiples_refused_for_block_stepping(num_moments):
    """1000 and 520 are the two F-01 cases: both even, neither a multiple of 16.

    1000 crashed under AddressSanitizer; 520 completed and returned numbers.
    """
    with pytest.raises(SystemExit) as e:
        kite.Calculation(_config()).conductivity_dc(
            direction="xy", num_points=100, num_moments=num_moments, num_random=1)
    assert "multiple of 16" in str(e.value)


def test_evenness_alone_does_not_satisfy_the_block_constraint():
    """The exact hole F-01 fell through."""
    assert 1000 % 2 == 0 and 1000 % 16 != 0
    with pytest.raises(SystemExit):
        kite.Calculation(_config()).conductivity_dc(
            direction="xy", num_points=100, num_moments=1000, num_random=1)


@pytest.mark.parametrize("bad", [0, -16, 3.5, "512", None])
def test_nonsense_counts_are_refused(bad):
    with pytest.raises(SystemExit):
        kite.Calculation(_config()).dos(num_points=100, num_moments=bad, num_random=1)


def test_the_message_suggests_a_usable_value():
    """An error a user cannot act on is half an error."""
    with pytest.raises(SystemExit) as e:
        kite.Calculation(_config()).conductivity_dc(
            direction="xy", num_points=100, num_moments=1000, num_random=1)
    message = str(e.value)
    assert "992" in message or "1008" in message, "suggest the nearest valid counts"

# --- every entry point, not just the two that were wired up first -----------
#
# The check and its table existed from the start; only `dos`, `arpes` and
# `conductivity_dc` ever called it, so the other seven exported a bad count
# silently and overran in C++. These parametrisations are the guard against
# that happening again one entry point at a time.

def _spin_calculation():
    """A Calculation with l0 registered, for the custom_* entry points."""
    calc = kite.Calculation(_config())
    calc.add_orbital_index("A", 0)
    calc.add_orbital_coupling("A", "A", 1.0, "l0")
    return calc


def _call(name, moments):
    """Invoke one entry point with this moment count."""
    from kite import custom
    calc = kite.Calculation(_config())
    if name == "dos":
        calc.dos(num_points=100, num_moments=moments, num_random=1)
    elif name == "ldos":
        calc.ldos(energy=[0.0], num_moments=moments, position=[[8, 8]],
                  sublattice=["A"])
    elif name == "arpes":
        calc.arpes(k_vector=np.array([[0.0, 0.0]]), weight=[1.0],
                   num_moments=moments, num_disorder=1)
    elif name == "custom_one":
        _spin_calculation().custom_one(
            stream_=custom.Vertex(moments, [[1.0, "l0"]]),
            num_random_=1, num_disorder_=1)
    elif name == "custom_one_local":
        _spin_calculation().custom_one_local(
            stream_=custom.Vertex(moments, [[1.0, "l0"]]),
            energy_=np.array([0.0]), position_=[[8, 8]], sublattice_=["A"])
    elif name == "conductivity_dc":
        calc.conductivity_dc(direction="xy", num_points=100,
                             num_moments=moments, num_random=1)
    elif name == "conductivity_optical":
        calc.conductivity_optical(direction="xx", num_points=100,
                                  num_moments=moments, num_random=1)
    elif name == "conductivity_optical_nonlinear":
        calc.conductivity_optical_nonlinear(direction="xxx", num_points=100,
                                            num_moments=moments, num_random=1)
    elif name == "custom_two":
        _spin_calculation().custom_two(
            stream_=[custom.Vertex(moments, [[1.0, "l0"]]),
                     custom.Vertex(moments, [[1.0, "l0"]])],
            num_random_=1, num_disorder_=1, num_points_=100, temperature_=0.01)
    elif name == "custom_two_local":
        _spin_calculation().custom_two_local(
            stream_=[custom.Vertex(moments, [[1.0, "l0"]]),
                     custom.Vertex(moments, [[1.0, "l0"]])],
            positions_=[[8, 8]])
    else:
        raise AssertionError(name)


PAIR = ["dos", "ldos", "arpes", "custom_one", "custom_one_local"]
BLOCK = ["conductivity_dc", "conductivity_optical",
         "conductivity_optical_nonlinear", "custom_two", "custom_two_local"]


@pytest.mark.parametrize("name", PAIR)
def test_every_pair_stepping_entry_point_refuses_an_odd_count(name):
    with pytest.raises(SystemExit) as e:
        _call(name, 63)
    assert "63" in str(e.value), f"{name}'s message must name the value given"


@pytest.mark.parametrize("name", BLOCK)
def test_every_block_stepping_entry_point_refuses_a_non_multiple(name):
    """200 is even, which is exactly what made this survive."""
    with pytest.raises(SystemExit) as e:
        _call(name, 200)
    assert "16" in str(e.value), (
        f"{name}'s message must say what to change 200 to")


@pytest.mark.parametrize("name", PAIR + BLOCK)
def test_no_entry_point_refuses_a_valid_count(name):
    """A check that over-rejects would be its own defect."""
    _call(name, 192)
