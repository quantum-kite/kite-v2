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
