"""The legacy fixtures are seeded, so a rerun reproduces its reference.

`SEED=ones`, which every legacy `test.py` sets, replaces the stochastic probe
vector with a vector of ones -- see the `seed=="ones"` branch in
`Src/Vector/KPM_Vector2D.cpp` and its 3D twin. It does not touch the random
numbers that *place* disorder: those come from `simul.rnd`, whose seeds are
`/Seed0` and `/Seed1` read out of the HDF5 file.

Those two fields are what this file guards, because zero is a trap. KITEx maps
them per thread as `(s0 != 0) * 2654435761 * (...)` in
`Src/Simulation/GlobalSimulation.cpp`, and `KPMRandom::init_random` seeds the
Mersenne twister deterministically only `if (seed_ > 0)`, drawing from
`std::random_device` otherwise. So `Seed0 = 0` does not mean "seed zero", it
means "do not seed", and a fixture written that way redraws its disorder on
every run.

Measured on `test016_KITEx_wy_cond`, which carries vacancies: with the seeds as
they ship, two runs of `Gammaxx` agree to 4.3e-19 against a norm of 1.1e-2 --
floating-point reduction noise and nothing else. The same configuration with
the seeds set to zero differs by 1.7e-4, which is four orders of magnitude
above the comparison tolerance of 1e-8 those tests use.
"""

import os
import re
import shutil
import subprocess

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

#: Every checked-in input a legacy test starts from.
FIXTURES = sorted((REPO / "tests").glob("*/configORIG.h5"))

#: The fixtures whose model places something at random, where an unseeded run
#: changes the Hamiltonian itself rather than only the probe vector.
RANDOMLY_PLACED = ("test011_KITEx_sq_defec", "test014_KITEx_gr_adat",
                   "test015_KITEx_wy_ss", "test016_KITEx_wy_cond")


def test_there_are_fixtures_to_check():
    assert len(FIXTURES) > 15, f"only {len(FIXTURES)} fixtures found; the glob is wrong"


@pytest.mark.parametrize("fixture", FIXTURES, ids=lambda p: p.parent.name)
def test_the_fixture_carries_both_seeds(fixture):
    with h5py.File(fixture, "r") as handle:
        for field in ("Seed0", "Seed1"):
            assert field in handle, (
                f"{fixture.parent.name} has no /{field}; KITEx reads it unconditionally "
                "and throws without it")


@pytest.mark.parametrize("fixture", FIXTURES, ids=lambda p: p.parent.name)
def test_neither_seed_is_zero(fixture):
    with h5py.File(fixture, "r") as handle:
        for field in ("Seed0", "Seed1"):
            value = int(handle[field][()])
            assert value > 0, (
                f"{fixture.parent.name} sets /{field} to {value}; zero is the "
                "non-reproducible setting, not a seed")


def test_zero_still_means_unseeded_in_the_source():
    """The guard above is only worth anything while this stays true."""
    mapping = (REPO / "Src" / "Simulation" / "GlobalSimulation.cpp").read_text()
    assert "(s0 != 0)" in mapping and "(s1 != 0)" in mapping, (
        "GlobalSimulation.cpp no longer zeroes the per-thread seed for seed 0; "
        "recheck what zero means before trusting this file")
    random = (REPO / "Src" / "Tools" / "Random.cpp").read_text()
    assert "if (seed_ > 0)" in random, (
        "KPMRandom::init_random no longer treats a positive seed as the "
        "deterministic case; recheck what zero means before trusting this file")


def test_ones_only_fixes_the_probe_vector():
    """Named here because it is the reason the seeds still matter."""
    for name in ("KPM_Vector2D.cpp", "KPM_Vector3D.cpp"):
        text = (REPO / "Src" / "Vector" / name).read_text()
        assert 'getenv("SEED")' in text and re.search(r'seed\s*==\s*"ones"', text), (
            f"{name} no longer recognises SEED=ones")
        assert "simul.rnd.init()" in text, (
            f"{name} no longer falls back to the seeded generator")


@pytest.mark.skipif(not KITEX.exists(), reason="needs a build to run")
@pytest.mark.parametrize("name", RANDOMLY_PLACED)
def test_a_seeded_fixture_reproduces_itself(name, tmp_path):
    source = REPO / "tests" / name / "configORIG.h5"
    outputs = []
    for run in ("first", "second"):
        work = tmp_path / run
        work.mkdir()
        target = work / "config.h5"
        shutil.copy(source, target)
        environment = dict(os.environ, SEED="ones")
        result = subprocess.run([str(KITEX), "config.h5"], cwd=work,
                                env=environment, capture_output=True, text=True)
        assert result.returncode == 0, f"{name} did not run: {result.stderr[-400:]}"
        outputs.append(target)

    collected = {}
    for target in outputs:
        with h5py.File(target, "r") as handle:
            found = {}
            handle["/Calculation"].visititems(
                lambda path, item: found.__setitem__(path, item[()])
                if isinstance(item, h5py.Dataset) else None)
            collected[target] = found

    first, second = (collected[target] for target in outputs)
    assert first, f"{name} wrote nothing under /Calculation"
    assert first.keys() == second.keys()
    for path, values in first.items():
        values, other = np.asarray(values), np.asarray(second[path])
        if not np.issubdtype(values.dtype, np.number):
            continue
        scale = max(np.abs(values).max(), 1.0)
        difference = np.abs(values - other).max()
        assert difference < 1e-12 * scale, (
            f"{name} reran {path} to a different answer: {difference:.3e} "
            f"against a scale of {scale:.3e}; the run is not reproducible")
