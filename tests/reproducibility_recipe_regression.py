"""What a KITEx run is reproducible from: the two seeds *and* `divisions`.

`Src/Simulation/GlobalSimulation.cpp:58-59` mixes each seed with the OpenMP
thread index --

    seed_h = (s0 != 0) * 2654435761 * (thread_id + 8191 * s0)

-- so the realisation is a property of the seeds and the thread count together.
Quoting a seed alone does not let anyone reproduce a figure.

That mixing is *one* reason `divisions` belongs in the recipe, and measurably
not the only one: with the `thread_id` term removed from the expression above,
a run at `divisions=[2,2]` still disagrees with one at `[1,1]`, because the
decomposition partitions the work differently. So the test below asserts that
`divisions` matters, which is what a reproduction recipe needs, and does not
assert *why*.

The thread count is *not* an environment setting. `GlobalSimulation.cpp:48`
calls `omp_set_num_threads(rglobal.n_threads)`, and `n_threads` is the product
of `divisions`, so `OMP_NUM_THREADS` is overridden and has no effect on the
result. That is worth asserting rather than assuming: the natural guess is the
opposite, and a recipe that told people to pin `OMP_NUM_THREADS` would be
pinning the wrong thing while leaving `divisions` free.

So: same seeds and same `divisions` is bit-identical; change either seed, or
change `divisions`, and the numbers move.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

name, nx, sh, sv = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['A', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'A', 'A', -1.0), ([0, 1], 'A', 'A', -1.0))

# Anderson disorder exercises seed_h; the stochastic trace exercises seed_v.
dis = kite.Disorder(lat)
dis.add_disorder('A', 'Uniform', 0.0, 0.4)

cfg = kite.Configuration(divisions=[nx, nx], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=False,
                         precision=1, spectrum_range=[-6, 6],
                         seed_h=sh, seed_v=sv)
calc = kite.Calculation(cfg)
calc.dos(num_points=512, num_moments=64, num_random=2, num_disorder=1)
kite.config_system(lat, cfg, calc, filename=name, disorder=dis)
'''

#: name -> (divisions along each axis, seed_h, seed_v, extra environment)
CASES = {
    "base":   (1, 7, 11, {}),
    "again":  (1, 7, 11, {}),
    "omp4":   (1, 7, 11, {"OMP_NUM_THREADS": "4"}),
    "div2":   (2, 7, 11, {}),
    "seed_h": (1, 8, 11, {}),
    "seed_v": (1, 7, 12, {}),
}


@pytest.fixture(scope="module")
def moments(tmp_path_factory):
    import os
    import h5py
    d = tmp_path_factory.mktemp("reproducibility")
    out = {}
    for name, (nx, sh, sv, extra) in CASES.items():
        h5 = d / f"{name}.h5"
        subprocess.run([sys.executable, "-c", _BUILDER, str(h5), str(nx),
                        str(sh), str(sv)], cwd=d, check=True,
                       capture_output=True, text=True)
        subprocess.run([str(KITEX), str(h5)], cwd=d, check=True,
                       capture_output=True, env=dict(os.environ, **extra))
        with h5py.File(h5, "r") as f:
            out[name] = np.asarray(f["/Calculation/dos/MU"][()]).ravel()
    return out


def _gap(moments, name):
    return np.abs(moments[name] - moments["base"]).max()


def test_the_same_seeds_and_divisions_reproduce_bit_for_bit(moments):
    assert _gap(moments, "again") == 0.0, (
        "two runs with the same seeds and the same divisions disagree, so "
        "nothing in this file can be reproduced from its inputs")


def test_the_thread_environment_does_not_enter(moments):
    """`omp_set_num_threads(prod(divisions))` overrides OMP_NUM_THREADS."""
    assert _gap(moments, "omp4") == 0.0, (
        "OMP_NUM_THREADS changed the result. KITEx calls omp_set_num_threads "
        "with the product of divisions, so if this ever starts mattering the "
        "reproduction recipe is wrong and has to name it too")


@pytest.mark.parametrize("name", ["div2", "seed_h", "seed_v"])
def test_changing_any_of_the_three_changes_the_answer(moments, name):
    """Each of divisions, seed_h and seed_v belongs in the recipe."""
    gap = _gap(moments, name)
    assert gap > 1e-6, (
        f"changing {name} left the moments unchanged (max difference "
        f"{gap:.3e}). Either it does not belong in the reproduction recipe, or "
        "this system is too small to show the difference")
