"""The two direction tables must agree with each other and with Python.

A conductivity direction travels as an integer: Python maps `'yz'` to 6, KITEx
turns 6 back into a string with `num2str2`, and KITE-tools does the same with its
own copy of that function. Three tables, in three files, that nothing tied
together.

In KITEx the string is not only a label. `SimulationCondDC.cpp` does

    std::string dir(num2str2(direction));
    std::string dirc = dir.substr(0,1) + "," + dir.substr(1,2);
    Gamma2D(..., process_string(dirc), ".../Gamma" + dir);

so it chooses the velocity pair *and* names the dataset. `num2str2` returned
`"yx"` for both case 5 and case 6, so `direction='yz'` computed yx and stored it
as `Gammayx`; KITE-tools, whose copy had `case 6: "yz"` correct, then looked for
`Gammayz`, failed, and told the user to make sure KITEx had processed the file --
which they had.

Checked as a whole table rather than at case 6, because the defect was a repeated
line, and the next repeated line will be somewhere else.
"""

import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

KITEX_SRC = REPO / "Src" / "Tools" / "queue.cpp"
TOOLS_SRC = REPO / "tools" / "Src" / "Tools" / "functions.cpp"
PY_SRC = REPO / "python" / "kite" / "__init__.py"

_CASE = re.compile(r'case\s+(\d+)\s*:\s*(?://[^\n]*\n\s*)*dir\s*=\s*"(\w+)"')


def _cpp_table(path):
    text = path.read_text()
    i = text.index("num2str2(int dir_num)")
    body = text[i:text.index("default:", i)]
    return {int(n): d for n, d in _CASE.findall(body)}


def _python_table():
    line = next(l for l in PY_SRC.read_text().splitlines()
                if "_avail_dir_full" in l and "{" in l)
    return {int(v): k for k, v in
            re.findall(r"'(\w+)':\s*(\d+)", line)}


def test_the_tables_were_parsed():
    """A regex that stopped matching would make every assertion below vacuous."""
    for name, table in (("KITEx", _cpp_table(KITEX_SRC)),
                        ("KITE-tools", _cpp_table(TOOLS_SRC)),
                        ("Python", _python_table())):
        assert len(table) == 9, f"{name} table has {len(table)} entries, expected 9"


def test_no_direction_string_is_repeated():
    """The defect was one string on two cases, which is what made it invisible."""
    for name, path in (("KITEx", KITEX_SRC), ("KITE-tools", TOOLS_SRC)):
        table = _cpp_table(path)
        seen = {}
        for code, text in sorted(table.items()):
            assert text not in seen, (
                f"{name}'s num2str2 returns '{text}' for both case {seen[text]} "
                f"and case {code}; one of them is picking the wrong velocity pair "
                "and writing under the wrong dataset name")
            seen[text] = code


def test_all_three_tables_agree():
    """Python chooses the integer; both C++ sides must read it the same way."""
    py, kitex, tools = _python_table(), _cpp_table(KITEX_SRC), _cpp_table(TOOLS_SRC)
    assert kitex == tools, (
        "KITEx and KITE-tools disagree on the direction codes: "
        f"{ {k: (kitex.get(k), tools.get(k)) for k in set(kitex) | set(tools) if kitex.get(k) != tools.get(k)} }")
    assert kitex == py, (
        "the C++ tables disagree with python/kite's _avail_dir_full: "
        f"{ {k: (kitex.get(k), py.get(k)) for k in set(kitex) | set(py) if kitex.get(k) != py.get(k)} }")


# --- and the behaviour the mismatch produced ---------------------------------

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt
lat = latt.Lattice(a1=[1, 0, 0], a2=[0, 1, 0], a3=[0, 0, 1])
lat.add_sublattices(('A', [0, 0, 0], 0.0))
lat.add_hoppings(([1, 0, 0], 'A', 'A', -1.0), ([0, 1, 0], 'A', 'A', -1.0),
                 ([0, 0, 1], 'A', 'A', -1.0))
cfg = kite.Configuration(divisions=[1, 1, 1], length=[16, 16, 16],
                         boundaries=['periodic'] * 3, is_complex=True,
                         precision=1, spectrum_range=[-7, 7],
                         seed_h=20261031, seed_v=20261031)
calc = kite.Calculation(cfg)
calc.conductivity_dc(num_points=50, num_moments=32, num_random=1,
                     num_disorder=1, direction=sys.argv[2], temperature=0.1)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


@pytest.mark.skipif(not (KITEX.exists() and KITE_TOOLS.exists()),
                    reason="needs a built KITEx and KITE-tools")
@pytest.mark.parametrize("direction", ["yz", "yx", "xz"])
def test_the_requested_direction_is_the_one_that_is_stored(tmp_path, direction):
    """KITEx must write Gamma<direction>, and KITE-tools must find it."""
    name = f"{direction}.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, name, direction],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), name], cwd=tmp_path, check=True, capture_output=True)

    with h5py.File(tmp_path / name, "r") as f:
        gammas = [k for k in f["Calculation/conductivity_dc"] if k.startswith("Gamma")]
    assert gammas == [f"Gamma{direction}"], (
        f"direction='{direction}' produced {gammas}; KITEx resolved the code to a "
        "different pair of velocities than the one asked for")

    out = subprocess.run([str(KITE_TOOLS), name, "--CondDC", "-N", "c.dat"],
                         cwd=tmp_path, check=True, capture_output=True, text=True)
    assert "was not found" not in out.stdout, (
        f"KITE-tools cannot find the moments for direction='{direction}':\n{out.stdout}")
    data = np.loadtxt(tmp_path / "c.dat")
    assert len(data) == 50 and np.isfinite(data).all()
