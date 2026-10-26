"""KITE-tools must reconstruct `custom_one_local` moments, and be wired to do it.

Until this reader existed, `custom_one_local` handed the user bare Chebyshev
moments at `/Calculation/CustomOneLocal/Gamma` with no shipped way to turn them
into a spectral density. That combination -- an entry point with no reader and
no example -- is why BUG-009 survived as long as it did: nothing exercised it
end to end.

Two kinds of check here.

*Registration* is static and mechanical. `shell_input` finds a flag by indexing
`keys_pos` with the flag's position in `valid_keys`, a hand-written integer. Get
it wrong, or insert a key ahead of it, and the flag is silently read from the
wrong place -- the same class of defect as BUG-011's mismatched sentinel. Both
are grep-level and both are checked.

*Reconstruction* runs the real binaries. The reference is an independent Jackson
sum written here in NumPy, against the moments as stored, so the test re-derives
the answer rather than re-running the implementation. The tolerance is set by the
six significant digits KITE-tools writes, not by what happens to pass.
"""

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"
PARSE = REPO / "tools" / "Src" / "Tools" / "parse_input.cpp"


# --------------------------------------------------------------------------
# registration
# --------------------------------------------------------------------------

def _valid_keys():
    text = PARSE.read_text()
    match = re.search(r"valid_keys\s*=\s*std::vector<std::string>\{(.*?)\};", text, re.S)
    assert match, "parse_input.cpp no longer declares valid_keys"
    return re.findall(r'"(--[A-Za-z0-9_]+)"', match.group(1))


def test_the_flag_is_registered():
    assert "--CustomOneLocal" in _valid_keys(), (
        "--CustomOneLocal is gone from valid_keys, so KITE-tools cannot see the flag")


def test_every_parser_indexes_its_own_key():
    """`keys_pos.at(N)` must match where that flag actually sits in valid_keys.

    The index is written by hand in each parse_* function. Inserting a key
    earlier in the list shifts every later one, and nothing complains: the
    parser just reads another flag's arguments.
    """
    keys = _valid_keys()
    text = PARSE.read_text()
    mismatched = []
    checked = 0
    for name, chunk in re.findall(r"void shell_input::parse_(\w+)\(int argc[^\n]*\n(.*?)\n\}", text, re.S):
        # Two forms in the file: keys_pos.at(6) directly, or `int j = 2;`
        # followed by keys_pos.at(j).
        index = re.search(r"keys_pos\.at\((\d+)\)", chunk)
        if index:
            used = int(index.group(1))
        else:
            named = re.search(r"int\s+j\s*=\s*(\d+);[\s\S]*?keys_pos\.at\(j\)", chunk)
            if not named:
                continue
            used = int(named.group(1))
        flag = "--" + name if "--" + name in keys else None
        if flag is None:
            # parse_lDOS -> --LDOS, parse_CondDC -> --CondDC, ...
            lowered = {k.lower(): k for k in keys}
            flag = lowered.get(("--" + name).lower())
        assert flag is not None, f"parse_{name} does not correspond to any key in valid_keys"
        checked += 1
        if keys.index(flag) != used:
            mismatched.append(f"parse_{name} uses keys_pos.at({used}) but {flag} is at {keys.index(flag)}")
    assert checked == len(keys), (
        f"only {checked} of {len(keys)} flags have a parser that indexes keys_pos; "
        "the scan missed some, so a mismatch could hide in the ones it skipped")
    assert not mismatched, (
        "a parser reads the wrong flag's arguments: " + "; ".join(mismatched))


def test_the_sentinels_agree_with_the_defaults():
    """Each override must compare against the value the parser assigns.

    BUG-011 is this check failing in conductivity_2order; the reader added here
    must not repeat it.
    """
    text = PARSE.read_text()
    chunk = re.search(r"void shell_input::parse_CustomOneLocal\(int argc[^\n]*\n(.*?)\n\}", text, re.S)
    assert chunk, "parse_CustomOneLocal is gone"
    assigned = dict(re.findall(r"(CustomOneLocal_\w+)\s*=\s*(-?[\d.]+);", chunk.group(1)))

    reader = (REPO / "tools" / "Src" / "Custom" / "customonelocal.cpp").read_text()
    compared = dict(re.findall(r"variables\.(CustomOneLocal_\w+)\s*!=\s*(-?[\d.]+)", reader))

    assert compared, "customonelocal.cpp compares against no sentinel at all"
    bad = sorted(f for f in compared if f in assigned and compared[f] != assigned[f])
    assert not bad, (
        f"{bad} are compared against a different value than parse_CustomOneLocal "
        "assigns, so those flags are ignored or the sentinel is used as a value")


# --------------------------------------------------------------------------
# reconstruction
# --------------------------------------------------------------------------

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt
from kite import custom

M = 256
lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
lat.add_sublattices(['A', [0, 0], +1.0], ['B', [0.5, 0.5], -1.0])
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([-1, 0], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0), ([-1, -1], 'A', 'B', -1.0))
cfg = kite.Configuration(divisions=[1, 1], length=[16, 16],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-8, 8])
calc = kite.Calculation(cfg)
calc.add_orbital_index('A', 0)
calc.add_orbital_index('B', 1)
calc.add_orbital_coupling('A', 'A', 0.5, 'l0')
calc.add_orbital_coupling('B', 'B', -0.5, 'l0')
calc.custom_one_local(stream_=custom.Vertex(M, [[1.0, 'l0']]),
                      energy_=np.array([0.0]), position_=[[3, 5]] * 2,
                      sublattice_=['A', 'B'])
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _jackson(M):
    n = np.arange(M)
    Mp = M + 1
    return ((Mp - n) * np.cos(np.pi * n / Mp)
            + np.sin(np.pi * n / Mp) / np.tan(np.pi / Mp)) / Mp


def _reconstruct(gamma, energy, scale, shift):
    """The Jackson-damped Chebyshev sum, written out rather than called.

    delta(m, E) carries the 2/(pi*sqrt(1-Etilde^2)) weight; the 1/(1+delta_m0)
    cancels the leading 2 for m=0. No 1/N: the seeding is a single basis ket.
    """
    moments = gamma.shape[0]
    scaled = (energy - shift) / scale
    theta = np.arccos(np.clip(scaled, -1.0, 1.0))
    damping = _jackson(moments)
    out = np.zeros((len(energy), gamma.shape[1]))
    for m in range(moments):
        weight = (2.0 / np.pi) * np.cos(m * theta) / np.sqrt(1.0 - scaled ** 2)
        out += np.outer(weight * damping[m] / (1.0 + (m == 0)) / scale, gamma[m].real)
    return out


@pytest.mark.skipif(not (KITEX.exists() and KITE_TOOLS.exists()),
                    reason="needs a built KITEx and KITE-tools")
def test_the_reader_reproduces_an_independent_reconstruction(tmp_path):
    import h5py

    h5 = tmp_path / "col.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(h5)], cwd=tmp_path,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), h5.name], cwd=tmp_path, check=True, capture_output=True)
    dat = tmp_path / "out.dat"
    subprocess.run([str(KITE_TOOLS), h5.name, "--CustomOneLocal", "-X",
                    "-E", "-7.5", "7.5", "300", "-N", dat.name],
                   cwd=tmp_path, check=True, capture_output=True)
    assert dat.exists(), "--CustomOneLocal produced no output file"

    with h5py.File(h5, "r") as handle:
        gamma = np.asarray(handle["/Calculation/CustomOneLocal/Gamma"]).T
        scale = float(np.asarray(handle["/EnergyScale"]).ravel()[0])
        shift = float(np.asarray(handle["/EnergyShift"]).ravel()[0])

    table = np.loadtxt(dat)
    energy, produced = table[:, 0], table[:, 1:]
    assert produced.shape[1] == 2, (
        f"expected one column per probed orbital, got {produced.shape[1]}")

    expected = _reconstruct(gamma, energy, scale, shift)
    scale_of_signal = np.abs(expected).max()
    assert scale_of_signal > 0.1, "the reference reconstruction is empty"
    # KITE-tools writes six significant digits, so the floor is relative to the
    # largest value on the line, not to each entry.
    assert np.abs(produced - expected).max() < 1e-4 * scale_of_signal, (
        "KITE-tools disagrees with a direct Jackson sum over the stored moments")


@pytest.mark.skipif(not (KITEX.exists() and KITE_TOOLS.exists()),
                    reason="needs a built KITEx and KITE-tools")
def test_a_diagonal_vertex_scales_the_local_density(tmp_path):
    """A|j> = A_jj|j> for a diagonal vertex, so the two sublattices must come
    out as exact opposites of each other on this symmetric lattice -- the A and
    B sites differ only by the sign of their on-site energy and of s_z."""
    import h5py

    h5 = tmp_path / "col.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(h5)], cwd=tmp_path,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), h5.name], cwd=tmp_path, check=True, capture_output=True)

    with h5py.File(h5, "r") as handle:
        gamma = np.asarray(handle["/Calculation/CustomOneLocal/Gamma"])

    assert gamma.shape[0] == 2, "expected one row of moments per probed orbital"
    assert np.abs(gamma.imag).max() == 0.0, (
        "a Hermitian diagonal vertex must give real moments")
    assert gamma[0, 0].real == pytest.approx(+0.5), (
        f"mu_0 on A is {gamma[0, 0].real}, not the +1/2 that s_z times "
        "normalisation requires")
    assert gamma[1, 0].real == pytest.approx(-0.5)
    # mu_n is even in the on-site sign for even n and odd for odd n, so the two
    # columns are related by that sign times the vertex's own sign: both flip,
    # giving equality on odd n and opposite signs on even n.
    even = np.abs(gamma[0, ::2].real + gamma[1, ::2].real).max()
    assert even < 1e-12, (
        f"the A and B columns are not opposite on even moments (max {even:.2e})")
