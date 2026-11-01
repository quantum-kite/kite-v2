"""What the magnetic-field example says about its field must be what it sets.

The field KITE simulates is not the number passed to `kite.Modification`. That
number is quantized to an integer multiple m of the smallest commensurate
field, and the smallest commensurate field is one flux quantum per *column* of
L_y unit cells. So the quantity the Hamiltonian sees is the flux per unit cell,
f = m / L_y, and the total through an L_x by L_y sample is m * L_x quanta -- a
factor of L_x larger than m, which is easy to quote by mistake.

It matters beyond bookkeeping, because f decides the regime. At small f the
continuum Dirac Landau levels are a good description; at f of order a tenth the
magnetic length is down to a lattice constant and the spectrum is Hofstadter's,
with only the zeroth level surviving as a Landau level. A docstring describing
one while the script sets up the other is wrong about the physics, not just
about a number.

The arithmetic here repeats the exporter's own (python/kite/__init__.py, the
`magnetic_field_min` block) rather than importing it, so a change to the
convention shows up as a disagreement instead of being tracked silently.
"""

import contextlib
import importlib.util
import io
import re
from pathlib import Path

import numpy as np
import pytest

import kite

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "dccond_graphene_magnetic_field.py"
SWEEP = REPO / "sweep-data" / "dccond_graphene_magnetic_field"

#: eV*s, as the exporter defines it.
HBAR = 6.58211899e-16
PHI0 = 2 * np.pi * HBAR


class _Captured(Exception):
    pass


def _configuration():
    """Run the example as far as its config_system call, writing nothing."""
    grabbed = {}
    real = kite.config_system

    def spy(lattice, configuration, calculation, **kwargs):
        grabbed.update(lattice=lattice, configuration=configuration, kwargs=kwargs)
        raise _Captured()

    kite.config_system = spy
    try:
        spec = importlib.util.spec_from_file_location("mag_example", EXAMPLE)
        module = importlib.util.module_from_spec(spec)
        with contextlib.redirect_stdout(io.StringIO()):
            spec.loader.exec_module(module)
            if not grabbed:
                try:
                    module.main()
                except _Captured:
                    pass
    finally:
        kite.config_system = real
    assert grabbed, "the example never reached config_system"
    return grabbed


@pytest.fixture(scope="module")
def field():
    grabbed = _configuration()
    config = grabbed["configuration"]
    lattice = grabbed["lattice"]
    modification = grabbed["kwargs"]["modification"]

    lx, ly = (int(n) for n in config.leng)
    vectors = np.asarray(lattice.vectors)[:, :2]
    area = abs(vectors[0, 0] * vectors[1, 1] - vectors[0, 1] * vectors[1, 0]) * 1e-18

    b_min = PHI0 / (ly * area)
    m = int(round(modification.magnetic_field / b_min))

    return {
        "lx": lx, "ly": ly, "m": m,
        "flux_per_cell": m / ly,
        "total_flux": m * lx,
        # l_B^2 = Phi_0/(2 pi B) = A_cell/(2 pi f). `area * 1e18` puts it back
        # in nm^2, so the square root is in nm; dividing by |a1| expresses it in
        # lattice constants, which is what the assertions below compare against.
        # Without that division the two agree only for a lattice whose |a1| is
        # 1 nm, which is true of the example this test drives and of nothing
        # else -- dos_dccond_haldane.py uses a = 0.24595 nm.
        "l_B": float(np.sqrt((area * 1e18) / (2 * np.pi * (m / ly)))
                     / np.linalg.norm(vectors[0])),
    }


def test_the_requested_field_lands_on_an_integer_multiple(field):
    """A requested field between two multiples would be silently rounded."""
    grabbed = _configuration()
    requested = grabbed["kwargs"]["modification"].magnetic_field
    config = grabbed["configuration"]
    ly = int(config.leng[1])
    vectors = np.asarray(grabbed["lattice"].vectors)[:, :2]
    area = abs(vectors[0, 0] * vectors[1, 1] - vectors[0, 1] * vectors[1, 0]) * 1e-18
    exact = requested / (PHI0 / (ly * area))
    assert abs(exact - round(exact)) < 0.02, (
        f"the requested field is {exact:.3f} minimum fields, which rounds to "
        f"{round(exact)}; the example would simulate a field it does not ask for")


def test_the_docstring_states_the_flux_the_script_sets(field):
    """Whatever numbers the prose gives must be the ones the configuration produces.

    Only the flux per cell is required -- it is what sets the regime. The others
    are checked if present, so the docstring can stay short without the check
    quietly passing on nothing.
    """
    text = EXAMPLE.read_text()

    fraction = re.search(r"f = (\d+)/(\d+)|(\d+)/(\d+) of a flux quantum", text)
    assert fraction, "the docstring no longer states the flux per unit cell"
    groups = [g for g in fraction.groups() if g]
    stated = int(groups[0]) / int(groups[1])
    assert stated == pytest.approx(field["flux_per_cell"]), (
        f"the docstring says {stated:.6f} flux quanta per cell; the configuration "
        f"gives {field['m']}/{field['ly']} = {field['flux_per_cell']:.6f}")

    length = re.search(r"magnetic length is ([\d.]+) a|l_B = ([\d.]+) a", text)
    if length:
        value = float(next(g for g in length.groups() if g))
        assert value == pytest.approx(field["l_B"], abs=0.01), (
            f"the docstring says l_B = {value} a; the configuration gives "
            f"{field['l_B']:.3f} a")

    total = re.search(r"(\d+) Phi_0 through the whole sample", text)
    if total:
        assert int(total.group(1)) == field["total_flux"], (
            f"the docstring says {total.group(1)} quanta through the sample; "
            f"the configuration gives m*Lx = {field['total_flux']}")


def test_the_field_is_in_the_lattice_regime_the_docstring_describes(field):
    """The prose describes Hofstadter physics; that has to be the regime set up.

    The continuum Dirac description needs the magnetic length to be long
    compared with the lattice constant. Here it is not, which is why the
    docstring says only the zeroth Landau level survives.
    """
    assert field["l_B"] < 3.0, (
        f"l_B = {field['l_B']:.2f} a is long enough that the continuum Dirac "
        "Landau-level picture may apply again, which is not what the docstring "
        "describes")

    # Where the first Dirac Landau level would fall: E_1 = t sqrt(2 pi sqrt(3) f).
    e1 = np.sqrt(2 * np.pi * np.sqrt(3) * field["flux_per_cell"])
    assert e1 > 0.5, (
        f"the first Dirac Landau level would sit at {e1:.2f} t, inside the linear "
        "region, so several levels would be resolvable and the docstring's "
        "account of a single surviving level would be wrong")


def test_the_measured_plateau_is_one_quantum_per_spin():
    """The plateau the docstring promises, from the example's own output."""
    path = SWEEP / "condDC.dat"
    if not path.exists():
        pytest.skip("needs sweep-data/; run tests/run_all_examples.py first")
    data = np.loadtxt(path)
    mu, sigma = data[:, 0], data[:, 1]

    window = (np.abs(mu) > 0.25) & (np.abs(mu) < 0.6)
    assert window.sum() > 10, "too few points in the plateau window"

    plateau = np.abs(sigma[window]).mean()
    assert plateau == pytest.approx(1.0, rel=0.05), (
        f"the Hall plateau is {plateau:.4f} e^2/h per spin, not the quantized 1 "
        "the docstring states")

    flatness = np.abs(sigma[window]).std() / plateau
    assert flatness < 0.05, f"the plateau is {flatness:.1%} ragged"
