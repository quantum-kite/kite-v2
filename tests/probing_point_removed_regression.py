"""`probing_point` is gone from the API, the output and the documentation.

It selected the site whose complex amplitude KITEx reports as `Id` -- the value
at one thread's domain centre, an artifact of how the lattice is split across
threads rather than a property of the wave packet. Nothing in C++ ever read the
parameter back, so removing it touches only the Python layer and the docs; the
`Id` dataset itself is still written by KITEx and is simply not surfaced.
"""

import inspect
from pathlib import Path

import pytest

import kite

REPO = Path(__file__).resolve().parents[1]


def test_the_parameter_is_not_in_the_signature():
    params = inspect.signature(kite.Calculation.gaussian_wave_packet).parameters
    assert "probing_point" not in params


def test_the_exported_file_carries_no_probing_point(tmp_path, monkeypatch):
    """Built from the shipped example, so this is the real export path."""
    h5py = pytest.importorskip("h5py")
    import sys
    sys.path.insert(0, str(REPO / "examples"))
    import spin_precession_simple as example

    monkeypatch.chdir(tmp_path)
    name = example.main(num_points=4, num_moments=8)

    with h5py.File(name, "r") as handle:
        fields = set(handle["Calculation/gaussian_wave_packet"].keys())
    assert "ProbingPoint" not in fields, f"still exported: {sorted(fields)}"


@pytest.mark.parametrize("page", [
    "docs/api/kite.md",
    "docs/documentation/editing_hdf_files.md",
])
def test_the_documentation_does_not_mention_it(page):
    text = (REPO / page).read_text()
    assert "probing_point" not in text and "ProbingPoint" not in text, (
        f"{page} still documents a parameter that no longer exists")
