"""The two Kubo-Bastin reconstructions must share one default integration grid.

Both integrate over energy on a uniform grid spanning [-0.99, 0.99] in KPM
reduced units, and both are reconstructions of the same kind of object, so a
default that drifts apart between them would make `conductivity_dc` and
`custom_two` disagree for a reason no user could see.

How fine that grid needs to be depends on the Green's-function broadening and
on the number of Chebyshev moments, so no single constant is right everywhere;
`-E` exists for work that needs more. The default only has to be sane out of
the box, and this file only pins that the two agree on it.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

SOURCES = {
    "custom_two": REPO / "tools" / "Src" / "Custom" / "customtwo.cpp",
    "conductivity_dc": REPO / "tools" / "Src" / "CondDC" / "conductivity_dc.cpp",
}


def _default_grid(path):
    match = re.search(r"NEnergies\s*=\s*(\d+)\s*;", path.read_text())
    assert match, f"{path.name} no longer sets a default NEnergies"
    return int(match.group(1))


def test_both_reconstructions_share_one_default():
    grids = {name: _default_grid(path) for name, path in SOURCES.items()}
    assert len(set(grids.values())) == 1, (
        f"the two reconstructions default to different integration grids: {grids}")


def test_the_default_grid_is_odd():
    """Simpson's rule needs an odd point count.

    Both files bump an even value to odd before building the grid, so an even
    default is corrected rather than wrong -- but it means the printed banner
    and the grid actually used disagree, which is worth not doing.
    """
    for name, path in SOURCES.items():
        points = _default_grid(path)
        assert points % 2 == 1, (
            f"{name} defaults to {points} integration points, which is even; "
            "the code will silently use a different number")
