"""`localized_wave_packet` is descoped this release and must stay undocumented.

RF-014 -- the localized wave packet and its lead/sample/transmission machinery --
is out of scope for this version. The implementation is deliberately NOT deleted:
`Calculation.localized_wave_packet()` still exists, KITEx still dispatches
`calc_localized_wavepacket`, and `Src/Simulation/SimulationLocalizedWavePacket.cpp`
still builds. Only the documentation is withdrawn.

Nothing is being hidden that was ever shipped as documented. The feature was
first written up in `3d67f7d` (2026-07-22, "Document the previously-undocumented
Calculation methods..."), and no commit before 2026-07 mentions it anywhere under
`docs/`. Removing those pages restores the state the project was in.

This guards both halves: the pages stay clear, and the code stays present.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

#: Every page that documented it, so a re-introduction is caught where it happened.
PAGES = [
    "docs/documentation/calculation.md",
    "docs/api/kite.md",
    "docs/documentation/examples/time_evolution.md",
]


@pytest.mark.parametrize("page", PAGES)
def test_no_page_documents_the_descoped_capability(page):
    text = (REPO / page).read_text()
    assert "localized_wave_packet" not in text, (
        f"{page} documents localized_wave_packet, which is descoped this "
        "release per RF-014. If it is being brought back into scope, that is a "
        "PM decision and this test should go with it")


def test_no_dangling_reference_survives_anywhere_in_the_docs():
    """A link definition left behind breaks the build only sometimes."""
    offenders = [p.relative_to(REPO) for p in (REPO / "docs").rglob("*.md")
                 if "localized_wave_packet" in p.read_text(errors="ignore")]
    assert not offenders, f"localized_wave_packet is still referenced in {offenders}"


def test_the_implementation_is_still_present():
    """Descoped means undocumented, not deleted. The code must survive."""
    api = (REPO / "python" / "kite" / "__init__.py").read_text()
    assert re.search(r"def localized_wave_packet\s*\(", api), (
        "Calculation.localized_wave_packet() has been deleted. RF-014 descopes "
        "the documentation only; removing the implementation is a different "
        "decision and was not taken")

    driver = REPO / "Src" / "Simulation" / "SimulationLocalizedWavePacket.cpp"
    assert driver.exists(), f"{driver.name} has been deleted; the implementation stays"

    dispatch = (REPO / "Src" / "Simulation" / "GlobalSimulation.cpp").read_text()
    assert "calc_localized_wavepacket" in dispatch, (
        "KITEx no longer dispatches calc_localized_wavepacket")
