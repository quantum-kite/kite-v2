"""The sweep's exit status distinguishes a broken example from one it skipped.

`tests/run_all_examples.py` is meant to be a gate: non-zero means something in
the repository is wrong. That only works if non-zero means exactly that.

It reports three different things. An example that failed or produced non-finite
output is a defect. A processor it could not pair to a producer is neither --
pairing is by the HDF5 filename a script happens to mention in its source, so a
processor taking that name as an argument has nothing to match, and it is simply
not run. Counting the second as a failure would make the gate red for a
limitation of its own heuristic; ignoring it would hide a processor that stopped
being checked. So the known ones are listed, and a new one is a failure.

The verdict is exercised against synthetic reports rather than by running the
sweep, which takes an hour and a half.
"""

import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "tests" / "run_all_examples.py"


@pytest.fixture(scope="module")
def runner():
    spec = importlib.util.spec_from_file_location("sweep_runner", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _verdict(runner, results):
    """The runner's own decision, not a copy of it."""
    bad, _unpaired, surprises = runner.verdict(results)
    return 1 if (bad or surprises) else 0


def test_a_clean_run_passes(runner):
    results = {"dos_graphene.py": {"status": "ok"},
               "band_structure_graphene.py": {"status": "ok (no KITEx stage)"}}
    assert _verdict(runner, results) == 0


def test_a_diverged_example_fails(runner):
    """Non-finite output is the signature the sweep exists to catch."""
    results = {"dos_graphene.py": {"status": "ok"},
               "hbn_optcond2_vacancies.py": {
                   "status": "diverged",
                   "diverged": ["nonlinear_cond.dat: 1024 non-finite values"]}}
    assert _verdict(runner, results) == 1


def test_a_failed_stage_fails(runner):
    results = {"dos_graphene.py": {"status": "ok"},
               "afm_zigzag_ribbon_process.py": {"status": "failed",
                                                "error": "FileNotFoundError"}}
    assert _verdict(runner, results) == 1


def test_a_known_unpaired_processor_does_not_fail(runner):
    """It was not run, but that is the pairing's limit, not a defect."""
    assert runner.KNOWN_UNPAIRED, "the exception list is empty"
    name = sorted(runner.KNOWN_UNPAIRED)[0]
    results = {"dos_graphene.py": {"status": "ok"},
               name: {"status": "unpaired", "kind": "processor"}}
    assert _verdict(runner, results) == 0


def test_a_processor_that_becomes_unpaired_does_fail(runner):
    """Otherwise a script could stop being checked and nothing would say so."""
    results = {"dos_graphene.py": {"status": "ok"},
               "some_new_process.py": {"status": "unpaired", "kind": "processor"}}
    assert _verdict(runner, results) == 1


def test_every_documented_exception_carries_a_reason(runner):
    for name, reason in runner.KNOWN_UNPAIRED.items():
        assert reason and reason.strip(), f"{name} is excused with no reason given"


def test_the_exceptions_name_scripts_that_ship(runner):
    """A stale exception silently excuses nothing, and hides that it is stale."""
    present = {p.name for p in (REPO / "examples").glob("*.py")}
    missing = sorted(n for n in runner.KNOWN_UNPAIRED if n not in present)
    assert not missing, (
        f"KNOWN_UNPAIRED names scripts that no longer ship: {missing}")
