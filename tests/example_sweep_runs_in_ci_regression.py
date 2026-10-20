"""Both halves of the example sweep run on every push.

A diverged example survived in the repository because the check that would have
caught it was one somebody had to remember to run. The sweep is only worth
having if it is automatic, so what is asserted here is not that the sweep works
-- other files do that -- but that CI actually runs it, on pushes and on pull
requests, without a filter that could quietly exclude it.

The sweep has two halves and they run in different places:

* The **static** half builds each example's lattice and compares the spectrum
  its own parameters produce against the `spectrum_range` the script declares.
  It executes nothing, takes seconds, and is collected into the ordinary test
  suite (`tests/example_spectrum_range_sweep.py`, plus the hopping-range sweep
  beside it), which CI runs as `pytest tests/`.
* The **behavioural** half runs every example through KITEx and KITE-tools and
  scans the output for non-finite values. That is roughly forty minutes of
  compute, so it is its own job rather than part of the test suite.

The workflow is parsed as text rather than with PyYAML, following
tests/docs_build_strict_regression.py: PyYAML is importable in a development
environment but declared in neither environment.yml nor pyproject.toml, so an
importorskip here would skip in CI -- exactly where the guard is needed, and it
would report as a pass.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / ".github" / "workflows" / "ci.yml"

#: The static half. Collected by pyproject.toml's `python_files` patterns.
STATIC_SWEEPS = ["example_spectrum_range_sweep.py", "shipped_model_hopping_range_sweep.py"]

#: The behavioural half.
RUNNER = "tests/run_all_examples.py"


def _workflow():
    assert WORKFLOW.exists(), f"{WORKFLOW} is missing; CI cannot run anything"
    return WORKFLOW.read_text()


def test_the_workflow_triggers_on_pushes_and_pull_requests():
    text = _workflow()
    trigger = text.split("jobs:")[0]
    assert re.search(r"^\s*push:", trigger, re.M), "the workflow does not run on push"
    assert re.search(r"^\s*pull_request:", trigger, re.M), (
        "the workflow does not run on pull requests, so a change can reach review "
        "without the sweep having run")


def test_the_behavioural_sweep_is_invoked():
    text = _workflow()
    calls = [line.strip() for line in text.splitlines() if RUNNER in line]
    assert calls, (
        f"{RUNNER} is never invoked in {WORKFLOW.name}; the behavioural sweep "
        "depends on someone remembering to run it, which is how a diverged "
        "example survived in the first place")

    for call in calls:
        assert "--only" not in call, (
            f"the sweep is invoked with --only, so it covers a subset: {call!r}")
        assert "--skip" not in call, (
            f"the sweep is invoked with --skip, so some example is excluded: {call!r}")


def test_the_behavioural_sweep_can_fail_the_build():
    """`continue-on-error` or a trailing `|| true` would make it advisory."""
    text = _workflow()
    block = text[text.index("example-sweep:"):]
    nxt = re.search(r"^  [a-z][a-z0-9-]*:$", block[len("example-sweep:"):], re.M)
    if nxt:
        block = block[:len("example-sweep:") + nxt.start()]

    assert RUNNER in block, (
        "the example-sweep job does not invoke the sweep runner")
    assert "continue-on-error" not in block, (
        "the example-sweep job is marked continue-on-error, so a diverged "
        "example would not fail the build")
    assert "|| true" not in block, (
        "the sweep's exit status is discarded, so a failure would not fail the build")


def test_the_test_suite_that_carries_the_static_sweep_runs():
    text = _workflow()
    runs = [line.strip() for line in text.splitlines()
            if re.search(r"pytest\s+tests/", line)]
    assert runs, "CI never runs `pytest tests/`, which carries the static sweep"
    for run in runs:
        assert "-k" not in run.split(), (
            f"pytest is invoked with a -k filter, which could exclude the sweep: {run!r}")


def test_the_static_sweeps_are_collected_by_pytest():
    """A sweep the glob does not match is a file nobody runs."""
    import fnmatch

    config = (REPO / "pyproject.toml").read_text()
    section = config[config.index("[tool.pytest.ini_options]"):]
    patterns = re.findall(r'"([^"]+)"', section[section.index("python_files"):
                                                section.index("\n", section.index("python_files"))])
    assert patterns, "pyproject.toml declares no python_files patterns"

    for sweep in STATIC_SWEEPS:
        assert (REPO / "tests" / sweep).exists(), f"tests/{sweep} is missing"
        assert any(fnmatch.fnmatch(sweep, p) for p in patterns), (
            f"tests/{sweep} matches none of {patterns}, so pytest never collects it")
