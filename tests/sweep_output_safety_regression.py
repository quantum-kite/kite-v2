"""The sweep must not delete what it did not create, nor log who ran it.

Two properties of `tests/run_all_examples.py`, both about `--outdir`, which is
taken from the command line and used for destructive work.

*Deletion.* Every example's working directory is removed with `shutil.rmtree`
before the run. `--outdir` was accepted unchecked, so `--outdir /etc` would have
deleted `/etc/<name>` for every example name that happened to exist. It is a
developer tool invoked by hand and nobody is attacking it, but a destructive
operation should not be one typo away. The sweep now marks its own output root
and refuses any directory that has content it did not put there.

*Disclosure.* `run.log` recorded the command line and the child's output
verbatim, both carrying `/Users/<name>/...`. The default output directory is
ignored by git, but `--outdir` is free, so the log can be written anywhere --
including into a public repository. This is a public scientific package, and the
name of whoever ran the sweep is not part of it.
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


# --- deletion -----------------------------------------------------------

def test_a_fresh_directory_is_claimed_and_marked(runner, tmp_path):
    root = runner.claim_output_root(tmp_path / "fresh")
    assert (root / runner.MARKER).exists(), "the sweep must mark what it owns"


def test_its_own_root_can_be_reclaimed(runner, tmp_path):
    """A second run over the same output directory is the normal case."""
    first = runner.claim_output_root(tmp_path / "again")
    assert runner.claim_output_root(tmp_path / "again") == first


def test_a_directory_it_did_not_create_is_refused(runner, tmp_path):
    occupied = tmp_path / "someone_elses"
    occupied.mkdir()
    (occupied / "important.txt").write_text("not the sweep's")

    with pytest.raises(SystemExit) as exc:
        runner.claim_output_root(occupied)

    assert runner.MARKER in str(exc.value), "the refusal should say what is missing"
    assert (occupied / "important.txt").exists(), "it must refuse, not empty the directory"


def test_an_empty_directory_is_allowed(runner, tmp_path):
    """Nothing to lose, so pointing at an empty directory is not an error."""
    empty = tmp_path / "empty"
    empty.mkdir()
    assert runner.claim_output_root(empty) == empty


def test_a_file_in_place_of_a_directory_is_refused(runner, tmp_path):
    notadir = tmp_path / "a_file"
    notadir.write_text("x")
    with pytest.raises(SystemExit):
        runner.claim_output_root(notadir)


# --- disclosure ---------------------------------------------------------

def test_the_repository_path_is_redacted(runner):
    assert str(REPO) not in runner.redact(f"$ python {REPO}/examples/dos_graphene.py")


def test_the_home_directory_is_redacted(runner):
    home = str(Path.home())
    assert home not in runner.redact(f"Traceback: {home}/miniforge3/envs/kite/bin/python")


def test_the_repository_is_redacted_before_the_home_that_contains_it(runner):
    """Otherwise the repo path becomes <home>/dev/... and stays machine-specific."""
    out = runner.redact(f"{REPO}/examples/x.py")
    assert out.startswith("<repo>"), out
