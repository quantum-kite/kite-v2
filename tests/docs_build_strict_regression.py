"""Every documentation build must use --strict, and the deploy must wait for it.

`mkdocs build --strict` turns every warning -- a broken internal link, a missing
nav target, an unresolved reference -- into a failure, so documentation rot
fails CI instead of reaching the published site. Without it mkdocs prints the
warning and exits 0.

Two halves, and the second is easy to lose. `mkdocs gh-deploy` has no --strict
equivalent, so unless the deploy job waits for the strict build the two race on
a push and a site that failed the check is published anyway -- which makes the
strict build advisory rather than a gate.

Both pipelines are surveyed, not just the GitHub one: `.gitlab-ci.yml`
publishes to GitLab Pages off the `docs` branch and is a separate path to the
same kind of mistake.

The workflows are parsed as text rather than with PyYAML, following
tests/documented_conda_install_regression.py. PyYAML is importable in a
development environment but is declared in neither environment.yml nor
pyproject.toml, so a `pytest.importorskip("yaml")` here would skip this whole
file in CI -- where `build-and-test` runs pytest without mkdocs, and so without
mkdocs's PyYAML. A guard that skips exactly where it is needed is worse than no
guard, because it reports as a pass.

The build itself is exercised by CI's `check-docs` job, which installs
mkdocs-material and runs the real command on every push and pull request. There
is deliberately no test duplicating that here: pytest runs in a job that has no
mkdocs, so such a test could only ever skip.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PIPELINES = [REPO / ".github" / "workflows" / "ci.yml",
             REPO / ".gitlab-ci.yml"]

MKDOCS_BUILD = re.compile(r"^\s*-?\s*(?:run:\s*)?(mkdocs\s+build\b.*)$")


def _build_commands():
    """Every `mkdocs build` invocation across every pipeline file present."""
    found = []
    for path in PIPELINES:
        if not path.exists():
            continue
        for number, line in enumerate(path.read_text().splitlines(), 1):
            match = MKDOCS_BUILD.search(line)
            if match:
                found.append((path.name, number, match.group(1).strip()))
    return found


def _job_block(text, job):
    """The lines of one top-level job, by indentation rather than by YAML."""
    lines = text.splitlines()
    try:
        start = next(i for i, l in enumerate(lines)
                     if re.match(rf"^  {re.escape(job)}:\s*$", l))
    except StopIteration:
        return None
    block = []
    for line in lines[start + 1:]:
        if line.strip() and not line.startswith("    "):
            break
        block.append(line)
    return block


def test_every_documentation_build_is_strict():
    """Covers both pipelines; without --strict mkdocs warns and exits 0 anyway."""
    commands = _build_commands()
    assert commands, "no pipeline runs `mkdocs build` at all; has CI changed?"
    lax = [(f, n, c) for f, n, c in commands if "--strict" not in c]
    assert not lax, (
        "these documentation builds do not use --strict, so a broken link "
        f"would be published rather than failing: {lax}")


def test_the_github_deploy_waits_for_the_strict_build():
    """Otherwise the strict check is advisory, not a gate.

    `mkdocs gh-deploy` takes no --strict flag, so the only thing stopping a
    failed strict build from being published is this dependency. The `if:` on
    the deploy job carries no always()/success() call, so Actions applies an
    implicit success() to its needs -- a failed check-docs skips the deploy.
    """
    text = (REPO / ".github" / "workflows" / "ci.yml").read_text()
    assert "gh-deploy" in text, "no job deploys the site; has the workflow changed?"
    assert "mkdocs build --strict" in text, "no job runs the strict build"

    deploy = _job_block(text, "deploy-docs")
    assert deploy is not None, "the deploy-docs job is gone or was renamed"
    assert any("gh-deploy" in line for line in deploy), (
        "deploy-docs no longer deploys; find the job that does and gate that")

    needs = [line for line in deploy if re.match(r"\s*needs:", line)]
    assert needs, (
        "deploy-docs publishes the site without waiting for anything; a site "
        "that failed the strict build would be deployed regardless")
    assert any("check-docs" in line for line in needs), (
        f"deploy-docs waits on {needs}, but not on the strict build")

    guard = [line for line in deploy if re.match(r"\s*if:", line)]
    assert not any("always()" in line for line in guard), (
        "deploy-docs uses always(), which overrides the implicit success() on "
        "needs and re-opens the race the dependency closes")
