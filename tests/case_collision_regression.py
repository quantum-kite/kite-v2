"""No two tracked paths may differ only by case.

macOS (APFS) and Windows are case-insensitive by default, so two paths that
differ only in case collapse into one on checkout and the repository cannot be
cloned intact. This is not hypothetical: a `Src/` and `src/` pair made exactly
that happen, and `src/` was renamed to `python/` to resolve it. Nothing stops a
recurrence, and the symptom appears only on a contributor's machine, never on
the Linux CI runner that would otherwise be the first to notice.

Runs standalone with no dependencies (`python tests/case_collision_regression.py`)
so CI can check it without building the project, and as part of the suite.
"""

import subprocess
import sys
from collections import defaultdict
from pathlib import Path

def repo_root():
    """Ask git rather than counting directories up from __file__.

    The checker is also run from a throwaway repository when verifying that it
    actually fails, where it does not sit in tests/.
    """
    out = subprocess.run(["git", "rev-parse", "--show-toplevel"],
                         cwd=str(Path(__file__).resolve().parent),
                         capture_output=True, text=True, check=True)
    return Path(out.stdout.strip())


def tracked_paths():
    out = subprocess.run(["git", "ls-files", "-z"], cwd=str(repo_root()),
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.split("\0") if p]


def find_collisions(paths):
    """Group paths that are equal once lowercased, keeping only real clashes.

    Compared per path component rather than on the whole string, because the
    collision is created by the directory entry: `Src/a.cpp` and `src/b.cpp` do
    not collide as strings but cannot coexist, since `Src` and `src` are one
    directory on a case-insensitive filesystem.

    Returns a sorted list of sorted tuples. Sorting matters twice over: the
    failure message is the whole value of this check, so it has to be stable,
    and `sorted()` over sets would silently do nothing -- sets compare by
    subset, which is a partial order, so Python's stable sort leaves them in
    whatever order the dict happened to produce.

    A collision is reported once, at the shallowest level it occurs. One
    renamed directory with a hundred files under it is one problem, and listing
    it a hundred times buries it.
    """
    seen = defaultdict(set)
    for path in paths:
        parts = path.split("/")
        for depth in range(len(parts)):
            key = "/".join(p.lower() for p in parts[:depth + 1])
            seen[key].add("/".join(parts[:depth + 1]))

    collisions = {}
    for key in sorted(seen):
        if len(seen[key]) < 2:
            continue
        if any(key.startswith(shallower + "/") for shallower in collisions):
            continue                      # already reported by its parent
        collisions[key] = tuple(sorted(seen[key]))
    return sorted(collisions.values())


def test_no_tracked_path_differs_only_by_case():
    collisions = find_collisions(tracked_paths())
    assert not collisions, (
        "these tracked paths differ only by case and cannot coexist on a "
        "case-insensitive filesystem:\n" +
        "\n".join("  " + " vs ".join(sorted(group)) for group in collisions))


def test_the_check_detects_a_collision():
    """The failure this guard exists for: the original Src/ vs src/ pair."""
    collisions = find_collisions([
        "Src/Simulation/Global.cpp",
        "src/kite/__init__.py",
        "docs/index.md",
    ])
    assert collisions == [("Src", "src")], collisions


def test_a_difference_below_the_top_level_is_also_caught():
    collisions = find_collisions(["a/Tools/x.cpp", "a/tools/y.cpp"])
    assert collisions == [("a/Tools", "a/tools")], collisions


def test_one_renamed_directory_is_reported_once():
    """Not once per file beneath it -- the failure message has to stay readable."""
    collisions = find_collisions([f"Src/f{i}.cpp" for i in range(5)] +
                                 [f"src/f{i}.py" for i in range(5)])
    assert collisions == [("Src", "src")], collisions


def test_the_process_exits_nonzero_on_a_colliding_index(tmp_path):
    """The acceptance criterion, end to end: the command CI actually runs.

    The colliding entry goes in through `git update-index --cacheinfo` rather
    than by creating the files, because on a case-insensitive filesystem --
    which is the whole point of this guard -- writing `Src/x` and `src/y` puts
    both under whichever directory was created first, and the collision never
    reaches the index.
    """
    def git(*args):
        return subprocess.run(["git", *args], cwd=str(tmp_path),
                              capture_output=True, text=True, check=True)

    git("init", "-q", "-b", "main")
    (tmp_path / "a.txt").write_text("x\n")
    blob = git("hash-object", "-w", "a.txt").stdout.strip()
    git("update-index", "--add", "--cacheinfo", f"100644,{blob},Src/Global.cpp")
    git("update-index", "--add", "--cacheinfo", f"100644,{blob},src/__init__.py")

    checker = tmp_path / "check.py"
    checker.write_text(Path(__file__).read_text())
    proc = subprocess.run([sys.executable, str(checker)],
                          cwd=str(tmp_path), capture_output=True, text=True)
    assert proc.returncode == 1, (
        f"the guard passed on a colliding index:\n{proc.stdout}\n{proc.stderr}")
    assert "Src vs src" in proc.stdout, proc.stdout


def test_the_process_exits_zero_on_this_repository():
    """Pins main()'s success path, so a broken checker cannot pass by default."""
    proc = subprocess.run([sys.executable, str(Path(__file__).resolve())],
                          capture_output=True, text=True)
    assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"


def test_paths_that_merely_share_a_prefix_are_not_collisions():
    assert find_collisions(["python/kite/lattice.py", "python/kite/custom.py"]) == []


def main():
    collisions = find_collisions(tracked_paths())
    if collisions:
        print("case-colliding tracked paths (cannot coexist on macOS or Windows):")
        for group in collisions:
            print("   " + " vs ".join(sorted(group)))
        return 1
    print("no case-colliding tracked paths")
    return 0


if __name__ == "__main__":
    sys.exit(main())
