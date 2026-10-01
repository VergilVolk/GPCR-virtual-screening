"""Output-path policy for new-molecule frozen runs.

The frozen `path_guard.guard_write_path` intentionally restricts writes to the
three historical v02 roots, so it cannot express a new molecule's outputs.  This
module supplies the equivalent policy for prospective runs:

  * the protected-root list is IMPORTED from the frozen guard, never copied, so
    it can never drift;
  * every output must live under a root declared in the run spec;
  * no output may fall inside any protected historical root;
  * no existing file may be overwritten.
"""

from __future__ import annotations

import os
from pathlib import Path

from project.pacer_fkg_v02.frozen_runner.engines import PROTECTED_RELATIVE_ROOTS, UnsafeWritePath


def _resolved(path: os.PathLike[str] | str) -> Path:
    return Path(path).expanduser().resolve(strict=False)


def _within(path: Path, root: Path) -> bool:
    try:
        return os.path.normcase(os.path.commonpath((str(path), str(root)))) == os.path.normcase(str(root))
    except ValueError:
        return False


def protected_roots(repo_root: Path) -> tuple[Path, ...]:
    return tuple(_resolved(repo_root / item) for item in PROTECTED_RELATIVE_ROOTS)


def guard_new_output(path: os.PathLike[str] | str, repo_root: Path, allowed_roots: tuple[Path, ...]) -> Path:
    """Resolve `path` or raise before any write happens."""
    candidate = _resolved(path)
    root = _resolved(repo_root)
    if any(_within(candidate, item) for item in protected_roots(root)):
        raise UnsafeWritePath(f"refusing to write inside a protected historical root: {candidate}")
    resolved_allowed = tuple(_resolved(item) for item in allowed_roots)
    if not resolved_allowed:
        raise UnsafeWritePath("run spec declares no allowed output roots")
    if not any(_within(candidate, item) for item in resolved_allowed):
        raise UnsafeWritePath(f"output is outside every run-spec root: {candidate}")
    return candidate


def guard_new_output_noclobber(path: os.PathLike[str] | str, repo_root: Path, allowed_roots: tuple[Path, ...]) -> Path:
    target = guard_new_output(path, repo_root, allowed_roots)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing artifact: {target}")
    return target
