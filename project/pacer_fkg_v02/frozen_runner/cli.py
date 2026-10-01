"""Single entry point for the frozen runner.

    python -m project.pacer_fkg_v02.frozen_runner.cli --spec <spec.json> --phase {1,2a,2b,all} --run
    python -m project.pacer_fkg_v02.frozen_runner.cli --spec <spec.json> --phase {1,2a,2b,all} --verify

Every phase re-verifies the frozen anchor before doing any work.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from project.pacer_fkg_v02.frozen_runner import engines
from project.pacer_fkg_v02.frozen_runner.spec import RunSpec, assert_frozen_spec, authenticate_inputs, load_run_spec

REPO_ROOT = Path(__file__).resolve().parents[3]


def _emit(payload: dict) -> None:
    print(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="frozen_runner")
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--phase", choices=("1", "2a", "2b", "all"), required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--checkpoint", type=Path, default=None)
    parser.add_argument("--geom2vec-source", type=Path, default=None)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args(argv)

    repo_root = args.repo_root.resolve()
    spec: RunSpec = load_run_spec(args.spec, repo_root)
    assert_frozen_spec(spec)

    preflight = {
        "run_id": spec.run_id,
        "candidate_id": spec.molecule.candidate_id,
        "frozen_numerics": engines.frozen_numerics_report()["ok"],
        "frozen_module_ledger": engines.frozen_module_ledger(repo_root),
        "input_authentication": authenticate_inputs(spec, repo_root),
    }

    from project.pacer_fkg_v02.frozen_runner import fkg_phase1, fkg_phase2a, fkg_phase2b

    results: dict[str, object] = {"preflight": preflight}
    phases = ["1", "2a", "2b"] if args.phase == "all" else [args.phase]
    for phase in phases:
        if phase == "1":
            if args.run:
                results["phase1"] = fkg_phase1.run(spec, repo_root=repo_root, device=args.device,
                                                   checkpoint=args.checkpoint,
                                                   geom2vec_source=args.geom2vec_source)
            else:
                results["phase1"] = fkg_phase1.verify(spec, repo_root=repo_root)
        elif phase == "2a":
            results["phase2a"] = (fkg_phase2a.run(spec, repo_root=repo_root) if args.run
                                  else fkg_phase2a.verify(spec, repo_root=repo_root))
        else:
            results["phase2b"] = (fkg_phase2b.run(spec, repo_root=repo_root) if args.run
                                  else fkg_phase2b.verify(spec, repo_root=repo_root))
    _emit(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
