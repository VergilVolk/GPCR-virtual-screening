"""CLI for the frozen DrugCLIP apply-only production path.

    # validate a Model A archive produced by the frozen CPU extractor
    python -m project.drugclip_freeze.apply validate --archive embeddings.npz

    # Model B projection-only re-scoring of the same frozen archive
    python -m project.drugclip_freeze.apply model-b --archive embeddings.npz --output out/modelb_scores.csv

    # frozen dual-model decision (OLD_TOP50 INTERSECT NEW_TOP25)
    python -m project.drugclip_freeze.apply decide --pairs pairs.csv --output out/dual_model.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from project.drugclip_freeze import dual_model, model_a, model_b

REPO_ROOT = Path(__file__).resolve().parents[2]


def _emit(payload: dict) -> None:
    print(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="drugclip_freeze")
    sub = parser.add_subparsers(dest="command", required=True)

    p_validate = sub.add_parser("validate", help="validate a frozen Model A archive")
    p_validate.add_argument("--archive", type=Path, required=True)
    p_validate.add_argument("--audit", type=Path, default=None)

    p_model_b = sub.add_parser("model-b", help="frozen projection-only apply")
    p_model_b.add_argument("--archive", type=Path, required=True)
    p_model_b.add_argument("--projections", type=Path, nargs="*", default=None)
    p_model_b.add_argument("--output", type=Path, required=True)
    p_model_b.add_argument("--repo-root", type=Path, default=REPO_ROOT)

    p_decide = sub.add_parser("decide", help="frozen dual-model decision")
    p_decide.add_argument("--pairs", type=Path, required=True,
                          help="CSV with candidate_id, old_rank, new_rank")
    p_decide.add_argument("--output", type=Path, required=True)

    p_verify = sub.add_parser("verify-shortlist", help="re-derive the frozen rule on the committed v02 shortlist")
    p_verify.add_argument("--shortlist", type=Path, required=True)

    args = parser.parse_args(argv)

    if args.command == "validate":
        _emit(model_a.validate_archive(args.archive, audit_path=args.audit))
        return 0

    if args.command == "model-b":
        paths = args.projections or model_b.default_projection_paths(args.repo_root)
        report = model_b.score(args.archive, paths)
        csv_path = model_b.write_csv(report, args.output)
        report_path = Path(args.output).with_suffix(".json")
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        _emit({"status": "MODEL_B_SCORED", "csv": str(csv_path), "json": str(report_path),
               "n_molecules": len(report["molecule_ids"]), "seeds": report["seeds"]})
        return 0

    if args.command == "decide":
        with Path(args.pairs).open(newline="", encoding="utf-8-sig") as handle:
            records = list(csv.DictReader(handle))
        decision = dual_model.decide(records)
        Path(args.output).write_text(json.dumps(decision, indent=2), encoding="utf-8")
        _emit({"status": "DUAL_MODEL_DECIDED", "output": str(args.output),
               "n_dual_top": decision["n_dual_top"], "rule": decision["rule"]})
        return 0

    _emit(dual_model.verify_committed_shortlist(args.shortlist))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
