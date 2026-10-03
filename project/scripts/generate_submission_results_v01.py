#!/usr/bin/env python3
"""Generate the standardized competition results file (附件5 第四节).

Main result entry point: reads frozen, committed artifacts and emits
``results.csv`` (UTF-8) with the required fields:

    候选编号 candidate_id | 所属赛道 track | SMILES | 关键预测指标 |
    对应模型与版本 | 三维结构文件 | 备注

All numbers carry their frozen provenance; nothing is re-scored here.
Rows are the frozen Stage3D MD shortlist when a freeze manifest is
supplied; otherwise the full PACER-200 routed ranking is emitted with
``md_review=pending`` so the file structure is reproducible pre-MD.

Usage:
    python project/scripts/generate_submission_results_v01.py \
        [--stage3d-manifest <freeze.json>] --output results.csv
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "project" / "artifacts" / "drugclip2023_m4_loto_pacer200_v01"

MODELS = (
    "drugclip2023_gpcr_loto(seed20260925/26/27, projection adapters); "
    "drugclip2026_13t_family_aug(seed20260925/26/27, projection adapters); "
    "frozen M4-safe router v01 (0.5/0.5, M4->2023)"
)
CLAIM = "高优先级预测 PAM 候选（计算假设）；未做 M4+ACh 功能实验，确认 PAM 数量=0"


def load_ranking() -> list[dict]:
    with (ART / "pacer200_m4_safe_routed_ranking.csv").open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage3d-manifest", type=Path, default=None,
                    help="Stage3D corrected freeze manifest (adds MD-shortlist rows)")
    ap.add_argument("--output", type=Path, default=ROOT / "results.csv")
    args = ap.parse_args()

    ranking = load_ranking()
    by_id = {r["pair_id"]: r for r in ranking}

    shortlist = None
    if args.stage3d_manifest and args.stage3d_manifest.exists():
        man = json.loads(args.stage3d_manifest.read_text(encoding="utf-8"))
        shortlist = man.get("frozen_poses", [])

    def row_for(cid: str, *, md_status: str, pose: str = "", gate: str = "") -> dict:
        r = by_id[cid]
        return {
            "candidate_id": cid,
            "track": "赛道三：GPCR (CHRM4) PAM 虚拟筛选",
            "smiles": r.get("canonical_smiles", ""),
            "key_metrics": (
                f"router_final_rank={r.get('final_rank', r.get('rank', ''))}; "
                f"score_2023_gpcr_loto={r.get('score_2023_gpcr_loto', '')}; "
                f"score_2026_13t_famaug={r.get('score_2026_13t_famaug', '')}; "
                f"structural_gate={gate or 'n/a'}; md_four_context={md_status}"
            ),
            "model_and_version": MODELS,
            "structure_file": pose,
            "remarks": CLAIM,
        }

    rows: list[dict] = []
    if shortlist:
        for p in shortlist:
            cid = p["candidate_id"]
            rows.append(row_for(
                cid,
                md_status="selected; production MD in progress (frozen pose)",
                pose=p["frozen_pose"]["path"],
                gate="PASS (cluster>=0.8, coverage>=0.5, contacts>=8)",
            ))
    else:
        for r in ranking:
            rows.append(row_for(r["pair_id"], md_status="pending"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    meta = {
        "schema": "pacer.submission.results.v1",
        "rows": len(rows),
        "mode": "stage3d_shortlist" if shortlist else "full_pacer200_pre_md",
        "claim_boundary": CLAIM,
        "provenance": {
            "ranking": "project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_safe_routed_ranking.csv",
            "router_audit": "pacer200_m4_safe_routed_ranking.audit.json",
            "stage3d_manifest": str(args.stage3d_manifest) if shortlist else None,
        },
    }
    Path(str(args.output).replace(".csv", ".meta.json")).write_text(
        json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"wrote {len(rows)} rows -> {args.output} (mode={meta['mode']})")


if __name__ == "__main__":
    main()
