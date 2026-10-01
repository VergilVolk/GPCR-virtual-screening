"""Run the CM00734 payload-identity regression and write its receipt.

    python project/integration_freeze/tools/run_regression.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.pacer_fkg_v02.frozen_runner import blocks, regression  # noqa: E402


def main() -> int:
    blocks_proof = blocks.prove_equivalence()
    report = regression.run_regression(REPO)
    report["block_constructor_equivalence"] = blocks_proof

    out = REPO / "project/integration_freeze/receipts/CM00734_PAYLOAD_REGRESSION_v01.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    summary = {k: v for k, v in report.items() if k != "differences"}
    print(json.dumps(summary, indent=2))
    print("\nscientific_payload_identical =", report["scientific_payload_identical"])
    print("block_equivalence_all_identical =", blocks_proof.get("all_identical"))
    return 0 if report["scientific_payload_identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
