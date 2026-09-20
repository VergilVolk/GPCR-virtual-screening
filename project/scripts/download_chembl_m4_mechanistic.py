"""Download ChEMBL M4 records carrying explicit allostery-model endpoints."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import requests


P = Path(__file__).resolve().parents[1]
OUT = P / "data" / "benchmarks" / "m4_pam_v1" / "mechanistic_chembl"
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TARGET = "CHEMBL1821"
TYPES = ("pKb", "log(activity)", "log(IA)")


def fetch(activity_type: str) -> list[dict]:
    offset, records = 0, []
    while True:
        response = requests.get(
            BASE,
            params={"target_chembl_id": TARGET, "standard_type": activity_type,
                    "limit": 1000, "offset": offset},
            headers={"User-Agent": "PACER-M4 academic benchmark"}, timeout=90,
        )
        response.raise_for_status()
        payload = response.json()
        records.extend(payload["activities"])
        next_page = payload["page_meta"].get("next")
        if not next_page:
            return records
        offset += int(payload["page_meta"]["limit"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_records, counts = [], {}
    for activity_type in TYPES:
        records = fetch(activity_type)
        counts[activity_type] = len(records)
        all_records.extend(records)
    # A ChEMBL activity may occur in multiple filtered requests only if its
    # standard type is malformed; make the uniqueness assertion explicit.
    by_id = {int(item["activity_id"]): item for item in all_records}
    if len(by_id) != len(all_records):
        raise ValueError("Duplicate activity IDs across endpoint requests")
    payload = {
        "target_chembl_id": TARGET,
        "query_standard_types": list(TYPES),
        "counts": counts,
        "activities": list(by_id.values()),
    }
    output = OUT / "raw_mechanistic_activities.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    audit = {
        "source": "ChEMBL REST API",
        "target": TARGET,
        "n_records": len(by_id),
        "counts": counts,
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }
    (OUT / "download_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
