"""Normalise CRLF checkouts back to the LF content the freeze manifest pins.

DEFECT THIS REPAIRS (declared, not hidden)
------------------------------------------
`V02_FREEZE_MANIFEST.json` records, for every artifact it protects, a raw
SHA256 **and** a byte count.  Those values were produced from the `-fkg-v02`
working tree, whose text files are LF.  A fresh Windows checkout under
`core.autocrlf=true` materialises the same Git blobs with CRLF, so the raw hash
and the byte count differ and `run_phase2_calibration.verify()` fails closed:

    RuntimeError: frozen artifact mismatch:
      .../project/pacer_fkg_v02/run_phase2_calibration.py

The fix is an ENVIRONMENT normalisation, not a content change: a file is
rewritten to LF only when

    lf_sha256(file) == recorded sha256   AND   lf_bytes(file) == recorded bytes

so the operation is self-proving.  Any file that fails both the raw and the
LF-normalised test is left untouched and reported as a hard error.

Git blob identity is unchanged, and under `core.autocrlf=true` `git status`
stays clean because Git normalises the working tree before comparing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

REPO_DEFAULT = Path(__file__).resolve().parents[3]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def collect_records(node: Any, found: list[dict[str, Any]]) -> None:
    if isinstance(node, dict):
        if {"path", "sha256", "bytes"} <= set(node):
            found.append(node)
        for value in node.values():
            collect_records(value, found)
    elif isinstance(node, list):
        for value in node:
            collect_records(value, found)


def normalise(repo_root: Path, *, apply: bool) -> dict[str, Any]:
    manifest_path = repo_root / "project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))

    records: list[dict[str, Any]] = []
    collect_records(manifest, records)
    seen: dict[str, dict[str, Any]] = {}
    for record in records:
        seen.setdefault(str(record["path"]), record)

    results = {"ok_raw": [], "ok_lf_already": [], "normalised": [], "unsatisfied": [], "missing": []}
    for rel, record in sorted(seen.items()):
        path = repo_root / rel
        expected_hash = str(record["sha256"]).lower()
        expected_bytes = int(record["bytes"])
        if not path.is_file():
            results["missing"].append(rel)
            continue
        raw = path.read_bytes()
        if sha256_bytes(raw) == expected_hash and len(raw) == expected_bytes:
            results["ok_raw"].append(rel)
            continue
        normalised = raw.replace(b"\r\n", b"\n")
        if sha256_bytes(normalised) == expected_hash and len(normalised) == expected_bytes:
            if apply:
                path.write_bytes(normalised)
                results["normalised"].append({"path": rel, "bytes_before": len(raw), "bytes_after": len(normalised),
                                              "lf_sha256": expected_hash})
            else:
                results["ok_lf_already"].append(rel)
            continue
        results["unsatisfied"].append({"path": rel, "disk_bytes": len(raw),
                                       "disk_sha256": sha256_bytes(raw),
                                       "expected_bytes": expected_bytes, "expected_sha256": expected_hash})

    results["n_records"] = len(seen)
    results["applied"] = bool(apply)
    results["all_satisfied"] = not results["unsatisfied"] and not results["missing"]
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_DEFAULT)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()
    report = normalise(args.repo_root.resolve(), apply=args.apply)
    text = json.dumps(report, indent=2)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report["all_satisfied"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
