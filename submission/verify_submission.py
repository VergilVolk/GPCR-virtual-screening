from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import pandas as pd

from src.frozen_replay import replay_frozen

ROOT = Path(__file__).resolve().parent
USED_WEIGHTS = {
    "seed20260925.projection.pt": "15ff4333c08761eefc5aaffd169c785e7f422c200eb187704e9f8ac8fce7e09a",
    "seed20260926.projection.pt": "00d365cdfe7e1711540990122bb067320c408cf891507bc6ca9a8c4579a1dc4b",
    "seed20260927.projection.pt": "b8aedfb507381c228ca09c02a250c9f98d101467b3ad64657a0f7db875fd0da6",
}


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    for name, expected in USED_WEIGHTS.items():
        if file_hash(ROOT / "models/used_module2" / name) != expected:
            raise ValueError(f"Model hash mismatch: {name}")

    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        result = replay_frozen(ROOT / "data/frozen", folder, folder / "results.csv")
        expected = pd.read_csv(ROOT / "results/results.csv")
        pd.testing.assert_frame_equal(result.fillna(""), expected.fillna(""),
                                      check_dtype=False, atol=1e-12)

    audit = {
        "status": "PASS",
        "verified_rows": {"module1_library": 2605, "module1_audit": 2592,
                          "module2_top200": 200, "module3_six_channel": 8,
                          "module4_four_context": 3, "integrated_results": 8},
        "used_weight_sha256": USED_WEIGHTS,
        "result_sha256": {
            path.name: file_hash(path)
            for path in sorted((ROOT / "results").glob("*.csv"))
        },
        "note": "Frozen-result replay and artifact identity check; docking and MD are not rerun.",
    }
    (ROOT / "logs/run_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print("PASS: frozen four-module replay, model hashes and final table")


if __name__ == "__main__":
    main()
