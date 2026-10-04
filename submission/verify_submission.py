from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import pandas as pd
from rdkit import Chem

from src.frozen_replay import replay_frozen
from src.screen_core import screen

ROOT = Path(__file__).resolve().parent
USED_WEIGHTS = {
    "drugclip2023_m4_loto_seed20260925.projection.pt": "a954c7c16fe857b5e3b02f1fbfdf7b2b1284d3559d86448b33736e16dd5d7398",
    "drugclip2023_m4_loto_seed20260926.projection.pt": "3cc92da9f066a36f9d65b47b7859e69ee7a456d48047a3e1bbd509283d50d4b9",
    "drugclip2023_m4_loto_seed20260927.projection.pt": "578ab8349fa8723f5609067e6511ebdfab8541dab8558e883bb7affeb7a8780c",
}
VISNET_SHA256 = "b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417"


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
    visnet = ROOT / "project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth"
    if file_hash(visnet) != VISNET_SHA256:
        raise ValueError("Module 4 ViSNet checkpoint hash mismatch")

    molecule_ids, scores, _ = screen(
        ROOT / "data/inference/module2_molecule_representations.npz",
        ROOT / "data/inference/module2_m4_pocket.npz",
        sorted((ROOT / "models/used_module2").glob("*.pt")),
    )
    library = pd.read_csv(ROOT / "data/inference/module2_library.csv")
    graph_to_id = {
        Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True): candidate_id
        for candidate_id, smiles in zip(library.molecule_id, library.smiles)
    }
    inferred = pd.DataFrame({"smiles": molecule_ids, "score_inferred": scores})
    inferred["gen_id"] = inferred.smiles.map(
        lambda value: graph_to_id[Chem.MolToSmiles(Chem.MolFromSmiles(value), canonical=True)])
    inferred = inferred.sort_values("score_inferred", ascending=False).reset_index(drop=True)
    inferred["rank_inferred"] = inferred.index + 1
    expected_top = pd.read_csv(ROOT / "results/02_drugclip_top200.csv")
    check = expected_top.merge(inferred[["gen_id", "score_inferred", "rank_inferred"]], on="gen_id")
    if len(check) != 200 or (check.final_rank != check.rank_inferred).any():
        raise ValueError("Module 2 forward inference did not reproduce the frozen top 200")
    if (check.score - check.score_inferred).abs().max() > 5e-7:
        raise ValueError("Module 2 forward scores differ from the frozen table")

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
        "module4_visnet_sha256": VISNET_SHA256,
        "result_sha256": {
            path.name: file_hash(path)
            for path in sorted((ROOT / "results").glob("*.csv"))
        },
        "module2_forward_max_abs_error": float((check.score - check.score_inferred).abs().max()),
        "note": "Module 2 is recomputed from representations and weights. Docking and MD are frozen-output replays.",
    }
    (ROOT / "logs/run_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print("PASS: frozen four-module replay, model hashes and final table")


if __name__ == "__main__":
    main()
