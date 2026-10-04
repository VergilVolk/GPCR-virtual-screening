from __future__ import annotations
import hashlib, json, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from predict import run_pipeline

ROOT = Path(__file__).resolve().parent
WEIGHTS = {
    "drugclip2023_m4_loto_seed20260925.projection.pt": "a954c7c16fe857b5e3b02f1fbfdf7b2b1284d3559d86448b33736e16dd5d7398",
    "drugclip2023_m4_loto_seed20260926.projection.pt": "3cc92da9f066a36f9d65b47b7859e69ee7a456d48047a3e1bbd509283d50d4b9",
    "drugclip2023_m4_loto_seed20260927.projection.pt": "578ab8349fa8723f5609067e6511ebdfab8541dab8558e883bb7affeb7a8780c",
}

def main() -> None:
    for name, expected in WEIGHTS.items():
        actual = hashlib.sha256((ROOT / "models" / name).read_bytes()).hexdigest()
        assert actual == expected, name
    rng = np.random.default_rng(20261004)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        ids = np.array(["TEST001", "TEST002", "TEST003"])
        pd.DataFrame({"candidate_id": ids, "smiles": [
            "CCOc1ccc(NC(=O)c2ccccc2)cc1", "COc1ccc(C(=O)Nc2nccs2)cc1", "CCN1CCC(c2ccccc2)CC1"
        ]}).to_csv(tmp / "library.csv", index=False)
        np.savez_compressed(tmp / "molecules.npz", molecule_ids=ids,
                            molecule_representations=rng.normal(size=(3, 512)).astype("float32"))
        np.savez_compressed(tmp / "pocket.npz",
                            pocket_representations=rng.normal(size=(1, 512)).astype("float32"))
        docking = {"candidate_id": ids}
        for i, channel in enumerate(["Glide_PDB", "Glide_BEmin", "Glide_BEavg", "Vina_PDB", "Vina_BEmin", "Vina_BEavg"]):
            docking[f"{channel}_raw"] = [-7-i*.1, -6-i*.1, -5-i*.1]
            docking[f"{channel}_rankpct"] = [.9, .6, .3]
        pd.DataFrame(docking).to_csv(tmp / "docking.csv", index=False)
        summary = {"source_status": "STAGE4_FULL_FROZEN_EVALUATION_COMPLETE",
                   "claim_boundary": "synthetic smoke test", "dataset": {"production_ns_per_trajectory": 1},
                   "candidates": {"TEST001": {"cluster": 0, "interpretation": ["synthetic smoke test"]}}}
        (tmp / "stage4.json").write_text(json.dumps(summary), encoding="utf-8")
        result = run_pipeline(tmp / "library.csv", tmp / "molecules.npz", tmp / "pocket.npz",
                              tmp / "docking.csv", tmp / "run", tmp / "run/results.csv",
                              tmp / "stage4.json")
        assert len(result) == 3 and result.six_channel_complete.all()
        assert result.four_context_status.eq("complete_unresolved").sum() == 1
    print("PASS: model hashes and four-module synthetic smoke test")

if __name__ == "__main__":
    main()
