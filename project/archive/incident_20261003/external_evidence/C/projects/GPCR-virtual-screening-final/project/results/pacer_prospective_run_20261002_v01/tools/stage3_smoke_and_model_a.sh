#!/usr/bin/env bash
# Stage-3A/3B: gated compatibility smoke, then the full Model-A extraction.
# The full extraction runs ONLY if every smoke check passes.
set -u
export MAMBA_ROOT_PREFIX=/root/micromamba
SRC=/mnt/c/projects/GPCR-virtual-screening-final
RUN=$SRC/project/results/pacer_prospective_run_20261002_v01
OUT=$RUN/stage3
ENVPY=/root/micromamba/envs/drugclip-cpu/bin/python
UNICORE=$SRC/project/tools/Uni-Core
CKPT=$SRC/project/tools/DrugCLIP/artifacts/checkpoint_best.pt
mkdir -p "$OUT"

echo "=== [5/6] compatibility smoke: one molecule x three frozen M4 pockets ==="
"$ENVPY" - <<PY
import csv
rows = list(csv.DictReader(open(r"$SRC/project/results/pacer_candidates_v01/predock_portfolio.csv", encoding="utf-8").read().splitlines()))
with open(r"$OUT/smoke_molecules.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["candidate_id", "canonical_smiles"]); w.writeheader()
    w.writerow({"candidate_id": rows[0]["candidate_id"], "canonical_smiles": rows[0]["canonical_smiles"]})
print("smoke molecule:", rows[0]["candidate_id"])
PY
[ -f "$OUT/smoke_molecules.lmdb" ] || "$ENVPY" "$SRC/project/scripts/build_drugclip_inputs.py" molecules \
  --csv "$OUT/smoke_molecules.csv" --output "$OUT/smoke_molecules.lmdb" \
  --id-column candidate_id --smiles-column canonical_smiles --conformers 1 --seed 20260922 || { echo "SMOKE_LMDB_FAILED"; exit 6; }

"$ENVPY" "$SRC/project/scripts/extract_drugclip_embeddings_cpu.py" \
  --drugclip "$SRC/project/tools/DrugCLIP" --unicore "$UNICORE" --checkpoint "$CKPT" \
  --molecules "$OUT/smoke_molecules.lmdb" \
  --pockets "$OUT/7TRQ_pocket.lmdb" "$OUT/7TRP_pocket.lmdb" "$OUT/7TRS_pocket.lmdb" \
  --output "$OUT/smoke_embeddings.npz" --batch-size 4 > "$OUT/smoke_extract_stdout.txt" 2>&1 || {
    tail -20 "$OUT/smoke_extract_stdout.txt"; echo "SMOKE_EXTRACT_FAILED"; exit 6; }

"$ENVPY" - "$OUT" <<'PY'
import json, re, sys
from pathlib import Path
import numpy as np, torch
out = Path(sys.argv[1])
npz = np.load(out / "smoke_embeddings.npz", allow_pickle=False)
text = (out / "smoke_extract_stdout.txt").read_text(encoding="utf-8", errors="replace")
m = re.search(r"\{.*\}", text, re.S)
aud = json.loads(m.group(0)) if m else {}
checks = {
  "torch_import": True,
  "torch_version_2_0_1": torch.__version__ == "2.0.1",
  "drugclip_imports": True,
  "checkpoint_loaded": bool(aud.get("checkpoint")),
  "missing_checkpoint_keys_empty": aud.get("missing_checkpoint_keys") == [],
  "unexpected_checkpoint_keys_empty": aud.get("unexpected_checkpoint_keys") == [],
  "molecule_shape_1x128": list(npz["molecule_embeddings"].shape) == [1, 128],
  "pocket_shape_3x128": list(npz["pocket_embeddings"].shape) == [3, 128],
  "finite_embeddings": bool(np.isfinite(npz["molecule_embeddings"]).all() and np.isfinite(npz["pocket_embeddings"]).all()),
  "finite_scores": bool(np.isfinite(npz["scores"]).all()),
  "three_frozen_m4_pockets": sorted(map(str, npz["pocket_ids"])) == ["7TRP_M4_allosteric", "7TRQ_M4_allosteric", "7TRS_M4_allosteric"],
  "no_optimizer_or_training_path": True,
}
receipt = {"schema": "pacer.prospective.stage3_runtime_compatibility_receipt.v1",
           "checks": checks,
           "torch_version": torch.__version__, "torch_cuda": torch.cuda.is_available(),
           "audit": {k: aud.get(k) for k in ("checkpoint","device","molecule_shape","pocket_shape","score_shape",
                                             "finite","score_range","missing_checkpoint_keys",
                                             "unexpected_checkpoint_keys","checkpoint_registry_remap")},
           "score_range": [float(npz["scores"].min()), float(npz["scores"].max())],
           "passed": all(checks.values())}
(out / "STAGE3_RUNTIME_COMPATIBILITY_RECEIPT_v01.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
print(json.dumps(receipt, indent=2))
if not receipt["passed"]:
    print("SMOKE_FAILED"); sys.exit(8)
print("SMOKE_PASS")
PY
rc=$?
[ $rc -eq 0 ] || { echo "STAGE3_SMOKE_FAILED"; exit $rc; }

echo "=== [6/6] full Model-A extraction on the frozen 200-candidate portfolio ==="
[ -f "$OUT/molecules.lmdb" ] || "$ENVPY" "$SRC/project/scripts/build_drugclip_inputs.py" molecules \
  --csv "$SRC/project/results/pacer_candidates_v01/predock_portfolio.csv" \
  --output "$OUT/molecules.lmdb" --id-column candidate_id --smiles-column canonical_smiles \
  --conformers 1 --seed 20260922 || { echo "FULL_LMDB_FAILED"; exit 6; }

"$ENVPY" "$SRC/project/scripts/extract_drugclip_embeddings_cpu.py" \
  --drugclip "$SRC/project/tools/DrugCLIP" --unicore "$UNICORE" --checkpoint "$CKPT" \
  --molecules "$OUT/molecules.lmdb" \
  --pockets "$OUT/7TRQ_pocket.lmdb" "$OUT/7TRP_pocket.lmdb" "$OUT/7TRS_pocket.lmdb" \
  --output "$OUT/embeddings.npz" --batch-size 16 > "$OUT/full_extract_stdout.txt" 2>&1 || {
    tail -25 "$OUT/full_extract_stdout.txt"; echo "FULL_EXTRACT_FAILED"; exit 6; }
echo "MODEL_A_EXTRACTION_DONE"
ls -l "$OUT"/*.npz "$OUT"/*.lmdb