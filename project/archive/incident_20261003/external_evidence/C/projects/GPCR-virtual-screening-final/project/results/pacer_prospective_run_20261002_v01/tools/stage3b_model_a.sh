#!/usr/bin/env bash
# Stage-3B: frozen DrugCLIP Model A apply-only path (official checkpoint).
set -euo pipefail
SRC=/mnt/c/projects/GPCR-virtual-screening-final
PY=/root/micromamba/envs/drugclip-cpu/bin/python
OUT=$SRC/project/results/pacer_prospective_run_20261002_v01/stage3
CKPT=$SRC/project/tools/DrugCLIP/artifacts/checkpoint_best.pt
mkdir -p "$OUT"

echo "=== asset identity gate ==="
"$PY" - <<PY
import hashlib, os, subprocess, sys
ck = r"$CKPT"
h = hashlib.sha256(open(ck, "rb").read()).hexdigest()
size = os.path.getsize(ck)
exp_h = "dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e"
assert size == 1183713459, f"checkpoint size {size}"
assert h == exp_h, f"checkpoint sha256 {h}"
for repo, commit in (("DrugCLIP", "7a3a3fa33673f8668c811790f2e4681c98af44ef"),
                     ("Uni-Core", "44f6386f4dcd7137fc1e5d5e768117d635d64a26")):
    head = subprocess.check_output(["git", "-C", f"$SRC/project/tools/{repo}", "rev-parse", "HEAD"], text=True).strip()
    assert head == commit, f"{repo} HEAD {head}"
print("ASSET_IDENTITY_OK", size, h[:16])
PY

echo "=== molecules lmdb from the Stage-2 portfolio ==="
"$PY" "$SRC/project/scripts/build_drugclip_inputs.py" molecules \
  --csv "$SRC/project/results/pacer_candidates_v01/predock_portfolio.csv" \
  --output "$OUT/molecules.lmdb" --id-column candidate_id --smiles-column canonical_smiles \
  --conformers 1 --seed 20260922

echo "=== pocket lmdb x3 (chain R) ==="
for pid in 7TRQ 7TRP 7TRS; do
  "$PY" "$SRC/project/scripts/build_drugclip_inputs.py" pocket \
    --pdb "$SRC/project/data/pdb/$pid.pdb" --output "$OUT/${pid}_pocket.lmdb" \
    --name "${pid}_M4_allosteric" --chain R
done

echo "=== frozen Model A apply-only extraction ==="
"$PY" "$SRC/project/scripts/extract_drugclip_embeddings_cpu.py" \
  --drugclip "$SRC/project/tools/DrugCLIP" --unicore "$SRC/project/tools/Uni-Core" \
  --checkpoint "$CKPT" --molecules "$OUT/molecules.lmdb" \
  --pockets "$OUT/7TRQ_pocket.lmdb" "$OUT/7TRP_pocket.lmdb" "$OUT/7TRS_pocket.lmdb" \
  --output "$OUT/embeddings.npz" --batch-size 16
echo "MODEL_A_EXTRACTION_DONE"
