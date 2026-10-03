#!/usr/bin/env bash
# Stage-3A completion: Uni-Core install, environment record, gated compatibility smoke.
# The full 200-molecule Model-A extraction runs ONLY if every smoke check passes.
set -u
export MAMBA_ROOT_PREFIX=/root/micromamba
MM=/root/miniforge3/bin/micromamba
SRC=/mnt/c/projects/GPCR-virtual-screening-final
RUN=$SRC/project/results/pacer_prospective_run_20261002_v01
OUT=$RUN/stage3
ENVPY=$MAMBA_ROOT_PREFIX/envs/drugclip-cpu/bin/python
UNICORE=$SRC/project/tools/Uni-Core
CKPT=$SRC/project/tools/DrugCLIP/artifacts/checkpoint_best.pt
mkdir -p "$OUT"

echo "=== [1/6] Uni-Core install (frozen step verbatim) ==="
if ! "$ENVPY" -c "import unicore" >/dev/null 2>&1; then
  pushd "$UNICORE" >/dev/null
  "$ENVPY" setup.py install --disable-cuda-ext >/tmp/unicore_install.log 2>&1 || {
    popd >/dev/null; tail -15 /tmp/unicore_install.log; echo "UNICORE_INSTALL_FAILED"; exit 7; }
  popd >/dev/null
fi
"$ENVPY" -c "import unicore; print('unicore OK')"

echo "=== [2/6] environment record ==="
"$MM" list -n drugclip-cpu --explicit > "$OUT/ENV_EXPLICIT_lock.txt" 2>/dev/null
"$MM" env export -n drugclip-cpu > "$OUT/ENV_EXPORT_final.yml" 2>/dev/null
"$ENVPY" - <<'PY'
import json, platform, sys
info = {"python": sys.version.split()[0], "platform": platform.platform()}
for m in ("torch", "rdkit", "numpy", "unicore"):
    try:
        mod = __import__(m); info[m] = getattr(mod, "__version__", "ok")
    except Exception as e:
        info[m] = f"FAIL {type(e).__name__}"
import torch; info["torch_cuda"] = torch.cuda.is_available()
try:
    import torch.version as tv
except Exception:
    pass
try:
    from torch.utils import collect_env  # noqa
except Exception:
    pass
print(json.dumps(info, indent=2))
PY
"$ENVPY" -c "import mkl" 2>/dev/null || true
grep -Ei '^(mkl|intel-openmp|pytorch|python) ' "$OUT/ENV_EXPLICIT_lock.txt" | head -8 || true

echo "=== [3/6] frozen asset identity re-check ==="
"$ENVPY" - <<PY
import hashlib, os, subprocess
ck = r"$CKPT"
h = hashlib.sha256(open(ck, "rb").read()).hexdigest()
assert os.path.getsize(ck) == 1183713459, "checkpoint size"
assert h == "dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e", h
for repo, commit in (("DrugCLIP", "7a3a3fa33673f8668c811790f2e4681c98af44ef"),
                     ("Uni-Core", "44f6386f4dcd7137fc1e5d5e768117d635d64a26")):
    head = subprocess.check_output(["git", "-C", f"$SRC/project/tools/{repo}", "rev-parse", "HEAD"], text=True).strip()
    assert head == commit, (repo, head)
print("ASSET_IDENTITY_OK")
PY

echo "=== [4/6] build the three frozen M4 pocket LMDBs (chain R) ==="
for pid in 7TRQ 7TRP 7TRS; do
  [ -f "$OUT/${pid}_pocket.lmdb" ] || "$ENVPY" "$SRC/project/scripts/build_drugclip_inputs.py" pocket \
    --pdb "$SRC/project/data/pdb/$pid.pdb" --output "$OUT/${pid}_pocket.lmdb" \
    --name "${pid}_M4_allosteric" --chain R || { echo "POCKET_BUILD_FAILED $pid"; exit 5; }
done
ls -l "$OUT"/*_pocket.lmdb

echo "=== [5/6] compatibility smoke: one molecule x three frozen M4 pockets ==="
"$ENVPY" - <<PY
import csv
rows = list(csv.DictReader(open(r"$SRC/project/results/pacer_candidates_v01/predock_portfolio.csv", encoding="utf-8").read().splitlines()))
with open(r"$OUT/smoke_molecules.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["candidate_id", "canonical_smiles"]); w.writeheader()
    w.writerow({"candidate_id": rows[0]["candidate_id"], "canonical_smiles": rows[0]["canonical_smiles"]})
print("smoke molecule:", rows[0]["candidate_id"])
PY
"$ENVPY" "$SRC/project/scripts/build_drugclip_inputs.py" molecules \
  --csv "$OUT/smoke_molecules.csv" --output "$OUT/smoke_molecules.lmdb" \
  --id-column candidate_id --smiles-column canonical_smiles --conformers 1 --seed 20260922 || { echo "SMOKE_LMDB_FAILED"; exit 6; }

"$ENVPY" "$SRC/project/scripts/extract_drugclip_embeddings_cpu.py" \
  --drugclip "$SRC/project/tools/DrugCLIP" --unicore "$UNICORE" --checkpoint "$CKPT" \
  --molecules "$OUT/smoke_molecules.lmdb" \
  --pockets "$OUT/7TRQ_pocket.lmdb" "$OUT/7TRP_pocket.lmdb" "$OUT/7TRS_pocket.lmdb" \
  --output "$OUT/smoke_embeddings.npz" --batch-size 4 || { echo "SMOKE_EXTRACT_FAILED"; exit 6; }

"$ENVPY" - "$OUT" <<'PY'
import json, sys
from pathlib import Path
import numpy as np
out = Path(sys.argv[1])
npz = np.load(out / "smoke_embeddings.npz", allow_pickle=False)
aud = json.loads((out / "smoke_embeddings.json").read_text(encoding="utf-8"))
checks = {
  "torch_2_0_1": __import__("torch").__version__ == "2.0.1",
  "drugclip_imports": True,
  "checkpoint_loaded": bool(aud.get("checkpoint")),
  "missing_checkpoint_keys_empty": aud.get("missing_checkpoint_keys") == [],
  "unexpected_checkpoint_keys_empty": aud.get("unexpected_checkpoint_keys") == [],
  "molecule_shape_1x128": list(npz["molecule_embeddings"].shape) == [1, 128],
  "pocket_shape_3x128": list(npz["pocket_embeddings"].shape) == [3, 128],
  "finite_embeddings": bool(np.isfinite(npz["molecule_embeddings"]).all() and np.isfinite(npz["pocket_embeddings"]).all()),
  "finite_scores": bool(np.isfinite(npz["scores"]).all()),
  "three_frozen_m4_pockets": sorted(map(str, npz["pocket_ids"])) == ["7TRP_M4_allosteric", "7TRQ_M4_allosteric", "7TRS_M4_allosteric"],
}
print(json.dumps({"smoke": checks, "audit": aud, "score_range": [float(npz["scores"].min()), float(npz["scores"].max())]}, indent=2))
if not all(checks.values()):
    print("SMOKE_FAILED"); sys.exit(8)
print("SMOKE_PASS")
PY
rc=$?
[ $rc -eq 0 ] || { echo "STAGE3_SMOKE_FAILED"; exit $rc; }

echo "=== [6/6] full Model-A extraction on the frozen 200-candidate portfolio ==="
"$ENVPY" "$SRC/project/scripts/build_drugclip_inputs.py" molecules \
  --csv "$SRC/project/results/pacer_candidates_v01/predock_portfolio.csv" \
  --output "$OUT/molecules.lmdb" --id-column candidate_id --smiles-column canonical_smiles \
  --conformers 1 --seed 20260922 || { echo "FULL_LMDB_FAILED"; exit 6; }

"$ENVPY" "$SRC/project/scripts/extract_drugclip_embeddings_cpu.py" \
  --drugclip "$SRC/project/tools/DrugCLIP" --unicore "$UNICORE" --checkpoint "$CKPT" \
  --molecules "$OUT/molecules.lmdb" \
  --pockets "$OUT/7TRQ_pocket.lmdb" "$OUT/7TRP_pocket.lmdb" "$OUT/7TRS_pocket.lmdb" \
  --output "$OUT/embeddings.npz" --batch-size 16 || { echo "FULL_EXTRACT_FAILED"; exit 6; }
echo "MODEL_A_EXTRACTION_DONE"
ls -l "$OUT"