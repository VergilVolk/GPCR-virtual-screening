#!/usr/bin/env bash
set -euo pipefail

# Full PACER binding-mainline run on the verified Science-2026 DrugCLIP weight.
# Usage: bash project/scripts/run_science2026_pacer_mainline_gpu.sh
# Override paths or batch size with environment variables when needed.

ROOT="${PACER_ROOT:-$(pwd)}"
DRUGCLIP="${DRUGCLIP_ROOT:-$ROOT/project/tools/DrugCLIP}"
UNICORE="${UNICORE_ROOT:-$ROOT/project/tools/Uni-Core}"
CHECKPOINT="${DRUGCLIP_CHECKPOINT:-$ROOT/project/data/external/drugclip_science2026/litpcba_identity_90.pt}"
DATA_ROOT="${LITPCBA_ROOT:-$ROOT/project/data/external/drugclip_science2026/LIT-PCBA/lit_pcba}"
SOURCE_ZIP="${LITPCBA_ZIP:-$ROOT/project/data/external/drugclip_science2026/LIT-PCBA.zip}"
OUTPUT="${SCIENCE2026_OUTPUT:-$ROOT/project/results/drugclip_science2026/full_litpcba}"
BATCH_SIZE="${DRUGCLIP_BATCH_SIZE:-64}"

python - <<'PY'
import torch
if not torch.cuda.is_available():
    raise SystemExit("CUDA gate failed: torch.cuda.is_available() is false")
print("cuda_device=", torch.cuda.get_device_name(0))
PY

python "$ROOT/project/scripts/extract_drugclip_science2026_references.py" \
  --drugclip "$DRUGCLIP" --unicore "$UNICORE" --checkpoint "$CHECKPOINT" \
  --data-root "$DATA_ROOT" --output "$OUTPUT" --device cuda \
  --batch-size "$BATCH_SIZE" --trusted-checkpoint

# Target-level outputs are resumable: completed NPZ+metrics pairs are skipped.
python "$ROOT/project/scripts/run_drugclip_science2026_full_litpcba.py" \
  --drugclip "$DRUGCLIP" --unicore "$UNICORE" --checkpoint "$CHECKPOINT" \
  --data-root "$DATA_ROOT" --source-zip "$SOURCE_ZIP" --output "$OUTPUT" \
  --device cuda --batch-size "$BATCH_SIZE" --trusted-checkpoint

# Three independently initialized strict target-LOO CGDA runs.
for seed in 20260925 20260926 20260927; do
  python "$ROOT/project/scripts/evaluate_drugclip_science2026_cgda_loso.py" \
    --embedding-root "$OUTPUT" \
    --output "$ROOT/project/results/drugclip_science2026/full15_cgda_seed${seed}.json" \
    --seed "$seed"
done

echo "PACER Science-2026 full-data mainline complete"
