#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
PROJECT="$ROOT/project"
BENCH="$PROJECT/data/benchmarks/gpcr_drugclip_screening_v01"
OUT="$PROJECT/results/gpcr_drugclip_screening_v01"
MAMBA="${MAMBA_EXE:-$HOME/.local/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
run_py() { "$MAMBA" run -n drugclip-cpu python "$@"; }

mkdir -p "$BENCH" "$OUT"

if [[ ! -f "$BENCH/audit.json" ]]; then
  run_py "$PROJECT/scripts/build_gpcr_drugclip_screening_benchmark.py" \
    --libraries "$PROJECT/tools/gpcr-am-ensemble-docking/compound_libraries" \
    --output-dir "$BENCH" --decoys-per-active 10
fi

if [[ ! -f "$OUT/published_docking_baselines.json" ]]; then
  run_py "$PROJECT/scripts/evaluate_gpcr_published_docking_baselines.py" \
    --pairs "$BENCH/screening_pairs.csv" \
    --scores-root "$PROJECT/tools/gpcr-am-ensemble-docking/docking_scores" \
    --output "$OUT/published_docking_baselines.json"
fi

if [[ ! -f "$OUT/ensemble_pockets.lmdb" ]]; then
  run_py "$PROJECT/scripts/build_drugclip_gpcr_ensemble_pockets.py" \
    --ensembles "$PROJECT/tools/gpcr-am-ensemble-docking/ensembles" \
    --output "$OUT/ensemble_pockets.lmdb" --metadata "$OUT/ensemble_pockets.csv"
fi

if [[ ! -f "$OUT/molecules.lmdb" ]]; then
  run_py "$PROJECT/scripts/build_drugclip_inputs.py" molecules \
    --csv "$BENCH/molecules.csv" --id-column canonical_molecule_id \
    --smiles-column canonical_smiles --conformers 1 --seed 20260925 \
    --skip-failures --workers "${DRUGCLIP_CONFORMER_WORKERS:-8}" \
    --output "$OUT/molecules.lmdb"
fi

if [[ ! -f "$OUT/ensemble_embeddings.npz" ]]; then
  run_py "$PROJECT/scripts/extract_drugclip_embeddings_cpu.py" \
    --drugclip "$PROJECT/tools/DrugCLIP" --unicore "$PROJECT/tools/Uni-Core" \
    --checkpoint "$PROJECT/tools/DrugCLIP/artifacts/checkpoint_best.pt" \
    --molecules "$OUT/molecules.lmdb" --pockets "$OUT/ensemble_pockets.lmdb" \
    --output "$OUT/ensemble_embeddings.npz" --batch-size 16 --progress-every 100
fi

run_py "$PROJECT/scripts/finetune_drugclip_gpcr_screening.py" \
  --representations "$OUT/ensemble_embeddings.npz" \
  --projection "$OUT/ensemble_embeddings.projection.pt" \
  --pairs "$BENCH/screening_pairs.csv" --output "$OUT/single_conformation_result.json"

run_py "$PROJECT/scripts/finetune_drugclip_gpcr_ensemble_screening.py" \
  --representations "$OUT/ensemble_embeddings.npz" \
  --projection "$OUT/ensemble_embeddings.projection.pt" \
  --pairs "$BENCH/screening_pairs.csv" --pocket-metadata "$OUT/ensemble_pockets.csv" \
  --output "$OUT/ensemble_result.json"

echo "Complete: $OUT"
