#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
PROJECT="$ROOT/project"
OUTPUT="${2:-$PROJECT/results/drugclip_m4_benchmark_v01}"
MAMBA="${MAMBA_EXE:-$HOME/.local/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

if [[ ! -x "$MAMBA" ]]; then
  echo "micromamba is required; set MAMBA_EXE if it is installed elsewhere" >&2
  exit 2
fi
CHECKPOINT="$PROJECT/tools/DrugCLIP/artifacts/checkpoint_best.pt"
if [[ ! -f "$CHECKPOINT" ]]; then
  echo "Missing official DrugCLIP checkpoint: $CHECKPOINT" >&2
  exit 3
fi

declare -A PDB_SHA256=(
  [7TRQ]="bccf53f6c49bcfa91c15e2151c66d6b0db82181f1536d92b314b9af69917becc"
  [7TRP]="b5b1af3a2974360dfb4fbad1f5930d9258b171eb9665004b0b7cf119e5517e7b"
  [7TRS]="53442edb5d1ce98d4fcbf9ae52d96f2b387c5ace1b913b40b54906227416fcc1"
)
mkdir -p "$PROJECT/data/pdb"
for pdb_id in 7TRQ 7TRP 7TRS; do
  pdb_file="$PROJECT/data/pdb/${pdb_id}.pdb"
  if [[ ! -f "$pdb_file" ]]; then
    curl --fail --location --retry 3 \
      "https://files.rcsb.org/download/${pdb_id}.pdb" --output "$pdb_file"
  fi
  actual_sha="$(sha256sum "$pdb_file" | awk '{print $1}')"
  if [[ "$actual_sha" != "${PDB_SHA256[$pdb_id]}" ]]; then
    echo "Unexpected SHA256 for $pdb_file: $actual_sha" >&2
    exit 4
  fi
done

mkdir -p "$OUTPUT"
run_py() { "$MAMBA" run -n drugclip-cpu python "$@"; }

run_py "$PROJECT/scripts/build_drugclip_inputs.py" molecules \
  --csv "$PROJECT/data/benchmarks/m4_pam_v1/pam_vs_inactive.csv" \
  --output "$OUTPUT/molecules.lmdb" \
  --id-column canonical_molecule_id --smiles-column canonical_smiles \
  --conformers 1 --seed 20260922

for pdb_id in 7TRQ 7TRP 7TRS; do
  run_py "$PROJECT/scripts/build_drugclip_inputs.py" pocket \
    --pdb "$PROJECT/data/pdb/${pdb_id}.pdb" \
    --output "$OUTPUT/${pdb_id}_pocket.lmdb" \
    --name "${pdb_id}_M4_allosteric" --chain R
done

run_py "$PROJECT/scripts/extract_drugclip_embeddings_cpu.py" \
  --drugclip "$PROJECT/tools/DrugCLIP" \
  --unicore "$PROJECT/tools/Uni-Core" \
  --checkpoint "$CHECKPOINT" \
  --molecules "$OUTPUT/molecules.lmdb" \
  --pockets "$OUTPUT/7TRQ_pocket.lmdb" "$OUTPUT/7TRP_pocket.lmdb" "$OUTPUT/7TRS_pocket.lmdb" \
  --output "$OUTPUT/embeddings.npz" --batch-size 16

run_py "$PROJECT/scripts/evaluate_drugclip_m4_baseline.py" \
  --embeddings "$OUTPUT/embeddings.npz" \
  --benchmark "$PROJECT/data/benchmarks/m4_pam_v1/pam_vs_inactive.csv" \
  --output "$OUTPUT/baseline.json"

echo "DrugCLIP M4 baseline complete: $OUTPUT/baseline.json"
