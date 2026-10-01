#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
PROJECT="$ROOT/project"
MAMBA="${MAMBA_EXE:-$HOME/.local/bin/micromamba}"
DRUGCLIP="$PROJECT/tools/DrugCLIP"
UNICORE="$PROJECT/tools/Uni-Core"
DRUGCLIP_COMMIT="7a3a3fa33673f8668c811790f2e4681c98af44ef"
UNICORE_COMMIT="44f6386f4dcd7137fc1e5d5e768117d635d64a26"

mkdir -p "$PROJECT/tools"
if [[ ! -d "$DRUGCLIP/.git" ]]; then
  git clone https://github.com/bowen-gao/DrugClip.git "$DRUGCLIP"
  git -C "$DRUGCLIP" checkout "$DRUGCLIP_COMMIT"
fi
if [[ ! -d "$UNICORE/.git" ]]; then
  git clone https://github.com/dptech-corp/Uni-Core.git "$UNICORE"
  git -C "$UNICORE" checkout "$UNICORE_COMMIT"
fi
[[ "$(git -C "$DRUGCLIP" rev-parse HEAD)" == "$DRUGCLIP_COMMIT" ]] || {
  echo "DrugCLIP checkout is not the pinned commit" >&2; exit 3;
}
[[ "$(git -C "$UNICORE" rev-parse HEAD)" == "$UNICORE_COMMIT" ]] || {
  echo "Uni-Core checkout is not the pinned commit" >&2; exit 3;
}
if [[ ! -x "$MAMBA" ]]; then
  echo "micromamba is required; install it once or set MAMBA_EXE" >&2
  exit 2
fi
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$($MAMBA shell hook --shell bash)"
micromamba create -y -f "$PROJECT/environment_drugclip_cpu.yml" || \
  micromamba update -y -n drugclip-cpu -f "$PROJECT/environment_drugclip_cpu.yml"
micromamba activate drugclip-cpu

python "$PROJECT/scripts/prepare_drugclip_cpu_checkout.py" \
  --checkout "$DRUGCLIP" \
  --audit "$DRUGCLIP/artifacts/cpu_patch_audit.json"

pushd "$UNICORE" >/dev/null
python setup.py install --disable-cuda-ext
popd >/dev/null

python - <<'PY'
import rdkit, torch, unicore
print({"torch": torch.__version__, "rdkit": rdkit.__version__, "cuda": torch.cuda.is_available()})
PY
