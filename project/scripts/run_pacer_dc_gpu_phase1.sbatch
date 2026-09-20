#!/bin/bash
#SBATCH --job-name=pacer-dc-p1
#SBATCH --array=0-17
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=project/results/pacer_dc_production_v01/slurm-%A_%a.out
#SBATCH --error=project/results/pacer_dc_production_v01/slurm-%A_%a.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

manifest="project/config/pacer_dc_gpu_phase1_manifest.csv"
row="$(awk -F, -v id="${SLURM_ARRAY_TASK_ID}" 'NR>1 && $1==id {print; exit}' "${manifest}")"
if [[ -z "${row}" ]]; then
  echo "No manifest row for array id ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
IFS=, read -r array_id system replica seed target_ns stage boundary <<<"${row}"

"${PACER_PYTHON:-python}" project/scripts/run_pacer_dc_production_md.py \
  "${system}" \
  --replica "${replica}" \
  --ns "${target_ns}" \
  --platform CUDA \
  --device-index "${PACER_DEVICE_INDEX:-0}" \
  --timestep-fs 2
