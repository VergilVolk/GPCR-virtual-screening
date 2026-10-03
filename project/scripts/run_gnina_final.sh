#!/bin/bash
# GNINA 7TRQ 全对接（6 并发，安全内存）+ 7TRP/7TRS score_only 重打分
set -u
export LIG=/mnt/d/CLC/project/results/pacer_rerun_v01/ligands
export OUTD=/mnt/d/CLC/project/results/pacer_rerun_v01
export GNINA=/mnt/d/CLC/project/tools/gnina/gnina_v1.1_linux
mkdir -p "$OUTD/gnina_tmp" "$OUTD/vina_poses"

dock() {
  idx=$1
  out="$OUTD/gnina_tmp/7TRQ_out$(printf %04d $idx).pdbqt"
  log="$OUTD/gnina_tmp/7TRQ_out$(printf %04d $idx).log"
  [ -f "$out" ] && return 0
  timeout 1200 "$GNINA" -r "$OUTD/vina_poses/7TRQ_R_meeko.pdbqt" -l "$LIG/lig$(printf %04d $idx).pdbqt" \
    --center_x 107.948 --center_y 85.737 --center_z 70.412 \
    --size_x 22 --size_y 22 --size_z 22 \
    --scoring vinardo --exhaustiveness 8 --num_modes 9 --seed 42 --cpu 1 \
    -o "$out" > "$log" 2>&1
}
rescore() {
  name=$1; idx=$2
  pose="$OUTD/vina_poses/${name}_pose$(printf %04d $idx).pdbqt"
  log="$OUTD/gnina_tmp/${name}_rescore$(printf %04d $idx).log"
  out="$OUTD/gnina_tmp/${name}_rescore$(printf %04d $idx).txt"
  [ -f "$out" ] && return 0
  [ -f "$pose" ] || return 0
  timeout 300 "$GNINA" -r "$OUTD/vina_poses/${name}_R_meeko.pdbqt" -l "$pose" \
    --score_only > "$log" 2>&1 && grep -E "CNNaffinity|minimizedAffinity|vinardo" "$log" > "$out" 2>/dev/null
}
export -f dock rescore

# 阶段1：7TRQ 全对接（6 并发）
seq 0 199 | xargs -P 6 -I@@ bash -c 'dock @@'
echo "7TRQ dock complete"

# 阶段2：三受体 score_only 重打分（6 并发，快）
for name in 7TRQ 7TRP 7TRS; do
  export name
  seq 0 199 | xargs -P 6 -I@@ bash -c 'rescore "$name" @@'
  echo "${name} rescore complete"
done
echo "GNINA_ALL_DONE"
