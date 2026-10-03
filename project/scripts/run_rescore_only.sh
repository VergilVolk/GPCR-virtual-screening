#!/bin/bash
# 仅 rescore 阶段（对接已完成）
set -u
export OUTD=/mnt/d/CLC/project/results/pacer_rerun_v01
export GNINA=/mnt/d/CLC/project/tools/gnina/gnina_v1.1_linux
rescore() {
  name=$1; idx=$2
  pose="$OUTD/vina_poses/${name}_pose$(printf %04d $idx).pdbqt"
  log="$OUTD/gnina_tmp/${name}_rescore$(printf %04d $idx).log"
  out="$OUTD/gnina_tmp/${name}_rescore$(printf %04d $idx).txt"
  [ -f "$out" ] && [ -s "$out" ] && return 0
  [ -f "$pose" ] || return 0
  timeout 600 "$GNINA" -r "$OUTD/vina_poses/${name}_R_meeko.pdbqt" -l "$pose" \
    --score_only > "$log" 2>&1
  grep -E "CNNaffinity|minimizedAffinity" "$log" > "$out" 2>/dev/null
  [ -s "$out" ] || rm -f "$out"
}
export -f rescore
for name in 7TRQ 7TRP 7TRS; do
  export name
  seq 0 199 | xargs -P 6 -I@@ bash -c 'rescore "$name" @@'
  echo "${name} rescore complete: $(ls $OUTD/gnina_tmp/${name}_rescore*.txt 2>/dev/null | wc -l)"
done
echo "RESCORE_ALL_DONE"
