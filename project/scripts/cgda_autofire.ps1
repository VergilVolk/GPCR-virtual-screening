# CGDA-on-Science-2026 auto-fire chain:
# wait for 15/15 target metrics -> build pickle root -> 3-seed CGDA (frozen
# small_hardrank recipe) -> 3-seed summary -> log. Detached-safe.
$ErrorActionPreference = 'Continue'
$proj = 'D:\CLC\project'
$full = "$proj\results\drugclip_science2026\full_litpcba"
$out  = "$proj\results\drugclip_science2026\cgda_full15_v01"
$root = "$proj\data\external\drugclip_science2026\cgda_root"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$log = "$out\autofire.log"
function W([string]$m) { "$(Get-Date -Format 'MM-dd HH:mm') $m" | Add-Content $log }

W 'watching for 15/15 target metrics'
while ((Get-ChildItem $full -Filter '*.metrics.json' -ErrorAction SilentlyContinue).Count -lt 15) {
  Start-Sleep -Seconds 300
}
W '15/15 reached; building CGDA pickle root'
python "$proj\scripts\build_cgda_root_science2026.py" `
  --full-litpcba $full --output-root $root --pockets-pkl "$root\pockets_map.pkl" `
  > "$out\build_root.stdout.log" 2> "$out\build_root.stderr.log"
W "build_root exit=$LASTEXITCODE"
if ($LASTEXITCODE -ne 0) { W 'ABORT: root build failed'; exit 1 }

foreach ($s in 20260925,20260926,20260927) {
  python "$proj\scripts\run_cgda_science2026.py" `
    --root $root --pocket-root $root --pocket-archive "$root\pockets_map.pkl" `
    --output "$out\cgda_small_hardrank_seed$s.json" `
    --rank 2 --experts 2 --epochs 5 `
    --hard-negative-fraction 0.5 --ranking-weight 0.2 --ranking-margin 0.1 `
    --seed $s `
    > "$out\cgda_seed$s.stdout.log" 2> "$out\cgda_seed$s.stderr.log"
  W "cgda seed$s exit=$LASTEXITCODE"
}
python "$proj\scripts\summarize_drugclip_cgda_multiseed.py" `
  --runs "$out\cgda_small_hardrank_seed20260925.json" "$out\cgda_small_hardrank_seed20260926.json" "$out\cgda_small_hardrank_seed20260927.json" `
  --output "$out\cgda_small_hardrank_3seed_summary.json" `
  > "$out\summary.stdout.log" 2> "$out\summary.stderr.log"
W "summary exit=$LASTEXITCODE"
W 'CGDA SCIENCE-2026 CHAIN COMPLETE'
