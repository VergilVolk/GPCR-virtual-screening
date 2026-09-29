# CGDA-on-Science-2026 auto-fire chain (official path):
# wait for 15/15 target metrics -> 3-seed evaluate_drugclip_science2026_cgda_loso
# (preregistered flagship recipe) -> 3-seed aggregate with paired bootstrap.
$ErrorActionPreference = 'Continue'
$proj = 'D:\CLC\project'
$full = "$proj\results\drugclip_science2026\full_litpcba"
$out  = "$proj\results\drugclip_science2026\cgda_full15_v01"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$log = "$out\autofire.log"
function W([string]$m) { "$(Get-Date -Format 'MM-dd HH:mm') $m" | Add-Content $log }

W 'watching for 15/15 target metrics'
while ((Get-ChildItem $full -Filter '*.metrics.json' -ErrorAction SilentlyContinue).Count -lt 15) {
  Start-Sleep -Seconds 300
}
W '15/15 reached'

foreach ($s in 20260925,20260926,20260927) {
  python "$proj\scripts\evaluate_drugclip_science2026_cgda_loso.py" `
    --embedding-root $full `
    --output "$out\cgda_small_hardrank_seed$s.json" `
    --minimum-targets 15 `
    --rank 2 --experts 2 --epochs 5 `
    --hard-negative-fraction 0.5 --ranking-weight 0.2 --ranking-margin 0.1 `
    --seed $s `
    > "$out\cgda_seed$s.stdout.log" 2> "$out\cgda_seed$s.stderr.log"
  W "cgda seed$s exit=$LASTEXITCODE"
}
python "$proj\scripts\aggregate_cgda_science2026_3seed.py" `
  --runs "$out\cgda_small_hardrank_seed20260925.json" "$out\cgda_small_hardrank_seed20260926.json" "$out\cgda_small_hardrank_seed20260927.json" `
  --output "$out\cgda_small_hardrank_3seed_summary.json" `
  > "$out\summary.stdout.log" 2> "$out\summary.stderr.log"
W "aggregate exit=$LASTEXITCODE"
W 'CGDA SCIENCE-2026 CHAIN COMPLETE'
