# Migration: 13-target LOSO from old-2023 checkpoint to Science-2026 checkpoint
# Usage:
#   migrate_13target_to_science2026.ps1 -Mode wait     (coordinator: wait for extractor, then spawn lanes)
#   migrate_13target_to_science2026.ps1 -Mode laneA    (preserve-weight 0.2, 3 seeds)
#   migrate_13target_to_science2026.ps1 -Mode laneB    (preserve-weight 0.0, 3 seeds)
#   migrate_13target_to_science2026.ps1 -Mode finalize (ECFP + 3-seed summaries)
param([ValidateSet('wait','laneA','laneB','finalize')][string]$Mode = 'wait')
$ErrorActionPreference = 'Continue'
$proj  = 'D:\CLC\project'
$out   = "$proj\results\drugclip_science2026\migration_13target_v01"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$gpcrNpz   = "$proj\results\gpcr_drugclip_screening_v01\science2026_ensemble_embeddings.npz"
$projPt    = "$proj\results\gpcr_drugclip_screening_v01\science2026_ensemble_embeddings.projection.pt"
if (-not (Test-Path $projPt)) { $projPt = "$proj\results\drugclip_science2026\litpcba_external_v01\science2026_90.projection.pt" }
$extNpz    = "$proj\results\drugclip_science2026\litpcba_external_v01\science2026_90.npz"
$gpcrPairs = "$proj\data\benchmarks\gpcr_drugclip_screening_v01\screening_pairs.csv"
$extPairs  = "$proj\data\benchmarks\litpcba_drugclip_external_v01\pairs.csv"
$script    = "$proj\scripts\migrate_13target_to_science2026.ps1"
function Write-Log([string]$m) { "$(Get-Date -Format 'MM-dd HH:mm') $m" | Add-Content "$out\runner.log" }

function Invoke-Loso([double]$preserve, [int]$seed) {
  $tag = "preserve$preserve`_seed$seed"
  python "$proj\scripts\evaluate_drugclip_13target_loso.py" `
    --gpcr-representations $gpcrNpz --gpcr-pairs $gpcrPairs `
    --external-representations $extNpz --external-pairs $extPairs `
    --projection $projPt --output "$out\loso_$tag.json" `
    --epochs 40 --seed $seed --preserve-weight $preserve `
    > "$out\loso_$tag.stdout.log" 2> "$out\loso_$tag.stderr.log"
  Write-Log "done loso_$tag exit=$LASTEXITCODE"
}

switch ($Mode) {
  'wait' {
    Write-Log 'waiting for extractor'
    while ($true) {
      if (Test-Path $gpcrNpz) {
        $s1 = (Get-Item $gpcrNpz).Length; Start-Sleep -Seconds 60
        if ((Get-Item $gpcrNpz).Length -eq $s1) { break }
      } else { Start-Sleep -Seconds 60 }
    }
    Write-Log "extractor done; projection=$projPt"
    Start-Process -FilePath powershell -ArgumentList '-NoProfile','-ExecutionPolicy','Bypass','-File',$script,'-Mode','laneA' -WindowStyle Hidden
    Start-Process -FilePath powershell -ArgumentList '-NoProfile','-ExecutionPolicy','Bypass','-File',$script,'-Mode','laneB' -WindowStyle Hidden
    while ((Get-ChildItem $out -Filter 'loso_preserve*_seed*.json').Count -lt 6) { Start-Sleep -Seconds 120 }
    Start-Process -FilePath powershell -ArgumentList '-NoProfile','-ExecutionPolicy','Bypass','-File',$script,'-Mode','finalize' -WindowStyle Hidden
    Write-Log 'lanes spawned; finalize armed'
  }
  'laneA' { foreach ($s in 20260925,20260926,20260927) { Invoke-Loso 0.2 $s } }
  'laneB' { foreach ($s in 20260925,20260926,20260927) { Invoke-Loso 0.0 $s } }
  'finalize' {
    python "$proj\scripts\evaluate_drugclip_13target_ecfp_loso.py" `
      --predictions "$out\loso_preserve0.2_seed20260925.predictions.csv" `
      --output "$out\ecfp_loso.json" `
      > "$out\ecfp.stdout.log" 2> "$out\ecfp.stderr.log"
    Write-Log "ecfp exit=$LASTEXITCODE"
    foreach ($p in '0.2','0') {
      python "$proj\scripts\summarize_drugclip_13target_multiseed.py" `
        --prediction "$out\loso_preserve${p}_seed20260925.predictions.csv" "$out\loso_preserve${p}_seed20260926.predictions.csv" "$out\loso_preserve${p}_seed20260927.predictions.csv" `
        --ecfp-predictions "$out\ecfp_loso.predictions.csv" `
        --output "$out\summary_preserve${p}_3seed.json" `
        > "$out\summary_p$p.stdout.log" 2> "$out\summary_p$p.stderr.log"
      Write-Log "summary preserve$p exit=$LASTEXITCODE"
    }
    Write-Log 'MIGRATION 13-TARGET COMPLETE'
  }
}
