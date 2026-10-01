# 20-target fine-tuning auto-chain: wait for the 7 subset-reps npz ->
# build extended set -> LOSO x3 seeds -> ensemble eval -> log.
$ErrorActionPreference = 'Continue'
$proj = 'D:\CLC\project'
$sub  = "$proj\results\drugclip_science2026\subset_reps"
$out  = "$proj\results\drugclip_science2026\extended20_finetune_v01"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$log = "$out\chain.log"
function W($m) { "$(Get-Date -Format 'MM-dd HH:mm') $m" | Add-Content $log }
$targets = 'ADRB2','ESR1_ago','ESR1_ant','OPRK1','IDH1','PPARG','TP53'

W 'waiting for 7 subset-reps npz'
while ((Get-ChildItem $sub -Filter '*.npz' -ErrorAction SilentlyContinue | Where-Object { $_.Name -notmatch 'references' }).Count -lt 7) {
  Start-Sleep -Seconds 180
}
W 'subset-reps complete'

python "$proj\scripts\build_extended20_finetune_set.py" `
  --pilot-npz "$proj\results\drugclip_science2026\litpcba_external_v01\science2026_90.npz" `
  --subset-root $sub --targets ($targets -join ',') `
  --output-npz "$out\extended16_external.npz" --output-pairs "$out\extended16_pairs.csv" `
  > "$out\build.stdout.log" 2> "$out\build.stderr.log"
W "build exit=$LASTEXITCODE"
if ($LASTEXITCODE -ne 0) { W 'ABORT'; exit 1 }

foreach ($s in 20260925,20260926,20260927) {
  python "$proj\scripts\evaluate_drugclip_13target_loso.py" `
    --gpcr-representations "$proj\results\gpcr_drugclip_screening_v01\science2026_ensemble_embeddings.npz" `
    --gpcr-pairs "$proj\data\benchmarks\gpcr_drugclip_screening_v01\screening_pairs.csv" `
    --external-representations "$out\extended16_external.npz" `
    --external-pairs "$out\extended16_pairs.csv" `
    --projection "$proj\results\drugclip_science2026\litpcba_external_v01\science2026_90.projection.pt" `
    --output "$out\loso20_seed$s.json" --epochs 40 --seed $s --preserve-weight 0 `
    > "$out\loso20_seed$s.stdout.log" 2> "$out\loso20_seed$s.stderr.log"
  W "loso20 seed$s exit=$LASTEXITCODE"
}
python "$proj\scripts\ensemble_eval_13target.py" `
  --predictions "$out\loso20_seed20260925.predictions.csv" "$out\loso20_seed20260926.predictions.csv" "$out\loso20_seed20260927.predictions.csv" `
  --output "$out\ensemble20_3seed.json" --bootstrap 2000 `
  > "$out\ensemble.stdout.log" 2> "$out\ensemble.stderr.log"
W "ensemble exit=$LASTEXITCODE"
W 'EXTENDED20 CHAIN COMPLETE'
