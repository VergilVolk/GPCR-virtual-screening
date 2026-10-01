# 20-target family-aug LOSO x3 seeds (detached, resumable per-seed via output existence check)
$ErrorActionPreference = 'Continue'
$o = 'D:\CLC\project\results\drugclip_science2026\family_aug_v01'
$fam = 'D:\CLC\project\data\drugclip_muscarinic_family_aug_v01'
$log = "$o\famaug20.log"
function W($m) { "$(Get-Date -Format 'HH:mm') $m" | Add-Content $log }
foreach ($s in 20260925, 20260926, 20260927) {
  if (Test-Path "$o\famaug20_ep80_seed$s.json") { W "skip seed$s"; continue }
  python D:\CLC\project\scripts\evaluate_drugclip_loso_family_aug.py `
    --gpcr-representations D:\CLC\project\results\gpcr_drugclip_screening_v01\science2026_ensemble_embeddings.npz `
    --gpcr-pairs D:\CLC\project\data\benchmarks\gpcr_drugclip_screening_v01\screening_pairs.csv `
    --external-representations D:\CLC\project\results\drugclip_science2026\extended20_finetune_v01\extended16_external.npz `
    --external-pairs D:\CLC\project\results\drugclip_science2026\extended20_finetune_v01\extended16_pairs.csv `
    --projection D:\CLC\project\results\drugclip_science2026\litpcba_external_v01\science2026_90.projection.pt `
    --family-npz "$fam\family_science2026.npz" --family-pairs "$fam\family_pairs_clean.csv" `
    --output "$o\famaug20_ep80_seed$s.json" --seed $s `
    > "$o\famaug20_seed$s.out.log" 2> "$o\famaug20_seed$s.err.log"
  W "seed$s exit=$LASTEXITCODE"
}
python D:\CLC\project\scripts\ensemble_eval_13target.py `
  --predictions "$o\famaug20_ep80_seed20260925.predictions.csv" "$o\famaug20_ep80_seed20260926.predictions.csv" "$o\famaug20_ep80_seed20260927.predictions.csv" `
  --output "$o\famaug20_ensemble.json" --bootstrap 2000 `
  > "$o\famaug20_ens.out.log" 2> "$o\famaug20_ens.err.log"
W "ensemble exit=$LASTEXITCODE"
W 'FAMAUG20 DONE'
