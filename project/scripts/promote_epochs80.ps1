# Promote sweep winner (epochs80) to 3 seeds + one combo probe (ep80+lr1e2).
$ErrorActionPreference = 'Continue'
$o = 'D:\CLC\project\results\drugclip_science2026\sweep_v01'
$py = 'D:\CLC\project\scripts\evaluate_drugclip_13target_loso.py'
$g = 'D:\CLC\project\results\gpcr_drugclip_screening_v01\science2026_ensemble_embeddings.npz'
$p = 'D:\CLC\project\results\drugclip_science2026\litpcba_external_v01\science2026_90.projection.pt'
$e = 'D:\CLC\project\results\drugclip_science2026\litpcba_external_v01\science2026_90.npz'
$gp = 'D:\CLC\project\data\benchmarks\gpcr_drugclip_screening_v01\screening_pairs.csv'
$ep = 'D:\CLC\project\data\benchmarks\litpcba_drugclip_external_v01\pairs.csv'
function W($m) { $m | Add-Content "$o\promote.log" }
foreach ($s in 20260926,20260927) {
  python $py --gpcr-representations $g --gpcr-pairs $gp `
    --external-representations $e --external-pairs $ep `
    --projection $p --output "$o\epochs80_seed$s.json" `
    --epochs 80 --seed $s --preserve-weight 0 `
    > "$o\epochs80_seed$s.stdout.log" 2> "$o\epochs80_seed$s.stderr.log"
  W "epochs80 seed$s exit=$LASTEXITCODE"
}
python $py --gpcr-representations $g --gpcr-pairs $gp `
  --external-representations $e --external-pairs $ep `
  --projection $p --output "$o\ep80lr1e2_seed20260925.json" `
  --epochs 80 --lr 0.01 --seed 20260925 --preserve-weight 0 `
  > "$o\ep80lr1e2.stdout.log" 2> "$o\ep80lr1e2.stderr.log"
W "ep80lr1e2 exit=$LASTEXITCODE"
# 3-seed ensemble for epochs80
python D:\CLC\project\scripts\ensemble_eval_13target.py `
  --predictions "$o\epochs80_seed20260925.predictions.csv" "$o\epochs80_seed20260926.predictions.csv" "$o\epochs80_seed20260927.predictions.csv" `
  --output "$o\epochs80_ensemble.json" --bootstrap 2000 `
  > "$o\epochs80_ens.stdout.log" 2> "$o\epochs80_ens.stderr.log"
W "ensemble exit=$LASTEXITCODE"
W 'PROMOTE DONE'
