# Recipe sweep for the 13-target LOSO on Science-2026 (1 seed each).
$ErrorActionPreference = 'Continue'
$o = 'D:\CLC\project\results\drugclip_science2026\sweep_v01'
New-Item -ItemType Directory -Force -Path $o | Out-Null
$py = 'D:\CLC\project\scripts\evaluate_drugclip_13target_loso.py'
$g = 'D:\CLC\project\results\gpcr_drugclip_screening_v01\science2026_ensemble_embeddings.npz'
$p = 'D:\CLC\project\results\drugclip_science2026\litpcba_external_v01\science2026_90.projection.pt'
$e = 'D:\CLC\project\results\drugclip_science2026\litpcba_external_v01\science2026_90.npz'
$gp = 'D:\CLC\project\data\benchmarks\gpcr_drugclip_screening_v01\screening_pairs.csv'
$ep = 'D:\CLC\project\data\benchmarks\litpcba_drugclip_external_v01\pairs.csv'
function W($m) { $m | Add-Content "$o\sweep.log" }
$variants = @(
  @{n='rank16';   a=@('--rank','16')},
  @{n='epochs80'; a=@('--epochs','80')},
  @{n='lr1e2';    a=@('--lr','0.01')},
  @{n='retr05';   a=@('--retrieval-weight','0.5')},
  @{n='temp10';   a=@('--temperature','0.1')}
)
foreach ($v in $variants) {
  python $py --gpcr-representations $g --gpcr-pairs $gp `
    --external-representations $e --external-pairs $ep `
    --projection $p --output "$o\$($v.n)_seed20260925.json" `
    --epochs 40 --seed 20260925 --preserve-weight 0 @($v.a) `
    > "$o\$($v.n).stdout.log" 2> "$o\$($v.n).stderr.log"
  W "$($v.n) exit=$LASTEXITCODE"
}
W 'SWEEP DONE'
