# Subset interim gate: when the 7 subset npz files exist ->
# references -> merge 8 full + 7 subset -> NEW-side CGDA x3 ->
# OLD-side subset root -> OLD-side CGDA x3 -> aggregates -> done.
$ErrorActionPreference = 'Continue'
$proj = 'D:\CLC\project'
$sub  = "$proj\results\drugclip_science2026\subset_litpcba"
$full = "$proj\results\drugclip_science2026\full_litpcba"
$merg = "$proj\results\drugclip_science2026\subset_merged15"
$oldroot = "$proj\data\external\drugclip_official\selected\lit_pcba"
$oldsb   = "$proj\data\external\drugclip_science2026\cgda_old_subset_root"
$out = "$proj\results\drugclip_science2026\cgda_subset_gate_v01"
New-Item -ItemType Directory -Force -Path $out, $merg | Out-Null
$log = "$out\chain.log"
function W([string]$m) { "$(Get-Date -Format 'MM-dd HH:mm') $m" | Add-Content $log }
$subsetTargets = 'OPRK1','VDR','FEN1','GBA','IDH1','KAT2A','ADRB2'

W 'waiting for 7 subset npz'
while ((Get-ChildItem $sub -Filter '*.npz' -ErrorAction SilentlyContinue | Where-Object { $_.Name -notmatch 'references' }).Count -lt 7) {
  Start-Sleep -Seconds 120
}
W 'subset npz complete'

# references for all 15 into merged dir (fast, ligands only)
python "$proj\scripts\extract_drugclip_science2026_references.py" `
  --drugclip "$proj\tools\DrugCLIP" --unicore "$proj\tools\Uni-Core" `
  --checkpoint "$proj\data\external\drugclip_science2026\litpcba_identity_90.pt" `
  --data-root "$proj\data\external\drugclip_science2026\LIT-PCBA\lit_pcba" `
  --output $merg --trusted-checkpoint `
  > "$out\references.stdout.log" 2> "$out\references.stderr.log"
W "references exit=$LASTEXITCODE"

# merge: 8 full npz + 7 subset npz
foreach ($f in (Get-ChildItem $full -Filter '*.npz')) { Copy-Item $f.FullName $merg -Force }
foreach ($t in $subsetTargets) { Copy-Item "$sub\$t.npz" $merg -Force }
$n = (Get-ChildItem $merg -Filter '*.npz' | Where-Object Name -notmatch 'references').Count
W "merged targets = $n"
if ($n -ne 15) { W 'ABORT: merged != 15'; exit 1 }

# NEW side: CGDA x3
foreach ($s in 20260925,20260926,20260927) {
  python "$proj\scripts\evaluate_drugclip_science2026_cgda_loso.py" `
    --embedding-root $merg --output "$out\new_cgda_seed$s.json" --minimum-targets 15 `
    --rank 2 --experts 2 --epochs 5 --hard-negative-fraction 0.5 `
    --ranking-weight 0.2 --ranking-margin 0.1 --seed $s `
    > "$out\new_seed$s.stdout.log" 2> "$out\new_seed$s.stderr.log"
  W "new cgda seed$s exit=$LASTEXITCODE"
}
python "$proj\scripts\aggregate_cgda_science2026_3seed.py" `
  --runs "$out\new_cgda_seed20260925.json" "$out\new_cgda_seed20260926.json" "$out\new_cgda_seed20260927.json" `
  --output "$out\new_cgda_3seed_summary.json" `
  > "$out\new_summary.stdout.log" 2> "$out\new_summary.stderr.log"
W "new aggregate exit=$LASTEXITCODE"

# OLD side: subset pickle root + CGDA x3 (original runner, LMDB pockets)
python "$proj\scripts\build_cgda_root_old_subset.py" `
  --old-root $oldroot --subset-npz-root $sub --output-root $oldsb `
  > "$out\oldroot.stdout.log" 2> "$out\oldroot.stderr.log"
W "old root build exit=$LASTEXITCODE"
foreach ($s in 20260925,20260926,20260927) {
  python "$proj\scripts\evaluate_drugclip_cgda_loso.py" `
    --root $oldsb `
    --pocket-root "$proj\data\external\drugclip_official\selected\pockets" `
    --pocket-archive "$proj\results\litpcba_drugclip_external_v01\official_full15_pockets.npz" `
    --output "$out\old_cgda_seed$s.json" `
    --rank 2 --experts 2 --epochs 5 --hard-negative-fraction 0.5 `
    --ranking-weight 0.2 --ranking-margin 0.1 --seed $s `
    > "$out\old_seed$s.stdout.log" 2> "$out\old_seed$s.stderr.log"
  W "old cgda seed$s exit=$LASTEXITCODE"
}
python "$proj\scripts\aggregate_cgda_science2026_3seed.py" `
  --runs "$out\old_cgda_seed20260925.json" "$out\old_cgda_seed20260926.json" "$out\old_cgda_seed20260927.json" `
  --output "$out\old_cgda_3seed_summary.json" `
  > "$out\old_summary.stdout.log" 2> "$out\old_summary.stderr.log"
W "old aggregate exit=$LASTEXITCODE"
W 'SUBSET CGDA INTERIM GATE COMPLETE'
