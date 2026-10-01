param()
$out = 'D:\CLC\project\results\drugclip_science2026\migration_13target_v01'
python D:\CLC\project\scripts\summarize_drugclip_13target_multiseed.py --prediction "$out\loso_preserve0.2_seed20260925.predictions.csv" --prediction "$out\loso_preserve0.2_seed20260926.predictions.csv" --prediction "$out\loso_preserve0.2_seed20260927.predictions.csv" --ecfp-predictions "$out\ecfp_loso.predictions.csv" --output "$out\summary_preserve0.2_3seed.json" > "$out\summary_p0.2.stdout.log" 2> "$out\summary_p0.2.stderr.log"
"summary p0.2 exit=$LASTEXITCODE" | Add-Content "$out\runner.log"
python D:\CLC\project\scripts\summarize_drugclip_13target_multiseed.py --prediction "$out\loso_preserve0_seed20260925.predictions.csv" --prediction "$out\loso_preserve0_seed20260926.predictions.csv" --prediction "$out\loso_preserve0_seed20260927.predictions.csv" --ecfp-predictions "$out\ecfp_loso.predictions.csv" --output "$out\summary_preserve0_3seed.json" > "$out\summary_p0.stdout.log" 2> "$out\summary_p0.stderr.log"
"summary p0 exit=$LASTEXITCODE" | Add-Content "$out\runner.log"
"SUMMARIES COMPLETE" | Add-Content "$out\runner.log"
