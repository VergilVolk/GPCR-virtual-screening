# LIT-PCBA remaining-targets watchdog: relaunch on silent death, exit when all 15 done
$wd = 'D:\CLC\project'
$dir = "$wd\results\drugclip_science2026\full_litpcba"
$targets = 'OPRK1,VDR,FEN1,GBA,IDH1,KAT2A,ADRB2'
while ($true) {
  $done = (Get-ChildItem $dir -Filter '*.metrics.json' -ErrorAction SilentlyContinue).Count
  if ($done -ge 15) { Add-Content "$dir\watchdog.log" "$(Get-Date -Format 'MM-dd HH:mm:ss') all 15 complete"; break }
  $running = Get-CimInstance Win32_Process -Filter "Name='python.exe'" -ErrorAction SilentlyContinue | Where-Object { $_.CommandLine -match 'run_drugclip_science2026_full_litpcba' }
  if (-not $running) {
    $p = Start-Process -FilePath python -ArgumentList @(
      'scripts/run_drugclip_science2026_full_litpcba.py',
      '--drugclip','tools/DrugCLIP','--unicore','tools/Uni-Core',
      '--checkpoint','data/external/drugclip_science2026/litpcba_identity_90.pt',
      '--data-root','data/external/drugclip_science2026/LIT-PCBA/lit_pcba',
      '--source-zip','data/external/drugclip_science2026/LIT-PCBA.zip',
      '--output','results/drugclip_science2026/full_litpcba',
      '--trusted-checkpoint','--targets',$targets
    ) -WorkingDirectory $wd -WindowStyle Hidden `
      -RedirectStandardOutput "$dir\remaining_stdout.log" `
      -RedirectStandardError "$dir\remaining_stderr.log" -PassThru
    Add-Content "$dir\watchdog.log" "$(Get-Date -Format 'MM-dd HH:mm:ss') relaunched PID $($p.Id) done=$done/15"
  }
  Start-Sleep -Seconds 300
}
