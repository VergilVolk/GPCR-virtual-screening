$ErrorActionPreference = 'Stop'
Set-Location 'D:\CLC'

$runRoot = 'D:\CLC\project\results\m4_gamd_ensemble'
$statusFile = Join-Path $runRoot 'BACKGROUND_RUN_STATUS.txt'
$doneFile = Join-Path $runRoot 'BACKGROUND_RUN_COMPLETE.txt'

function Invoke-Step {
    param([string]$Name, [string[]]$Arguments)
    $stamp = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'
    Set-Content -LiteralPath $statusFile -Value "$stamp RUNNING $Name"
    & python @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Name failed with exit code $LASTEXITCODE"
    }
}

try {
    Invoke-Step 'primary_validation_docking' @(
        'project/scripts/dock_m4_gamd_ensemble.py', '--set', 'validation',
        '--validation-scope', 'primary', '--clusters', '0-9', '--workers', '8'
    )
    Invoke-Step 'primary_validation_analysis' @('project/scripts/analyze_m4_gamd_ensemble.py')
    Invoke-Step 'candidate_ensemble_docking' @(
        'project/scripts/dock_m4_gamd_ensemble.py', '--set', 'candidates',
        '--clusters', '0-9', '--workers', '8'
    )
    Invoke-Step 'candidate_pareto_ranking' @('project/scripts/rank_m4_gamd_candidates.py')
    $stamp = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'
    Set-Content -LiteralPath $statusFile -Value "$stamp COMPLETE"
    Set-Content -LiteralPath $doneFile -Value "$stamp COMPLETE"
}
catch {
    $stamp = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'
    Set-Content -LiteralPath $statusFile -Value "$stamp FAILED $($_.Exception.Message)"
    throw
}
