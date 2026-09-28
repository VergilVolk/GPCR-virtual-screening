"""Render existing frozen Experiment B results; never loads source features."""
import argparse
import csv
import json
from pathlib import Path

import experiment_a as A


def report(run):
    run = Path(run)
    def read(name):
        return json.loads((run / name).read_text(encoding='utf-8'))
    cfg, completion = read('frozen_config.json'), read('completion.json')
    if A.file_sha256(run / 'frozen_config.json') != completion['config_sha256']:
        raise ValueError('Configuration hash mismatch')
    datasets = {split: read(split + '_results.json') for split in ('R2', 'R3')}
    replay = read('reproducibility.json')
    lines = ['# Experiment B v01 results', '',
             'Local residual masked reconstruction only. R3 is historically inspected, descriptive engineering validation.', '',
             'All 40 source archives authenticated through Experiment A before training. No encoder inference or backbone updates.', '',
             'Configuration SHA256: `' + completion['config_sha256'] + '`.', '',
             '## Frozen configuration', '',
             'Seeds 17, 43, 101; all 2,000 R2 frames for fitting, all 2,000 R3 frames for held-out evaluation. '
             'Five chronological windows per replica × condition trajectory; no frame or window random splitting.', '',
             'Unscaled FP64 decomposition with FP32 local training features. Every one of 109 memberships is masked '
             'on every frame. G never enters the local model or objective. Shared decoder family: target embedding '
             '270×8, Linear(136,32), tanh, Linear(32,128). Adam 0.001, constant schedule, five epochs, batch 25; '
             'no early stopping, model selection, or tuning. See frozen_config.json for the full predeclaration.', '',
             '| Comparator | Trainable parameters | Readout parameters |', '|---|---:|---:|']
    for model, counts in cfg['parameter_counts'].items():
        lines.append(f"| {model} | {counts['total']} | {counts['total']-counts['decoder']} |")
    lines += ['| R2 residue mean | 0 (34,560 fitted mean values) | — |', '| zero | 0 | — |', '',
              'The mean predictor is accumulated in FP64 from the same FP32 R2 local features. '
              'Decoder initialization is identical across comparators within each seed.', '']
    rows = []
    for split, data in datasets.items():
        lines += [f'## {split} comparator results', '',
                  'No additional token dropout. NRE = MSE / target mean square; relative-simple = MSE / '
                  'R2-residue-mean predictor MSE on this split. Lower errors are better. Entropy uses natural logs; '
                  'effective count is mean exp(per-task entropy).', '',
                  '| Seed | Comparator | MSE | NRE | Relative-simple | Improvement vs simple | Entropy | Effective residues |',
                  '|---:|---|---:|---:|---:|---:|---:|---:|']
        for row in data['summary']:
            rows.append({'split': split, **row})
            if row['dropout']:
                continue
            fmt = lambda v: '—' if v is None else f'{v:.8g}'
            lines.append(f"| {row['seed']} | {row['model']} | " + ' | '.join(fmt(row[k]) for k in
                         ('mse', 'normalized_error', 'relative_simple', 'improvement_over_simple', 'entropy', 'effective_count')) + ' |')
        lines += ['', f'## {split} dropout sensitivity', '',
                  'Fixed removal of the earliest visible membership tokens; floor(25% or 50% × visible count). '
                  'No retraining. Predictor baselines are unchanged. Full metrics, including entropy under dropout, are in JSON/CSV.', '',
                  '| Seed | Comparator | MSE 0% | MSE 25% | MSE 50% |', '|---:|---|---:|---:|---:|']
        for seed in cfg['seeds']:
            for mode in cfg['comparators']:
                selected = [r for r in data['summary'] if r['seed'] == seed and r['model'] == mode]
                lines.append(f'| {seed} | {mode} | ' + ' | '.join(f"{r['mse']:.8g}" for r in selected) + ' |')
        lines += ['', f'## {split} chronological-block sensitivity', '',
                  'MSE by ordered window, pooled over the four condition trajectories; no additional dropout. '
                  'Blocks are temporally correlated, not independent replicates. JSON additionally contains '
                  'every block/dropout combination and each trajectory/window cell.', '',
                  '| Seed | Comparator | Window 0 | Window 1 | Window 2 | Window 3 | Window 4 |',
                  '|---:|---|---:|---:|---:|---:|---:|']
        for seed in cfg['seeds']:
            for mode in cfg['comparators']:
                selected = [r for r in data['chronological_blocks'] if r['seed'] == seed and r['model'] == mode and r['dropout'] == 0]
                lines.append(f'| {seed} | {mode} | ' + ' | '.join(f"{r['mse']:.8g}" for r in selected) + ' |')
        lines.append('')
    lines += ['## Observed comparison', '',
              'All seeds are shown below. Positive percentages mean attention reduced MSE relative to mean pooling. '
              'These are descriptive differences, not tests of statistical significance.', '',
              '| Split | Seed | Attention MSE reduction vs mean (%) | Attention / simple-predictor MSE |',
              '|---|---:|---:|---:|']
    for split, data in datasets.items():
        for seed in cfg['seeds']:
            entries = {r['model']: r for r in data['summary'] if r['seed'] == seed and r['dropout'] == 0}
            reduction = 100 * (1 - entries['attention']['mse'] / entries['mean']['mse'])
            lines.append(f"| {split} | {seed} | {reduction:.6g} | {entries['attention']['relative_simple']:.6g} |")
    lines += ['', 'The R2-residue-mean predictor has lower MSE than every trained decoder on both splits. '
              'Attention modestly improves over the matched mean-pooling decoder in this run, but does not '
              'outperform the simple predictor. Static weighting is worse than mean pooling on R3 for all three seeds. '
              'This fixed-budget experiment does not establish biological usefulness or the best attainable '
              'performance of any readout.', '',
              '## Reproducibility and tests', '',
              '| Fit | Exact full training replay | Exact loss history | Checkpoint roundtrip |', '|---|---|---|---|']
    for name, result in replay.items():
        lines.append(f"| {name} | {result['bitwise_replay']} | {result['loss_history_equal']} | {result['checkpoint_roundtrip']} |")
    if len({r['state_sha256'] for r in replay.values()}) != len(replay):
        raise ValueError('Unexpected duplicate states across fits')
    lines += ['', 'All nine fitted states have distinct hashes. All three predefined seeds are included. '
              'Replay checks cover this CPU/runtime, not cross-platform bitwise equivalence.', '',
              'Before real training: `python -B -m unittest discover -s project/encoder_readout_v01/tests -v` '
              'passed all 22 tests (11 Experiment A regressions and 11 Experiment B synthetic contracts). '
              'B tests cover uniform equivalence, normalized weights, masked-target absence and zero target gradient, '
              'G exclusion, decomposition, seed determinism, ordering/indexing, malformed inputs, identical decoder '
              'initialization, split rejection and trajectory/window constraints.', '',
              '## Scope and limitations', '',
              'Interpret these metrics only as local residual readout information preservation under the frozen objective '
              'and training budget. Shared full-frame centering creates dependence between visible R and masked H; '
              'the direct target R token and the separate G branch are excluded. Temporally correlated frames, '
              'overlapping memberships and historically inspected R3 do not provide independent biological confirmation. '
              'The five-epoch budget may underfit. No biological efficacy claim is made.', '',
              '## Git change summary', '',
              'Added experiment_b.py, tests/test_experiment_b.py, EXPERIMENT_B.md and report_b.py. '
              'Experiment A and all pre-existing source files are unchanged. Results/checkpoints are in the new '
              'versioned result directory, which is ignored by repository policy. No commit was created.', '']
    with (run / 'REPORT.md').open('x', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    with (run / 'comparators.csv').open('x', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(run / 'REPORT.md')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    report(parser.parse_args().run)
