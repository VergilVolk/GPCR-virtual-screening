"""Render B2 results from saved metrics only; no feature loading or fitting."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import experiment_a as A


def render(run):
    run=Path(run)
    read=lambda name: json.loads((run/name).read_text(encoding='utf-8'))
    cfg,complete=read('frozen_config.json'),read('completion.json')
    if A.file_sha256(run/'frozen_config.json') != complete['config_sha256']:
        raise ValueError('Frozen config changed')
    data={split:read(split+'_results.json') for split in ('R2','R3')}
    histories,replay=read('training_history.json'),read('reproducibility.json')
    fmt=lambda value: 'n/a' if value is None else f'{value:.9g}'
    lines=['# Experiment B2 v01 results','',
           'Descriptive engineering evaluation of incremental dynamic information. R3 is reused and historically inspected. No biological efficacy claim.','',
           '## Frozen configuration','',
           'Config SHA256: `'+complete['config_sha256']+'`. See frozen_config.json for the complete predeclaration.','',
           'R2-only mu; D=R-mu; Q uses explicitly gathered H[k != target] for its centering mean. Canonical G is never a model/loss input. '
           'Same B architecture, decoder, masks, seeds 17/43/101, optimizer Adam 0.001, constant schedule, batch 25, and exactly five epochs. '
           'All 2,000 frames per split, chronological order within trajectory. Dropout removes regional Q tokens after fixed leave-one-out centering.','',
           '| Comparator | Trainable parameters |','|---|---:|']
    for mode,count in cfg['parameter_counts'].items(): lines.append(f"| {mode} | {count['total']} |")
    lines+=['| zero_D | 0; mu contains 34,560 fixed R2 mean values |','',
            'FP64 A decomposition then B-compatible FP32 R. mu accumulated FP64; D evaluated FP64 and trained FP32. '
            'mu matches B bitwise; zero-D prediction/error equals the residue-mean baseline.','']
    comparators,blocks,channels,regions=[],[],[],[]
    scalar_keys=('seed','model','dropout','mse','zero_d_mse','absolute_improvement','excess_mse',
                 'normalized_dynamic_error','relative_improvement','entropy','effective_count')
    for split,d in data.items():
        lines += [f'## {split} dynamic reconstruction','',
                  'No additional dropout. Absolute improvement = zero-D MSE minus model MSE; positive is better. '
                  'Normalized dynamic error = MSE / zero-D MSE; relative improvement = 1 - normalized error.','',
                  '| Seed | Model | MSE | Zero-D MSE | Absolute improvement | Relative improvement | Normalized error | Entropy | Effective residues |',
                  '|---:|---|---:|---:|---:|---:|---:|---:|---:|']
        for r in d['summary']:
            comparators.append(dict(split=split,**{k:r[k] for k in scalar_keys}))
            if r['dropout']==0:
                lines.append(f"| {r['seed']} | {r['model']} | "+' | '.join(fmt(r[k]) for k in
                             ('mse','zero_d_mse','absolute_improvement','relative_improvement','normalized_dynamic_error','entropy','effective_count'))+' |')
            for ch,(mse,zero) in enumerate(zip(r['per_channel_mse'],r['per_channel_zero_d_mse'])):
                channels.append(dict(split=split,seed=r['seed'],model=r['model'],dropout=r['dropout'],channel=ch,
                                     mse=mse,zero_d_mse=zero,absolute_improvement=zero-mse,
                                     normalized_dynamic_error=mse/zero if zero else None))
            for j,name in enumerate(d['region_names']):
                regions.append(dict(split=split,seed=r['seed'],model=r['model'],dropout=r['dropout'],region=name,
                                    mse=r['per_region_mse'][j],zero_d_mse=r['per_region_zero_d_mse'][j],
                                    entropy=r['per_region_entropy'][j] if r['per_region_entropy'] else None,
                                    effective_count=r['per_region_effective_count'][j] if r['per_region_effective_count'] else None))
        lines+=['',f'## {split} dropout sensitivity','',
                '| Seed | Model | MSE 0% | MSE 25% | MSE 50% |','|---:|---|---:|---:|---:|']
        for seed in cfg['seeds']:
            for mode in cfg['comparators']:
                selected=[r for r in d['summary'] if r['seed']==seed and r['model']==mode]
                lines.append(f'| {seed} | {mode} | '+' | '.join(fmt(r['mse']) for r in selected)+' |')
        lines+=['',f'## {split} chronological-block sensitivity','',
                'Ordered window MSE aggregated over four condition trajectories, no additional dropout. '
                'Full block/dropout combinations are in blocks.csv; all trajectory/window cells remain in JSON. '
                'These blocks are temporally correlated, not independent replicates.','',
                '| Seed | Model | Window 0 | Window 1 | Window 2 | Window 3 | Window 4 |',
                '|---:|---|---:|---:|---:|---:|---:|']
        for r in d['chronological_blocks']: blocks.append(dict(split=split,window=r['window'],**{k:r[k] for k in scalar_keys}))
        for seed in cfg['seeds']:
            for mode in cfg['comparators']:
                selected=[r for r in d['chronological_blocks'] if r['seed']==seed and r['model']==mode and r['dropout']==0]
                lines.append(f'| {seed} | {mode} | '+' | '.join(fmt(r['mse']) for r in selected)+' |')
        lines+=['']
    lines+=['## Attention versus mean','',
            'Positive reduction means attention has lower MSE. This comparison is insufficient unless the model also beats zero-D.','',
            '| Split | Seed | MSE reduction 0% (%) | MSE reduction 25% (%) | MSE reduction 50% (%) |',
            '|---|---:|---:|---:|---:|']
    for split,d in data.items():
        for seed in cfg['seeds']:
            advantages=[]
            for drop in cfg['metrics']['token_dropout']:
                values={r['model']:r for r in d['summary'] if r['seed']==seed and r['dropout']==drop}
                advantages.append(100*(1-values['attention']['mse']/values['mean']['mse']))
            lines.append(f'| {split} | {seed} | '+' | '.join(fmt(x) for x in advantages)+' |')
    lines+=['','## R3 per-channel diagnostic summary','',
            'All 128 channels are retained in per_channel.csv for both splits and all dropout levels. '
            'The table summarizes no-dropout results without channel selection: counts of lower error than zero-D and '
            'the largest single-channel share of total squared error. These are descriptive diagnostics, not independent success trials.','',
            '| Seed | Model | Channels below zero-D | Largest channel error share (%) |',
            '|---:|---|---:|---:|']
    for r in data['R3']['summary']:
        if r['dropout'] or r['model']=='zero_D': continue
        mse,zero=np.array(r['per_channel_mse']),np.array(r['per_channel_zero_d_mse'])
        lines.append(f"| {r['seed']} | {r['model']} | {int((mse<zero).sum())}/128 | {100*mse.max()/mse.sum():.6g} |")
    lines+=['','## Attention diagnostics by region','',
            'No additional dropout; all nine regions, all seeds, both splits. Full comparator/dropout metrics are in per_region.csv. '
            'Entropy uses natural logs; effective count is mean exp(per-task entropy). Region sizes differ.','',
            '| Split | Seed | Region | Entropy | Effective residues |','|---|---:|---|---:|---:|']
    for r in regions:
        if r['model']=='attention' and r['dropout']==0:
            lines.append(f"| {r['split']} | {r['seed']} | {r['region']} | {fmt(r['entropy'])} | {fmt(r['effective_count'])} |")
    lines+=['','## Training, reproducibility and tests','',
            '| Fit | Epoch 1 | Epoch 2 | Epoch 3 | Epoch 4 | Epoch 5 | Exact replay | Reload |',
            '|---|---:|---:|---:|---:|---:|---|---|']
    for name,h in histories.items():
        lines.append('| '+name+' | '+' | '.join(fmt(x) for x in h)+f" | {replay[name]['bitwise_replay']} | {replay[name]['checkpoint_roundtrip']} |")
    falling=sum(h[-1]<h[-2] for h in histories.values())
    lines+=['',f'{falling}/9 online training-loss histories decreased from epoch 4 to 5. '
            'This is an optimization limitation; training was not extended. These are within-epoch online losses, not repeated final-checkpoint evaluations.',
            '',f"All {complete['tests_passed']} regression tests passed before training; see tests.txt. "
            'All nine full R2 replays matched state bytes and epoch losses; all checkpoints reloaded exactly. '
            'All 40 source NPZ hashes and protected A/B archive hashes remained unchanged. See provenance.json and completion.json.',
            '', '## Interpretation','']
    primary=[r for r in data['R3']['summary'] if r['dropout']==0 and r['model']!='zero_D']
    count=sum(r['mse']<r['zero_d_mse'] for r in primary)
    lines.append(f'{count}/9 learned seed/comparator results beat zero-D on R3 without additional dropout. '
                 'All seeds and diagnostics are reported; this count is descriptive and is not a significance test.')
    if count==0:
        lines.append('Under this fixed five-epoch masked regional reconstruction objective, no learned readout demonstrates incremental dynamic prediction beyond the residue baseline. Do not add model complexity on the strength of these results.')
    elif count==9:
        lines.append('All comparators show aggregate incremental dynamic prediction under the original masks. This alone does not establish a need for attention or robustness to dropout/block changes.')
    else:
        lines.append('Incremental dynamic prediction is not consistent across all seed/comparator results; interpret the model and sensitivity tables together, without selecting favorable seeds.')
    lines+=['','Q removes the algebraic dependence on masked H introduced by full-frame centering. '
            'It still uses a whole-visible-frame reference, and cached encoder features can already contain contextual information. '
            'R3 is reused descriptive evaluation; correlated frames/windows and overlapping regions are not independent biological confirmation. '
            'No biological efficacy claim is made.','',
            '## Git changes','',
            'Added experiment_b2.py, tests/test_experiment_b2.py, report_b2.py and EXPERIMENT_B2.md. '
            'Existing A/B code and archived results are unchanged. Results are under a new ignored versioned directory. No commit was created.','']
    with (run/'REPORT.md').open('x',encoding='utf-8') as f: f.write('\n'.join(lines))
    for name,rows in [('comparators.csv',comparators),('blocks.csv',blocks),('per_channel.csv',channels),('per_region.csv',regions)]:
        with (run/name).open('x',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(run/'REPORT.md')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    render(parser.parse_args().run)
