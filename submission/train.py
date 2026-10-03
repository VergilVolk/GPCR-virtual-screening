#!/usr/bin/env python3
"""一键训练入口（附件5 二-3）：微调 M4 筛选投影权重。

两层入口设计：
  --full-data   在主仓库全部数据上训练部署权重（复现 models/ 下的权重）
  --loso        13 靶留出评估（复现论文基准数字，耗时较长）

训练实际执行主仓库脚本（保持单一事实源）：
  project/scripts/evaluate_drugclip_loso_family_aug.py
"""
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).parent
CLC = Path(__file__).resolve().parents[1]          # 主仓库根（提交包位于其下）

SCRIPT = CLC / 'project' / 'scripts' / 'evaluate_drugclip_loso_family_aug.py'
R = CLC / 'project' / 'results'
D = CLC / 'project' / 'data' / 'benchmarks'
COMMON = [
    '--gpcr-representations', str(R / 'gpcr_drugclip_screening_v01' / 'science2026_ensemble_embeddings.npz'),
    '--gpcr-pairs', str(D / 'gpcr_drugclip_screening_v01' / 'screening_pairs.csv'),
    '--external-representations', str(R / 'drugclip_science2026' / 'litpcba_external_v01' / 'science2026_90.npz'),
    '--external-pairs', str(D / 'litpcba_drugclip_external_v01' / 'pairs.csv'),
    '--projection', str(R / 'drugclip_science2026' / 'litpcba_external_v01' / 'science2026_90.projection.pt'),
    '--family-npz', str(D / 'drugclip_muscarinic_family_aug_v01' / 'family_science2026.npz'),
    '--family-pairs', str(D / 'drugclip_muscarinic_family_aug_v01' / 'family_pairs_clean.csv'),
]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--full-data', action='store_true', help='全量数据训练部署权重')
    g.add_argument('--loso', action='store_true', help='13 靶留出基准评估')
    p.add_argument('--seed', type=int, default=20260925)
    p.add_argument('--epochs', type=int, default=80)
    a = p.parse_args()
    if not SCRIPT.exists():
        sys.exit(f'主仓库脚本不存在: {SCRIPT}（提交包需位于主仓库内运行）')
    out = ROOT / 'logs'
    out.mkdir(exist_ok=True)
    if a.full_data:
        cmd = ['python', str(SCRIPT), *COMMON, '--seed', str(a.seed), '--epochs', str(a.epochs),
               '--full-only', '--save-full-checkpoint',
               str(out / f'retrained_seed{a.seed}.projection.pt')]
    else:
        cmd = ['python', str(SCRIPT), *COMMON, '--seed', str(a.seed), '--epochs', str(a.epochs),
               '--output', str(out / f'loso_seed{a.seed}.json')]
    print(' '.join(cmd), flush=True)
    sys.exit(subprocess.call(cmd))


if __name__ == '__main__':
    main()
