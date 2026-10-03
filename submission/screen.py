#!/usr/bin/env python3
"""一键 M4 虚拟筛选：分子/口袋表征 + 投影权重 -> 标准化候选清单 results.csv

演示模式（自带 50 分子示例数据，无任何外部依赖文件）：
    python screen.py --demo

全量模式（主仓库 28,519 分子库表征）：
    python screen.py --full --clc-root D:/CLC

输出：results/results.csv（附件5第四节标准字段）。
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))

from src.screen_core import screen, write_results_csv

MODEL_TAG = "DrugCLIP(2023基座)+muscarinic-triplet投影 v01 (external_result.seed20260924-26)"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--demo', action='store_true', help='使用自带 50 分子示例')
    g.add_argument('--full', action='store_true', help='使用主仓库全量 28,519 分子库')
    p.add_argument('--clc-root', type=Path, default=Path('D:/CLC'),
                   help='主仓库根目录（--full 模式用）')
    p.add_argument('--pooling', choices=['mean', 'max'], default='mean')
    p.add_argument('--top-n', type=int, default=50)
    a = p.parse_args()

    models = sorted((ROOT / 'models').glob('external_result.seed*.projection.pt'))
    if not models:
        sys.exit('models/ 下未找到 external_result 权重')

    if a.demo:
        mol_npz = ROOT / 'data' / 'demo' / 'demo_molecules.npz'
        pocket_npz = ROOT / 'data' / 'demo' / 'demo_m4_pockets.npz'
    else:
        r = a.clc_root / 'project' / 'results' / 'gpcr_drugclip_screening_v01'
        mol_npz = r / 'ensemble_embeddings.npz'
        pocket_npz = r / 'ensemble_embeddings.npz'

    ids, scores, per_seed = screen(mol_npz, pocket_npz, models, a.pooling)
    out = write_results_csv(ROOT / 'results' / 'results.csv', ids, scores,
                            MODEL_TAG, a.top_n)
    print(f'screened {len(ids)} molecules x {len(models)} seeds ({a.pooling}-pool)')
    print(f'results -> {out}')
    top = np.argsort(-scores)[:5]
    for r, i in enumerate(top, 1):
        print(f'  {r}. {ids[i][:58]:<58} {scores[i]:.4f}')


if __name__ == '__main__':
    main()
