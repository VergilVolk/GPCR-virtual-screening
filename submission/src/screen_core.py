"""M4 虚拟筛选核心逻辑。

流程（附件5 要求的 推理/筛选 环节）：
    1. 载入分子 / 口袋的 DrugCLIP 512 维表征（npz）；
    2. 载入 3 个随机种子的投影头权重（2023 基座 triplet 线，M4 同测试集
       对决胜出配置：ROC 0.7136 / BEDROC20 0.1802 / PR 0.1547）；
    3. 每种子对 10 个 GaMD M4 口袋构象打分后取 mean（head-to-head 实测
       mean 全面优于 max，见主仓库 m4r_head2head_2023vs2026.json）；
    4. 三种子等权平均得最终分数，降序输出候选清单。
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import torch

from src.proj import Proj

TRACK = "AI模型与计算方法-虚拟筛选(CHRM4/M4R正构-变构候选)"


def load_npz(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    return {k: z[k] for k in z.files}


def screen(mol_npz: Path, pocket_npz: Path, model_paths: list[Path],
           pooling: str = "mean") -> tuple[list[str], np.ndarray, np.ndarray]:
    """返回 (molecule_ids, 集成分数, 每种子分数[n_seed, n_mol])。"""
    mols, pockets = load_npz(mol_npz), load_npz(pocket_npz)
    ids = [str(v) for v in mols['molecule_ids']]
    mol = torch.as_tensor(mols['molecule_representations'].astype(np.float32))
    pocket = torch.as_tensor(pockets['pocket_representations'].astype(np.float32))
    per_seed = []
    with torch.inference_mode():
        for mp in model_paths:
            d = torch.load(mp, map_location='cpu', weights_only=False)
            pm, pp = Proj(d['mol_project']).eval(), Proj(d['pocket_project']).eval()
            sc = (pm(mol) @ pp(pocket).T).numpy()          # n_mol x n_pockets
            per_seed.append(sc.max(axis=1) if pooling == 'max' else sc.mean(axis=1))
    per_seed = np.stack(per_seed)
    return ids, per_seed.mean(axis=0), per_seed


def write_results_csv(out_path: Path, ids: list[str], scores: np.ndarray,
                      model_tag: str, top_n: int = 50) -> Path:
    """按附件5第四节输出标准化候选清单（UTF-8 results.csv）。"""
    import csv
    order = np.argsort(-scores)[:top_n]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['候选编号', '所属赛道', 'SMILES', '关键预测指标_结合分数(余弦)',
                    '对应模型与运行版本', '备注'])
        for rank, i in enumerate(order, 1):
            w.writerow([f'CAND-{rank:04d}', TRACK, ids[i], f'{scores[i]:.6f}',
                        model_tag,
                        '10xGaMD M4口袋mean池化; 3种子(20260924-26)集成; '
                        '结合检索候选, 非PAM功能结论, 需四上下文MD复核'])
    return out_path
