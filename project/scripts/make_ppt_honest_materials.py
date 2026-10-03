"""诚实版 PPT 材料生产线：模块2/模块3 真实输出 + 四张图。

数据源（全部来自事故分支归档，逐字节可溯）：
  - corrected_v1 双模型路由表（SHA bda22c93，权威）
  - framewise_scores.csv（2,000 任务 × 每构象 Vina 分）
  - structural_gate / diverse_shortlist（模块3 冻结输出）
  - 生成库 generated_pam_analogs.csv（2,605 分子）
"""
import subprocess, io, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BR = 'origin/incident/pacer-xr-protocol-deviation-20261003'
OUT = Path(r'D:\CLC\project\results\pacer_ppt_honest_v01')
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'axes.grid': True,
                     'grid.alpha': .3, 'axes.spines.top': False, 'axes.spines.right': False,
                     'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'DejaVu Sans'],
                     'axes.unicode_minus': False})

def show(path):
    # 长路径下 git show 失效：先 ls-tree 拿 blob hash，再 cat-file
    h = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', BR, '--', path],
                       capture_output=True, text=True).stdout.split()
    if not h:
        raise RuntimeError('no ls-tree entry: ' + path)
    r = subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h[2]], capture_output=True)
    if not r.stdout:
        raise RuntimeError('empty blob: ' + path)
    return r.stdout

E = 'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening/project/results'

routed = pd.read_csv(io.BytesIO(show(f'{E}/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv')))
framewise = pd.read_csv(io.BytesIO(show(f'{E}/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/framewise_scores.csv')))
gate = pd.read_csv(io.BytesIO(show(f'{E}/pacer_stage3d_corrected_20261002_v01/pacer200_structural_gate.csv')))
short = pd.read_csv(io.BytesIO(show(f'{E}/pacer_stage3d_corrected_20261002_v01/pacer_stage3d_diverse_shortlist.csv')))
gen = pd.read_csv(io.BytesIO(show('project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening-final/project/results/generated/generated_pam_analogs.csv')))

print('routed', routed.shape, '| framewise', framewise.shape, '| gate', gate.shape,
      '| short', short.shape, '| gen', gen.shape)
routed.to_csv(OUT / 'module2_routed200_full.csv', index=False)

# ============ 模块2 输出：top-100 库（真实 M4-safe 绑定分） ============
top100 = routed.sort_values('final_rank').head(100)[
    ['pair_id', 'canonical_smiles', 'score_2023_gpcr_loto', 'score_2026_13t_famaug',
     'pacer_binding_score', 'final_rank']].copy()
top100.to_csv(OUT / 'module2_top100_library.csv', index=False)

# ============ 模块3 输出：真实结构门 + 140 骨架库 + top12 切片 ============
# short 表已自带全部 gate 指标，无需再并 gate
m3 = short.copy()
m3.to_csv(OUT / 'module3_diverse140_with_gates.csv', index=False)
m3.head(12).to_csv(OUT / 'module3_top12_slice.csv', index=False)

FINAL3 = ['PACER0010', 'PACER0073', 'PACER0027']
f3 = m3[m3.candidate_id.isin(FINAL3)].sort_values('diverse_rank')
f3.to_csv(OUT / 'final3_candidates.csv', index=False)

C = '#2563eb'; C2 = '#dc2626'; C3 = '#059669'

# ---- 图1：生成分子分析（2,605）----
fig, ax = plt.subplots(2, 4, figsize=(11, 5))
for a, col, unit in zip(ax[0], ['mw', 'logp', 'tpsa', 'aromatic_rings'],
                        ['MW (Da)', 'logP', 'TPSA (Å²)', '芳香环数']):
    a.hist(gen[col].dropna(), bins=40, color=C, alpha=.85)
    a.set_xlabel(unit); a.set_ylabel('分子数')
ax[0, 0].set_title('分子量', fontsize=10)
ax[0, 1].set_title('脂溶性 logP', fontsize=10)
ax[0, 2].set_title('极性表面积 TPSA', fontsize=10)
ax[0, 3].set_title('芳香环数', fontsize=10)
for a, col, unit in zip(ax[1], ['hba', 'hbd', 'aromatic_rings', 'tpsa'],
                        ['H 键受体数', 'H 键供体数', '', '']):
    if unit:
        a.hist(gen[col].dropna(), bins=30, color=C3, alpha=.85)
        a.set_xlabel(unit); a.set_ylabel('分子数')
ax[1, 2].axis('off'); ax[1, 3].axis('off')
ax[1, 2].text(0.02, 0.75,
              'BRICS 片段重组生成\n17,146 次 join 尝试\n5,000 唯一分子\n2,605 通过生成规则\n(2,395 未通过)',
              fontsize=11, va='top',
              bbox=dict(boxstyle='round', fc='#f0f9ff', ec=C, alpha=.9))
ax[1, 3].text(0.02, 0.95,
              '药化过滤（Stage2）\n排除 2,405：\n13 已知重复 / 931 域违规\n15 风险 / 1,446 配额\n→ 200 入栈\n(local 160 + exploratory 40)',
              fontsize=9.5, va='top',
              bbox=dict(boxstyle='round', fc='#fff7ed', ec='#ea580c', alpha=.9))
fig.suptitle('模块1 输出：M4 PAM 候选生成与药化筛选（2,605 个生成分子性质分布）', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, .95])
fig.savefig(OUT / 'fig1_generation_analysis.png'); plt.close(fig)

# ---- 图2：模块2 双模型路由 ----
fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
x = routed.score_2023_gpcr_loto; y = routed.score_2026_13t_famaug
sel = routed.final_rank <= 100
ax[0].scatter(x[~sel], y[~sel], s=14, c='#94a3b8', alpha=.6, label='未入选 (100)')
ax[0].scatter(x[sel], y[sel], s=18, c=C, alpha=.85, label='模块2 输出 top-100')
for pid in FINAL3:
    r = routed[routed.pair_id == pid]
    xs, ys = float(r.score_2023_gpcr_loto.iloc[0]), float(r.score_2026_13t_famaug.iloc[0])
    ax[0].scatter(xs, ys, s=60, c=C2, marker='*')
    ax[0].annotate(pid, (xs, ys), fontsize=7.5, xytext=(3, 3), textcoords='offset points')
ax[0].set_xlabel('2023 基座 M4-safe 结合分（主路由）'); ax[0].set_ylabel('Science2026 famaug 分（第二意见）')
ax[0].legend(fontsize=7.5, loc='lower left'); ax[0].set_title('双模型打分（200 分子）', fontsize=10)
ax[1].hist(routed.pacer_binding_score, bins=40, color=C, alpha=.85)
ax[1].set_xlabel('M4-safe 绑定分（百分位）'); ax[1].set_ylabel('分子数')
ax[1].set_title('绑定分分布与 top-100 截断', fontsize=10)
ax[1].axvline(top100.pacer_binding_score.min(), color=C2, ls='--', lw=1.2,
              label=f'截断 = {top100.pacer_binding_score.min():.3f}')
ax[1].legend(fontsize=8)
r = routed.sort_values('final_rank')
ax[2].plot(range(1, 201), r.score_2023_gpcr_loto.values, '-', color=C, lw=1.4)
ax[2].axvline(100, color=C2, ls='--', lw=1.2)
ax[2].set_xlabel('M4-safe 排名'); ax[2].set_ylabel('2023 结合分')
ax[2].set_title('排名曲线（top-100 = 模块2 输出库）', fontsize=10)
fig.suptitle('模块2：基于 DrugCLIP 迁移学习的 M4 智能靶点筛选（双权重路由）', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(OUT / 'fig2_module2_routing.png'); plt.close(fig)

# ---- 图3：模块3 系综对接 + 结构门 + 漏斗 ----
fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
ax[0].hist(gate.best_vina, bins=40, color=C, alpha=.85, label='best (BEmin 型)')
ax[0].hist(gate.mean_vina, bins=40, color=C3, alpha=.65, label='mean (BEavg 型)')
ax[0].set_xlabel('Vina 结合能 (kcal/mol)'); ax[0].set_ylabel('分子数')
ax[0].legend(fontsize=8); ax[0].set_title('200 分子 × 10 GaMD 构象系综打分', fontsize=10)
ax[1].scatter(gate.median_native_pocket_coverage, gate.best_vina, s=15, c=C, alpha=.7)
ax[1].axvline(0.5, color=C2, ls='--', lw=1.2); ax[1].axhline(-8, color=C2, ls=':', lw=1)
ax[1].set_xlabel('原生口袋覆盖中位数'); ax[1].set_ylabel('best Vina (kcal/mol)')
ax[1].set_title(f'结构门四项全部通过：{len(gate)}/200 PASS', fontsize=10)
stages = ['生成通过\n2,605', '入栈\n200', '系综对接\n2,000 任务\n200 PASS', '骨架代表\n140', '进入 MD\n3']
nums = [2605, 200, 200, 140, 3]
ax[2].bar(range(5), np.log10(nums), color=[C3, C, C, C, C2], alpha=.85)
for i, (n, s) in enumerate(zip(nums, stages)):
    ax[2].text(i, np.log10(n) + .08, f'{n}', ha='center', fontsize=10, fontweight='bold')
ax[2].set_xticks(range(5)); ax[2].set_xticklabels(stages, fontsize=7.5)
ax[2].set_ylabel('log₁₀ 分子数'); ax[2].set_title('虚筛漏斗（每级数字均可溯源）', fontsize=10)
fig.suptitle('模块3：静态-动态构象组合分子对接（Vina × 10 GaMD 构象 + 四项结构门）', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(OUT / 'fig3_module3_ensemble.png'); plt.close(fig)

# ---- 图4：最终三候选 ----
from rdkit import Chem
from rdkit.Chem import Draw
rows = f3 if len(f3) else m3[m3.candidate_id.isin(FINAL3)]
mols, legends = [], []
for _, r_ in rows.iterrows():
    m = Chem.MolFromSmiles(r_.canonical_smiles)
    if m:
        mols.append(m)
        legends.append(f"{r_.candidate_id}  绑定排名#{int(r_.final_rank)}  骨架代表#{int(r_.diverse_rank)}\n"
                       f"2023分 {r_.score_2023_gpcr_loto:.3f} | best Vina {r_.best_vina:.1f} kcal/mol")
if mols:
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(420, 320), legends=legends)
    img.save(OUT / 'fig4_final3_molecules.png')

print('figures saved to', OUT)
print('final3:'); print(f3[['candidate_id', 'final_rank', 'diverse_rank',
                            'score_2023_gpcr_loto', 'best_vina']].to_string(index=False))

