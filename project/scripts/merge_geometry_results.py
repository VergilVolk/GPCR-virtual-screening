# -*- coding: utf-8 -*-
"""合并三路候选的几何证据结果 (从 _tmp_geomdock 的 pose 重算, 不重新 docking).

来源:
  fragment top 300  (peamvs_candidates_ranked.csv 前 300)
  lstm    top 100  (ranked CSV 中 source==lstm 前 100)
  gpt     16       (gpt_generated_pam_analogs.csv, 无 QSAR 分数)

输出: results/qsar/peamvs_candidates_geometry_all.csv
"""
import sys, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, r'D:\CLC\project\scripts')
import dock_candidates_geometry as dcg

STRUCT = dcg.STRUCT
OUT = dcg.OUT
TMP = STRUCT / '_tmp_geomdock'
dcg.REC_ATOMS = dcg.load_receptor(STRUCT / '7trq_R_receptor.pdbqt')

# ---- 重建 key -> smiles/source/QSAR ----
df = pd.read_csv(OUT / 'peamvs_candidates_ranked.csv')
frag = df.head(300).copy()
frag['source'] = 'fragment'
lstm = df[df['source'] == 'lstm'].head(100).copy()
lstm['source'] = 'lstm'
gpt = pd.read_csv(r'D:\CLC\project\results\generated\gpt_generated_pam_analogs.csv').head(16).copy()
gpt['source'] = 'gpt'
for c in ['pEC50_pred', 'PAM_prob', 'PEAMVS_score']:
    if c not in gpt.columns:
        gpt[c] = np.nan
meta = pd.concat([frag, lstm, gpt], ignore_index=True)
meta['key'] = meta['smiles'].apply(lambda s: hashlib.md5(s.encode()).hexdigest()[:10])

# ---- 遍历 pose 重算几何特征 ----
rows = []
for out_pdbqt in TMP.glob('*.out.pdbqt'):
    key = out_pdbqt.name.replace('.out.pdbqt', '')
    m = meta[meta['key'] == key]
    if m.empty:
        continue
    lig_atoms = dcg.parse_pose(out_pdbqt)
    if not lig_atoms:
        continue
    feats = dcg.geometric_features(lig_atoms, dcg.REC_ATOMS)
    r = m.iloc[0]
    rows.append({
        'smiles': r['smiles'], 'source': r['source'],
        'pEC50_pred': r.get('pEC50_pred'), 'PAM_prob': r.get('PAM_prob'),
        'PEAMVS_score': r.get('PEAMVS_score'),
        **feats,
    })

g = pd.DataFrame(rows).drop_duplicates(subset='smiles')
# 几何证据合成分
g['z_n_contacts'] = dcg.zscore(g['n_contacts'])
g['z_hub'] = dcg.zscore(g['hub_contacts'])
g['geom_score'] = 0.6 * g['z_n_contacts'] + 0.4 * g['z_hub']
# PEAM-VS v2 融合 (fragment/lstm 有 QSAR; gpt 仅几何)
has_qsar = g['PEAMVS_score'].notna()
g['PEAMVS_v2'] = np.nan
g.loc[has_qsar, 'PEAMVS_v2'] = (0.8 * dcg.zscore(g.loc[has_qsar, 'PEAMVS_score'])
                                + 0.2 * g.loc[has_qsar, 'geom_score'])
g = g.sort_values(['PEAMVS_v2', 'geom_score'], ascending=[False, False], na_position='last')
g['rank_final'] = range(1, len(g) + 1)
g.to_csv(OUT / 'peamvs_candidates_geometry_all.csv', index=False, encoding='utf-8-sig')

summary = {
    'n_total': len(g),
    'by_source': g['source'].value_counts().to_dict(),
    'hub_contact_rate': round(float((g['hub_contacts'] > 0).mean()), 3),
    'mean_n_contacts': round(float(g['n_contacts'].mean()), 1),
    'top15': g.head(15)[['rank_final', 'source', 'PEAMVS_score', 'n_contacts', 'hub_contacts', 'geom_score', 'PEAMVS_v2']].to_dict('records'),
    'gpt_rank': g[g['source'] == 'gpt'][['rank_final', 'n_contacts', 'hub_contacts', 'geom_score']].to_dict('records'),
}
(OUT / 'peamvs_geometry_all_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))
