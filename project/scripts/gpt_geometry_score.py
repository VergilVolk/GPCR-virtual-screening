# -*- coding: utf-8 -*-
"""GPT 生成 16 候选: 几何接触证据打分 (无 QSAR 分数, 仅几何证据排名)."""
import sys
sys.path.insert(0, r'D:\CLC\project\scripts')
import pandas as pd
import dock_candidates_geometry as dcg

dcg.REC_ATOMS = dcg.load_receptor(dcg.STRUCT / '7trq_R_receptor.pdbqt')
tmp = dcg.STRUCT / '_tmp_geomdock'
tmp.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(r'D:\CLC\project\results\generated\gpt_generated_pam_analogs.csv')
rows = []
for i, r in enumerate(df.itertuples()):
    aff, feats = dcg.prepare_and_dock(r.smiles, tmp)
    if feats is None:
        print(f'  [{i+1}/{len(df)}] FAIL {r.smiles[:50]}', flush=True)
        continue
    rows.append({'smiles': r.smiles, 'source': 'gpt',
                 'vina_affinity_ref': round(aff, 2) if aff else None, **feats})
    print(f'  [{i+1}/{len(df)}] n_contacts={feats["n_contacts"]} hub={feats["hub_contacts"]} {r.smiles[:50]}', flush=True)

g = pd.DataFrame(rows)
g['z_n_contacts'] = dcg.zscore(g['n_contacts'])
g['z_hub'] = dcg.zscore(g['hub_contacts'])
g['geom_score'] = 0.6 * g['z_n_contacts'] + 0.4 * g['z_hub']
g = g.sort_values('geom_score', ascending=False).reset_index(drop=True)
g['rank_geom'] = g.index + 1
g.to_csv(r'D:\CLC\project\results\qsar\gpt_geometry.csv', index=False, encoding='utf-8-sig')
print('saved', len(g))
print(g[['rank_geom', 'n_contacts', 'hub_contacts', 'geom_score']].head(16).to_string())
