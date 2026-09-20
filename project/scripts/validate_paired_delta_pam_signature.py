"""Probe-matched paired-delta external validation of an M4 PAM structural signature."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import test_consensus_pam_signature as sigmod
from validate_external_static_pam_signature import ca_coords


PROJECT = Path(__file__).resolve().parents[1]
EXT = PROJECT / "data" / "pdb" / "m4_external_structures"
ENS = PROJECT / "results" / "structure" / "ensemble"
OUT = PROJECT / "results" / "pacer_paired_delta_pam_signature_v04"


def distances(path, chain):
    return sigmod.pair_features(ca_coords(path, chain))


def model(q, p, k, threshold):
    dq, dp = q-k, p-k
    mask = (dq*dp > 0) & (np.abs(dq) >= threshold) & (np.abs(dp) >= threshold)
    d = (dq+dp)/2
    agree = np.minimum(np.abs(dq), np.abs(dp))/np.maximum(np.abs(dq), np.abs(dp))
    v = d[mask]*np.sqrt(agree[mask]); denom=float(v@v)
    def score(x): return float(((x[mask]-k[mask])*np.sqrt(agree[mask]))@v/denom)
    return score, mask


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    q = sigmod.pair_features(np.stack([v[1] for _,v in sorted(sigmod.static_ca(ENS/'7TRQ_R.pdb').items())
                                       if _ in sigmod.STATIC_RESIDS]))
    p = sigmod.pair_features(np.stack([v[1] for _,v in sorted(sigmod.static_ca(ENS/'7TRP_R.pdb').items())
                                       if _ in sigmod.STATIC_RESIDS]))
    structures = {
        '7TRK_iperoxo_only': distances(EXT/'7TRK.pdb','R'),
        '7TRQ_VU0467154_PAM': q,
        '7TRP_LY2033298_PAM': p,
        '7V69_iperoxo_only_external': distances(EXT/'7V69.pdb','R'),
        '7V68_LY2119620_PAM_external': distances(EXT/'7V68.pdb','R'),
        '7V6A_compound110_allosteric_agonist': distances(EXT/'7V6A.pdb','R'),
    }
    k=structures['7TRK_iperoxo_only']
    rows=[]; checks=[]
    for threshold in (.05,.10,.20):
        score,mask=model(q,p,k,threshold)
        vals={name:score(x) for name,x in structures.items()}
        for name,val in vals.items(): rows.append({'threshold_A':threshold,'structure':name,'score':val,'features':int(mask.sum())})
        delta=vals['7V68_LY2119620_PAM_external']-vals['7V69_iperoxo_only_external']
        specific=vals['7V68_LY2119620_PAM_external']>vals['7V6A_compound110_allosteric_agonist']
        checks.append({'threshold_A':threshold,'external_paired_delta':delta,
                       'paired_delta_positive':delta>0,'PAM_above_allosteric_agonist':specific})
    table=pd.DataFrame(rows); check=pd.DataFrame(checks)
    table.to_csv(OUT/'scores.csv',index=False); check.to_csv(OUT/'checks.csv',index=False)
    audit={'probe_matched_development':True,'independent_paired_delta_support':bool(check.paired_delta_positive.all()),
           'PAM_specificity_support':bool(check.PAM_above_allosteric_agonist.all()),
           'claim_boundary':'Probe-matched static PAM increment; no efficacy inference.'}
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    (OUT/'REPORT.md').write_text('# PACER Paired-Delta PAM Signature v0.4\n\n'+table.to_markdown(index=False,floatfmt='.4f')+
        '\n\n'+check.to_markdown(index=False,floatfmt='.4f')+'\n\n```json\n'+json.dumps(audit,indent=2)+'\n```\n',encoding='utf-8')
    print(table.to_string(index=False));print(check.to_string(index=False));print(json.dumps(audit,indent=2))


if __name__=='__main__': main()
