"""Recover byte-exact author reference blobs at the already documented frozen commit."""
import hashlib, json, urllib.request
from pathlib import Path

P = Path(__file__).resolve().parents[1]
COMMIT = '44798c841ee77230b1f89fe41074e41b458c7070'
ROOT = P / 'results/pacer_xr_rerun_v01/references'
FILES = {
    'Glide_PDB': ('PDB/M4R_PDB_Glide_scores.csv', 'glide_gscore', 'e53bc82a61c3d808c00d0418d5e837922c881d6a'),
    'Glide_BEmin': ('Ensemble/M4R_Ensemble_Glide_BEmin_ranked.csv', 'BE_min', 'e4c00d306c7f14e3e6a8600408c243a164552693'),
    'Glide_BEavg': ('Ensemble/M4R_Ensemble_Glide_BEavg_ranked.csv', 'BE_avg', '454f531b32f16e6a04ba260c37c4c71671c1dd04'),
    'Vina_PDB': ('PDB/M4R_PDB_Vina_scores.csv', 'vina_score', 'c35622bb512049579eb826e5f786686b1dc8c3fe'),
    'Vina_BEmin': ('Ensemble/M4R_Ensemble_Vina_BEmin_ranked.csv', 'BE_min', '9e6d585d2df4e591103d8ec7d7c61845c73a15ee'),
    'Vina_BEavg': ('Ensemble/M4R_Ensemble_Vina_BEavg_ranked.csv', 'BE_avg', '911d5c33f54877c4357e8a44d07a2e74612edb22'),
}

def main():
    import numpy as np
    import pandas as pd
    records = {}
    for channel, (relative, column, blob) in FILES.items():
        target = ROOT / relative
        url = f'https://raw.githubusercontent.com/tylerdt1/gpcr-am-ensemble-docking/{COMMIT}/docking_scores/M4R/{relative}'
        data = target.read_bytes() if target.exists() else urllib.request.urlopen(url, timeout=180).read()
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if actual != blob:
            raise ValueError(f'REFERENCE_RESTORE_BLOCKED: {channel} Git blob mismatch')
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists(): target.write_bytes(data)
        frame = pd.read_csv(target)
        values = pd.to_numeric(frame[column], errors='raise').to_numpy(float)
        finite = values[np.isfinite(values)]
        if not len(finite): raise ValueError(f'Empty reference: {channel}')
        records[channel] = dict(path=str(target.relative_to(P)).replace('\\','/'), column=column,
            sha256=hashlib.sha256(data).hexdigest(), git_blob=blob, rows=len(values), finite_rows=len(finite),
            nonfinite_rows=len(values)-len(finite), source_url=url, source_commit=COMMIT)
        print(channel + ' FOUND + VERIFIED ' + records[channel]['sha256'], flush=True)
    lock = dict(reference_source='tylerdt1/gpcr-am-ensemble-docking', reference_commit=COMMIT,
        provenance_document='docs/THOMPSON_MIAO_2026_EXACT_REPRODUCTION.md', channels=records)
    target = P / 'config/pacer_xr_references_v01.json'
    if target.exists() and json.loads(target.read_text()) != lock:
        raise ValueError('Existing frozen reference lock differs; refusing overwrite')
    target.write_text(json.dumps(lock, indent=2)+'\n', encoding='utf-8')

if __name__ == '__main__': main()
