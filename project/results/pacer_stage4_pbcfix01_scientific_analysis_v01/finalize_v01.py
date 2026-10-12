"""Verify delivered tables against block rows and record complete output hashes."""
from pathlib import Path
import csv,json,hashlib,subprocess,itertools
import numpy as np
from PIL import Image
P=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(n): return list(csv.DictReader((P/n).open(encoding='utf-8-sig')))
def dump(n,d): (P/n).write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
version=json.loads((P/'ANALYSIS_CODE_VERSION_v01.json').read_text())
assert version['script_sha256']==sha(P/'analyze_v01.py')
assert version['report_renderer_sha256']==sha(P/'write_report_v01.py')
inputs=read('INPUT_SHA256_v01.csv')
assert all(sha(Path(x['path']))==x['sha256'] for x in inputs)
regions=read('STAGE4_PBCFIX01_REGION_COMPARISON_v01.csv')
blocks=read('BLOCK_MAGNITUDES_v01.csv')
assert len(regions)==162 and len(blocks)==4860
groups={}
for x in blocks:
    key=tuple(x[k] for k in ['candidate','branch','contrast','region','replica'])
    groups.setdefault(key,[]).append(x)
for r in regions:
    for i in [1,2,3]:
        key=tuple(r[k] for k in ['candidate','branch','contrast','region'])+(str(i),)
        g=groups[key]; assert sorted(int(x['block']) for x in g)==list(range(1,11))
        vals=np.array([float(x['norm']) for x in g])
        for suffix,expected in [('min',vals.min()),('max',vals.max()),('mean',vals.mean()),
                                ('median',np.median(vals)),('q10',np.quantile(vals,.1)),('q90',np.quantile(vals,.9))]:
            assert np.isclose(float(r[f'R{i}_block_norm_{suffix}']),expected,rtol=1e-12,atol=1e-14)
figs=read('FIGURE_PROVENANCE_v01.csv'); assert len(figs)==7
for f in figs:
    p=P/(f['figure']+'.png')
    with Image.open(p) as im: im.verify()
    pdf=(P/(f['figure']+'.pdf')).read_bytes()
    assert pdf.startswith(b'%PDF-') and b'%%EOF' in pdf[-1024:]
report=(P/'STAGE4_PBCFIX01_SCIENTIFIC_ANALYSIS_v01.md').read_text(encoding='utf-8')
assert version['script_sha256'] in report
assert all(x in report for x in ['PACERGEN02432','PACERGEN00123','PACERGEN00462','SIGNED_DRIFT','STATE_MOTION'])
review={'status':'PASS','input_hashes_rechecked':len(inputs),'all_input_hashes_unchanged':True,
        'block_distributions_replayed_from_delivered_csv':486,'png_pdf_pairs_verified':7,
        'source_code_hashes_match_version_record':True,'no_replica_exclusion':True,
        'tracked_git_diff':subprocess.check_output(['git','diff','--stat'],cwd=P.parents[2],text=True).strip(),
        'visual_review':'Heatmap, magnitude, specificity and region dependence layouts inspected as PNG; PDF saved from same figures',
        'report_sha256':sha(P/'STAGE4_PBCFIX01_SCIENTIFIC_ANALYSIS_v01.md'),
        'finalizer_sha256':sha(Path(__file__))}
dump('FINAL_DELIVERY_VALIDATION_v01.json',review)
files=[p for p in sorted(P.iterdir()) if p.is_file() and p.name!='OUTPUT_SHA256_v01.csv']
with (P/'OUTPUT_SHA256_v01.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['filename','bytes','sha256']); w.writeheader()
    w.writerows({'filename':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in files)
print(json.dumps(review,indent=2))
