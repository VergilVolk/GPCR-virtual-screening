"""从临时 pose 文件恢复 Vina 分数（无需重跑）。"""
import glob, os, re
import pandas as pd

R = r'D:\CLC\project\results\pacer_rerun_v01'
PAT = re.compile(r'REMARK VINA RESULT:\s+(-?\d+\.?\d*)')

def recover(name, tmp_pattern):
    tmp = sorted(glob.glob(os.environ['TEMP'] + '\\' + tmp_pattern))
    if not tmp:
        print(name, '无临时目录'); return
    tmp = tmp[-1]
    scores = {}
    for f in glob.glob(tmp + r'\out*.pdbqt'):
        m = PAT.search(open(f, errors='replace').read())
        if m:
            idx = int(os.path.basename(f)[3:7])   # out0007 -> 7
            scores[idx] = float(m.group(1))
    csv = f'{R}\\module3_{name}_results.csv'
    d = pd.read_csv(csv)
    col = f'vina_{name}'
    d[col] = [scores.get(i) for i in range(len(d))]
    d.to_csv(csv, index=False)
    print(f'{name}: 恢复 {d[col].notna().sum()}/200 分数（来自 {len(scores)} 个 pose）')

recover('7TRP', 'vina7TRP_*')
recover('7TRS', 'vina7TRS_*')
