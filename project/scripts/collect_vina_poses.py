"""汇集 Vina pose + 受体 pdbqt 到共享目录（WSL 侧 GNINA 访问），再启动 GNINA 终链。"""
import glob, os, re, shutil, subprocess, tempfile
R = r'D:\CLC\project\results\pacer_rerun_v01'
VP = R + r'\vina_poses'
os.makedirs(VP, exist_ok=True)

# 受体三份
for n in ('7TRQ', '7TRP', '7TRS'):
    src = rf'D:\CLC\project\results\structure\ensemble\{n}_R_meeko.pdbqt'
    shutil.copy(src, rf'{VP}\{n}_R_meeko.pdbqt')

# Vina pose（各受体临时目录里的首个 MODEL = best pose）
pat = re.compile(r'MODEL 1.*?ENDMDL', re.S)
for name, tpat in (('7TRQ', 'vina7trq_*'), ('7TRP', 'vina7TRP_*'), ('7TRS', 'vina7TRS_*')):
    tmp = sorted(glob.glob(tempfile.gettempdir() + '\\' + tpat))
    if not tmp:
        print(name, '无临时 pose'); continue
    n_done = 0
    for f in glob.glob(tmp[-1] + '\\out*.pdbqt'):
        idx = int(re.search(r'out(\d+)', os.path.basename(f)).group(1))
        dst = rf'{VP}\{name}_pose{idx:04d}.pdbqt'
        if os.path.exists(dst):
            continue
        txt = open(f, errors='replace').read()
        m = pat.search(txt)
        if m:
            body = m.group(0).replace('MODEL 1', 'MODEL 1', 1) + '\n'
            open(dst, 'w').write('REMARK from Vina best pose\n' + body)
            n_done += 1
    print(name, 'pose 拷贝', n_done, '累计', len(glob.glob(rf'{VP}\{name}_pose*.pdbqt')))
