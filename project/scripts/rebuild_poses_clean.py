import glob, os, re, shutil, tempfile

R = r'D:\CLC\project\results\pacer_rerun_v01'
VP = R + r'\vina_poses'

def extract_pose(txt):
    m = re.search(r'MODEL 1(.*?)ENDMDL', txt, re.S)
    if not m:
        return None
    keep = []
    for ln in m.group(1).splitlines():
        s = ln.strip()
        if s.startswith(('REMARK', 'MODEL', 'ENDMDL')):
            continue
        keep.append(ln)
    return '\n'.join(keep) + '\n'

for name, tpat in (('7TRQ', 'vina7trq_*'), ('7TRP', 'vina7TRP_*'), ('7TRS', 'vina7TRS_*')):
    tmp = sorted(glob.glob(tempfile.gettempdir() + '\\' + tpat))
    if not tmp:
        print(name, 'no temp'); continue
    n = 0
    for f in glob.glob(tmp[-1] + '\\out*.pdbqt'):
        idx = int(re.search(r'out(\d+)', os.path.basename(f)).group(1))
        dst = rf'{VP}\{name}_pose{idx:04d}.pdbqt'
        body = extract_pose(open(f, errors='replace').read())
        if body:
            open(dst, 'w').write(body)
            n += 1
    print(name, 'rebuilt', n)
