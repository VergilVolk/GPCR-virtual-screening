"""Read-only SFTP acceptance sample; credentials are prompted, never saved.

Requires paramiko in the invoking environment or --ssh-runtime. Downloads
the first COMPLETE job's first two committed segments and all 12 topologies.
"""
import argparse
import getpass
import hashlib
import json
import shlex
import sys
from pathlib import Path, PurePosixPath


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--host", default="36.103.199.153")
    p.add_argument("--user", default="ubuntu")
    p.add_argument("--ssh-runtime", type=Path)
    p.add_argument("--destination", type=Path, required=True)
    a = p.parse_args()
    if a.ssh_runtime:
        sys.path.insert(0, str(a.ssh_runtime.resolve()))
    import paramiko
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    # OpenSSH's verified host-key file is shared; never accept unknown hosts.
    known = Path.home() / ".ssh/known_hosts"
    if known.exists():
        client.load_host_keys(str(known))
    password = getpass.getpass("SSH password (not saved): ")
    client.connect(a.host, username=a.user, password=password, timeout=15,
                   allow_agent=False, look_for_keys=False)
    del password
    remote = r'''
import csv, hashlib, json
from pathlib import Path
root=Path('/data/pacer_stage4_rerun_v01_deployment/pacer_stage4_rerun_v01_pbcfix01')
def read(p): return json.loads(p.read_text())
def relative(p):
    p=Path(p)
    if p.is_absolute(): p=p.relative_to(root)
    assert '..' not in p.parts
    return p
rows=list(csv.DictReader((root/'manifests/STAGE4_RERUN_JOB_MATRIX.csv').open()))
inventory=[]
chosen=None
for row in rows:
    base=root/relative(row['output_path'])
    progress=read(base/'progress.json') if (base/'progress.json').is_file() else {}
    receipts=sorted(base.glob('segments/[0-9]*/receipt.json'))
    ends=[]
    for rp in receipts:
        r=read(rp)
        assert r['end_step']==int(rp.parent.name) and r['start_step']==r['end_step']-50000
        assert r['identity']==progress['identity']
        ends.append(r['end_step'])
    assert ends==list(range(50000,50000*(len(ends)+1),50000))
    inventory.append(dict(job_id=row['job_id'],status=progress.get('status'),completed_steps=progress.get('completed_steps'),committed_segments=len(ends)))
    if chosen is None and progress.get('status')=='COMPLETE' and len(ends)==100: chosen=row
assert chosen is not None, 'No COMPLETE 100-segment trajectory available'
files=set()
for line in (root/'SHA256SUMS_PBCFIX01.txt').read_text().splitlines(): files.add(relative(line.split('  ',1)[1]))
files.add(Path('SHA256SUMS_PBCFIX01.txt'))
files.update(map(Path,['authorization.json','evidence/source_snapshot.json','logs/server_preflight_pacer_stage4_rerun_v01_pbcfix01.json']))
for row in rows:
    files.update([Path('systems')/row['system_id']/'minimized.pdb',Path('systems')/row['system_id']/'build_audit.json'])
sid=chosen['system_id']
build=read(root/'systems'/sid/'build_audit.json')
files.update(Path('systems')/sid/name for name in build['outputs'])
files.add(Path('smoke')/sid/'audit.json')
prep=Path('equilibration')/sid/('replica_%02d'%int(chosen['replica']))
files.add(prep/'ready_pacer_stage4_rerun_v01_pbcfix01.json')
for label in ['nvt','npt','release_01','release_02','release_03','release_04']:
    files.add(prep/label/'receipt.json')
    files.update(prep/label/name for name in read(root/prep/label/'receipt.json')['outputs'])
base=relative(chosen['output_path'])
files.add(base/'progress.json')
for end in [50000,100000]:
    files.update(base/'segments'/('%08d'%end)/name for name in ['trajectory.dcd','state.csv','checkpoint.chk','receipt.json'])
records=[]
for rel in sorted(files):
    path=root/rel
    before=path.stat()
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    after=path.stat()
    assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns), 'Source changing: '+str(rel)
    records.append(dict(path=rel.as_posix(),bytes=after.st_size,sha256=h.hexdigest()))
print(json.dumps(dict(server_root=str(root),inventory=inventory,selected_job=chosen,files=records)))
'''
    try:
        _, stdout, stderr = client.exec_command("python3 -c " + shlex.quote(remote), timeout=120)
        payload = stdout.read()
        error = stderr.read().decode()
        if stdout.channel.recv_exit_status() != 0:
            raise RuntimeError(error)
        report = json.loads(payload)
        root = a.destination.resolve()
        root.mkdir(parents=True, exist_ok=True)
        with client.open_sftp() as sftp:
            for index, record in enumerate(report["files"], 1):
                rel = PurePosixPath(record["path"])
                if rel.is_absolute() or ".." in rel.parts:
                    raise ValueError("unsafe remote member")
                dest = (root / Path(*rel.parts)).resolve()
                if not dest.is_relative_to(root):
                    raise ValueError("destination path escape")
                dest.parent.mkdir(parents=True, exist_ok=True)
                def digest(path):
                    h = hashlib.sha256()
                    with path.open("rb") as f:
                        for b in iter(lambda: f.read(8*1024*1024), b""):
                            h.update(b)
                    return h.hexdigest()
                if dest.exists():
                    if dest.stat().st_size == record["bytes"] and digest(dest) == record["sha256"]:
                        continue
                    raise FileExistsError(f"existing download differs; preserve and use a fresh destination: {dest}")
                partial = dest.with_name(dest.name + ".download.partial")
                sftp.get(report["server_root"] + "/" + rel.as_posix(), str(partial))
                if partial.stat().st_size != record["bytes"] or digest(partial) != record["sha256"]:
                    raise RuntimeError(f"download SHA256 mismatch: {rel}")
                partial.rename(dest)
                print(f"download {index}/{len(report['files'])} {rel}", flush=True)
        report["local_root"] = str(root)
        report["status"] = "REAL_FILES_DOWNLOADED_SHA256_VERIFIED_NOT_YET_FKG_ACCEPTED"
        with (root / "REMOTE_SAMPLE_DOWNLOAD_v01.json").open("x", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        print(json.dumps({"status": report["status"], "selected_job": report["selected_job"]["job_id"],
                          "files": len(report["files"]), "bytes": sum(r["bytes"] for r in report["files"])}))
    finally:
        client.close()


if __name__ == "__main__":
    main()
