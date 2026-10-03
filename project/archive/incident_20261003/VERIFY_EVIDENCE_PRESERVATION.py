"""Verify preserved bytes and ancestry only; no scientific recomputation."""
import csv,hashlib,json,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
def extpath(p):return Path("\\\\?\\"+str(p.resolve()))
def sha(p):
 h=hashlib.sha256()
 with extpath(p).open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
def git(*args):return subprocess.check_output(["git","-c","core.quotepath=false",*args],cwd=REPO)
manifest=list(csv.DictReader((HERE/"EXTERNAL_EVIDENCE_MANIFEST_SHA256.csv").open(encoding="utf-8-sig",newline="")))
index={}
for raw in git("ls-files","-s","-z").split(b"\0"):
 if not raw:continue
 head,path=raw.split(b"\t",1)
 index[path.decode()]=head.decode().split()[1]
p=subprocess.Popen(["git","cat-file","--batch"],cwd=REPO,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
count=0
for row in manifest:
 rel=row["archived_relative_path"];source=Path(row["source_absolute_path"]);dest=REPO/rel
 expected=row["sha256"];size=int(row["size_bytes"])
 if sha(source)!=expected or sha(dest)!=expected:raise RuntimeError("Source/copy SHA mismatch: "+rel)
 if extpath(dest).stat().st_size!=size:raise RuntimeError("Copy size mismatch: "+rel)
 if rel not in index:raise RuntimeError("Evidence not staged: "+rel)
 p.stdin.write((index[rel]+"\n").encode());p.stdin.flush()
 header=p.stdout.readline().decode().strip().split()
 if len(header)!=3 or header[1]!="blob" or int(header[2])!=size:raise RuntimeError("Git blob size mismatch: "+rel)
 h=hashlib.sha256();remaining=size
 while remaining:
  b=p.stdout.read(min(8*1024*1024,remaining))
  if not b:raise RuntimeError("Unexpected blob EOF")
  h.update(b);remaining-=len(b)
 if p.stdout.read(1)!=b"\n":raise RuntimeError("Invalid batch boundary")
 if h.hexdigest()!=expected:raise RuntimeError("Staged Git blob modified bytes: "+rel)
 count+=1
p.stdin.close();p.wait()
modified=[]
baseline=json.loads((HERE/"TRACKED_WORKTREE_SHA256_BEFORE_ARCHIVE.json").read_text(encoding="utf-8"))
allowed={"project/docs/RECENT_EXECUTION_ARCHIVE_20261002_20261003.md","project/docs/INCIDENT_PACER_XR_PROTOCOL_DEVIATION_20261003.md"}
for row in baseline:
 if row["path"] in allowed:continue
 if sha(REPO/row["path"])!=row["sha256"]:modified.append(row["path"])
if modified:raise RuntimeError("Existing tracked file changed: "+str(modified))
commitrows=list(csv.DictReader((HERE/"COMMIT_REACHABILITY_AUDIT.csv").open(encoding="utf-8-sig",newline="")))
unreachable=[]
for row in commitrows:
 r=subprocess.run(["git","merge-base","--is-ancestor",row["commit_sha"],"HEAD"],cwd=REPO)
 if r.returncode:unreachable.append(row["commit_sha"])
if unreachable:raise RuntimeError("Unreachable original commit: "+str(unreachable))
same_index=(HERE/"TRACKED_INDEX_BEFORE_MERGES.txt").read_bytes()==(HERE/"TRACKED_INDEX_AFTER_MERGES.txt").read_bytes()
if not same_index:raise RuntimeError("Merge-stage index changed")
result={"status":"PASS","original_commit_count":len(commitrows),"reachable_from_head":len(commitrows),
 "copied_evidence_verified_against_source_worktree_and_staged_blob":count,
 "preexisting_tracked_files_verified_unchanged":len(baseline)-len(allowed),
 "allowed_document_updates":sorted(allowed),"scientific_result_files_modified":False,
 "ancestry_merge_tracked_index_identical":same_index,
 "no_scientific_jobs_or_tests_executed":True}
(HERE/"PRESERVATION_VERIFICATION_RECEIPT.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
print(json.dumps(result),flush=True)
