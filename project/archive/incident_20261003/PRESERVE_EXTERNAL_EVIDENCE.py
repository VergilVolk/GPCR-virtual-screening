"""Evidence preservation only: enumerate, stream hash, and byte-copy existing files.
No project imports, scientific analysis, training, inference, docking or MD.
All writes are confined to this archive directory; sources are read-only.
"""
import csv, hashlib, json, re, shutil, subprocess, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
P = Path("C:/projects")
ROOTS = [
 P/"GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01",
 P/"GPCR-virtual-screening-final/project/results/generated",
 P/"GPCR-virtual-screening-final/project/results/pacer_candidates_v01",
 P/"GPCR-virtual-screening-final/project/data/generated",
 P/"GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01",
 P/"GPCR-virtual-screening/project/results/pacer_candidates_v01",
 P/"GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01",
 P/"GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01",
 P/"PACER_STAGE4_MD_backup",
 P/"GPCR-virtual-screening/project/results/pacer_stage4_prospective_fkg_v02_v01",
 P/"GPCR-virtual-screening/project/cache/pacer_stage4_prospective_fkg_v02_v01",
 P/"GPCR-virtual-screening-fkg-v02/project/results/pacer_fkg_v02_longmd_v01/calibration",
]
LIGHT = {".log",".txt",".json",".jsonl",".csv",".tsv",".yaml",".yml",".md",
 ".manifest",".receipt",".provenance",".out",".output",".err",".py",".sh",".ps1",
 ".pdbqt",".sdf",".smi",".toml"}
BINARY = {".dcd",".xtc",".trr",".nc",".chk",".ckpt",".checkpoint",".npy",".npz",
 ".pt",".pth",".safetensors",".lmdb",".lmdb-lock"}
def git(*args):
 return subprocess.check_output(["git","-c","core.quotepath=false",*args],cwd=REPO).decode("utf-8",errors="strict")
def digest(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
def writecsv(name, rows, fields):
 with (HERE/name).open("w",encoding="utf-8",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def normalize(raw, src):
 s=raw.replace("\\","/")
 s=re.sub(r"/+","/",s)
 if s.startswith("/mnt/c/"):s="C:/"+s[7:]
 if s.startswith("project/"):
  parts=src.as_posix().split("/project/",1)
  if len(parts)!=2:return None
  s=parts[0]+"/"+s
 if not re.match(r"^[A-Za-z]:/",s):return None
 p=Path(s)
 if not p.is_relative_to(P):return None
 return p
def walkrecords(obj):
 if isinstance(obj,dict):
  if isinstance(obj.get("path"),str):yield obj
  for v in obj.values():yield from walkrecords(v)
 elif isinstance(obj,list):
  for v in obj:yield from walkrecords(v)

files={p.resolve() for root in ROOTS if root.exists() for p in root.rglob("*") if p.is_file()}
files.add(P/"GPCR-virtual-screening/project/results/STAGE1_STAGE4_RESULTS_AUDIT_REPORT_v01.md")
# Traverse direct artifact records in current execution JSONs; historical freeze
# is followed once only so old calibration is not presented as a new run.
records=[];missing=[]
for src in sorted(files):
 if src.suffix.lower()!=".json":continue
 try: data=json.loads(src.read_text(encoding="utf-8-sig"))
 except (UnicodeError,json.JSONDecodeError):continue
 for rec in walkrecords(data):
  q=normalize(rec["path"],src)
  if q is None:continue
  if q.is_file():files.add(q.resolve());records.append({"record_file":str(src),"asset":str(q)})
  elif rec.get("sha256") and (rec.get("bytes") is not None):
   missing.append({"record_file":src.as_posix(),"recorded_path":rec["path"],
     "resolved_local_path":q.as_posix(),"size_bytes":rec["bytes"],"sha256":rec["sha256"],
     "status":"reference not locally available; not regenerated"})

prior=json.loads((HERE/"PRESERVATION_INPUT.json").read_text(encoding="utf-8"))
prior_by={x["path"]:x for x in prior["external_prior_inventory"]}
missing_prior=[x for x in prior_by if not Path(x).is_file()]
if missing_prior:raise RuntimeError("Missing prior evidence: "+str(missing_prior))
tracked=set(git("ls-files").splitlines())
# Reachable Git objects allow exact identification of NPZ/etc already preserved
# by the ancestry-only merges; these are recorded separately, never overwritten.
objects={}
for line in git("rev-list","--objects","HEAD").splitlines():
 if " " in line:
  oid,name=line.split(" ",1);objects.setdefault(name,set()).add(oid)
def repo_relative(path):
 parts=path.as_posix().split("/project/",1)
 return "project/"+parts[1] if len(parts)==2 else None
def git_blob_id(path,size):
 h=hashlib.sha1();h.update(("blob "+str(size)+"\0").encode())
 with path.open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
manifest=[];large=[];already=[];n=0;start=time.time()
for src in sorted(files):
 if not src.is_file():raise RuntimeError("Missing enumerated source "+str(src))
 before=src.stat();size=before.st_size;sha=digest(src)
 key=src.as_posix()
 if key in prior_by and sha!=prior_by[key]["sha256"]:
  raise RuntimeError("Prior evidence changed before archival: "+key)
 relative=repo_relative(src)
 exclusion=None
 if size>=50*1024*1024:exclusion="50 MiB archive safety cap; no new LFS policy"
 elif src.suffix.lower() in BINARY:exclusion="binary trajectory/checkpoint/model/cache/intermediate; indexed without copying"
 elif src.suffix.lower() not in LIGHT:exclusion="raw MD structure/state XML or runtime auxiliary asset; indexed without copying"
 if exclusion:
  oid=git_blob_id(src,size) if relative in objects else None
  if oid and oid in objects[relative]:
   already.append({"source_absolute_path":key,"git_relative_path":relative,"git_blob":oid,
    "size_bytes":size,"sha256":sha,"status":"exact bytes already in reachable history; no duplicate archive copy"})
  else:large.append({"absolute_path":key,"size_bytes":size,"sha256":sha,"reason_not_committed":exclusion})
 else:
  # Separate directory roots, original file names retained. No overwrites.
  rel=Path("project/archive/incident_20261003/external_evidence/C/projects")/src.relative_to(P)
  # Windows long-path prefix preserves full source layering and file names.
  dest=Path("\\\\?\\"+str((REPO/rel).resolve()))
  dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():
   if dest.stat().st_size!=size or digest(dest)!=sha:
    raise RuntimeError("Existing archive differs; refusing overwrite "+str(dest))
  else:shutil.copyfile(src,dest)
  if dest.stat().st_size!=size or digest(dest)!=sha:raise RuntimeError("Copy mismatch "+key)
  if "GPCR-virtual-screening-fkg-v02" in key:category="historical_frozen_dependency"
  elif "PACER_STAGE4_MD_backup" in key:category="stage4_md_execution_evidence"
  elif "pacer_stage4_prospective" in key:category="stage4_fkg_evidence"
  elif "stage3d" in key:category="stage3_selection_evidence"
  elif "drugclip_corrected" in key:category="stage3_drugclip_evidence"
  elif "m4_gamd_ensemble" in key:category="stage3_docking_pose_or_execution_evidence"
  elif "candidates_v01" in key:category="stage2_portfolio_or_historical_candidate_evidence"
  else:category="stage1_to_stage3_execution_or_dependency_evidence"
  manifest.append({"source_absolute_path":key,"archived_relative_path":rel.as_posix(),
   "size_bytes":size,"sha256":sha,"category":category})
 if src.stat().st_size!=size or src.stat().st_mtime_ns!=before.st_mtime_ns:
  raise RuntimeError("Source changed during read/copy "+key)
 n+=1
 if n%100==0:print(json.dumps({"processed":n,"total":len(files),"copied":len(manifest),
  "indexed":len(large),"elapsed_seconds":round(time.time()-start)}),flush=True)

writecsv("EXTERNAL_EVIDENCE_MANIFEST_SHA256.csv",manifest,
 ["source_absolute_path","archived_relative_path","size_bytes","sha256","category"])
writecsv("EXTERNAL_LARGE_ASSET_INDEX.csv",large,
 ["absolute_path","size_bytes","sha256","reason_not_committed"])
writecsv("ALREADY_PRESERVED_BINARY_ASSETS.csv",already,
 ["source_absolute_path","git_relative_path","git_blob","size_bytes","sha256","status"])
# Preserve recorded-but-unavailable reference identities separately; never
# substitute a derived file or fabricate hashes.
unique_missing={tuple(x.items()):x for x in missing}
preserved_by_sha={x["sha256"]:x for x in manifest+large}
for row in unique_missing.values():
 match=preserved_by_sha.get(row["sha256"])
 if match and int(match["size_bytes"])==int(row["size_bytes"]):
  row["matching_preserved_source"]=match.get("source_absolute_path",match.get("absolute_path"))
  row["archived_relative_path"]=match.get("archived_relative_path","")
  row["status"]="recorded relative location absent; exact bytes found in separately declared historical frozen root and copied/indexed"
 else:
  row["matching_preserved_source"]="";row["archived_relative_path"]=""
writecsv("RECORDED_NONLOCAL_ASSET_REFERENCES.csv",list(unique_missing.values()),
 ["record_file","recorded_path","resolved_local_path","size_bytes","sha256","status","matching_preserved_source","archived_relative_path"])
summary={"enumerated_files":len(files),"copied_files":len(manifest),"large_or_raw_files_not_in_git":len(large),
 "binary_files_already_in_reachable_history":len(already),"recorded_nonlocal_reference_count":len(unique_missing),
 "recorded_nonlocal_resolved_elsewhere_count":sum(bool(x["matching_preserved_source"]) for x in unique_missing.values()),
 "unresolved_recorded_reference_count":sum(not x["matching_preserved_source"] for x in unique_missing.values()),
 "prior_text_inventory_count":len(prior_by),"prior_text_all_preserved":all(
  p in {x["source_absolute_path"] for x in manifest} for p in prior_by),
 "copied_bytes":sum(x["size_bytes"] for x in manifest),"indexed_bytes":sum(x["size_bytes"] for x in large),
 "source_content_modified":False,"scientific_computation_executed":False,
 "no_lfs_policy_change":True,"archive_text_policy":"external_evidence/** -text; no normalization"}
(HERE/"PRESERVATION_COPY_RECEIPT.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
print(json.dumps(summary),flush=True)
