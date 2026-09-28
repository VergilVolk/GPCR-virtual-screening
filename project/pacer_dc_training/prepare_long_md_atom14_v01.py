from pathlib import Path
import argparse, csv, hashlib, json, tempfile
import numpy as np
import MDAnalysis as mda

AA3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
"HIS":"H","HSE":"H","HSD":"H","HSP":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M",
"PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

ATOM14={
"A":"N CA C O CB".split(),"R":"N CA C O CB CG CD NE CZ NH1 NH2".split(),
"N":"N CA C O CB CG OD1 ND2".split(),"D":"N CA C O CB CG OD1 OD2".split(),
"C":"N CA C O CB SG".split(),"Q":"N CA C O CB CG CD OE1 NE2".split(),
"E":"N CA C O CB CG CD OE1 OE2".split(),"G":"N CA C O".split(),
"H":"N CA C O CB CG ND1 CD2 CE1 NE2".split(),"I":"N CA C O CB CG1 CG2 CD1".split(),
"L":"N CA C O CB CG CD1 CD2".split(),"K":"N CA C O CB CG CD CE NZ".split(),
"M":"N CA C O CB CG SD CE".split(),"F":"N CA C O CB CG CD1 CD2 CE1 CE2 CZ".split(),
"P":"N CA C O CB CG CD".split(),"S":"N CA C O CB OG".split(),
"T":"N CA C O CB OG1 CG2".split(),"W":"N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2".split(),
"Y":"N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH".split(),"V":"N CA C O CB CG1 CG2".split()
}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

p=argparse.ArgumentParser()
p.add_argument("--dcd",type=Path,required=True)
p.add_argument("--topology",type=Path,required=True)
p.add_argument("--graph",type=Path,required=True)
p.add_argument("--out",type=Path,required=True)
a=p.parse_args()

if a.out.exists():
    raise FileExistsError(f"REFUSE_OVERWRITE: {a.out}")
a.out.mkdir(parents=True)

with tempfile.TemporaryDirectory() as td:
    clean=Path(td)/"topology.pdb"
    with a.topology.open(errors="replace") as fi, clean.open("w") as fo:
        for line in fi:
            if not line.startswith("CONECT"): fo.write(line)

    u=mda.Universe(str(clean),str(a.dcd))
    rec=u.select_atoms("chainID E and protein")
    residues=list(rec.residues)

    seq="".join(AA3.get(r.resname,"?") for r in residues)
    graph=json.loads(a.graph.read_text(encoding="utf-8"))
    gseq="".join(n["site"][0] for n in graph["nodes"])

    if len(residues)!=270 or seq!=gseq:
        raise ValueError("receptor/graph ordering mismatch")

    idx=np.full((270,14),-1,dtype=np.int64)
    mapped=0

    for ri,(res,aa) in enumerate(zip(residues,seq)):
        actual={}
        for atom in res.atoms:
            name=atom.name
            if res.resname=="ILE" and name=="CD": name="CD1"
            actual.setdefault(name,[]).append(atom.index)

        for si,name in enumerate(ATOM14[aa]):
            hits=actual.get(name,[])
            if len(hits)!=1:
                raise ValueError(
                    f"expected exactly one heavy atom: residue={ri} "
                    f"{res.resname} atom={name} hits={len(hits)}"
                )
            idx[ri,si]=hits[0]
            mapped+=1

    n=len(u.trajectory)
    atom14=np.zeros((n,270,14,3),dtype=np.float32)
    max_ca_jump=0.0
    ca=rec.select_atoms("name CA")

    valid=idx>=0
    take=np.where(valid,idx,0)

    for fi,ts in enumerate(u.trajectory):
        if ts.dimensions is None or not np.all(np.isfinite(ts.dimensions)):
            raise ValueError(f"invalid box frame {fi}")
        pos=ts.positions
        if not np.isfinite(pos).all():
            raise ValueError(f"nonfinite coordinates frame {fi}")

        atom14[fi]=pos[take]*valid[...,None]

        d=np.linalg.norm(np.diff(ca.positions,axis=0),axis=1)
        max_ca_jump=max(max_ca_jump,float(d.max()))

    if max_ca_jump>=8.0:
        raise ValueError(f"PBC/backbone continuity failure: {max_ca_jump}")

    np.save(a.out/"apo_R1.atom14.npy",atom14)

    with (a.out/"apo_R1.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["name","seqres"]); w.writerow(["apo_R1",seq])

    audit={
        "status":"PASS",
        "frames":n,
        "residues":270,
        "mapped_heavy_atoms":mapped,
        "atom14_shape":list(atom14.shape),
        "dtype":str(atom14.dtype),
        "sequence_graph_match":seq==gseq,
        "max_ca_neighbor_distance_A":max_ca_jump,
        "dcd":str(a.dcd),
        "topology":str(a.topology),
        "topology_sha256":sha(a.topology),
        "atom14_sha256":sha(a.out/"apo_R1.atom14.npy"),
    }
    (a.out/"atom14_prep_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print("ATOM14_PREP_PASS")
    for k in ("frames","residues","mapped_heavy_atoms","atom14_shape","dtype",
              "sequence_graph_match","max_ca_neighbor_distance_A","atom14_sha256"):
        print(k, audit[k])
