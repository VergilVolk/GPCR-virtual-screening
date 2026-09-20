# -*- coding: utf-8 -*-
"""Export the final PACER hypotheses as annotated 2D SDF files."""
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
P=Path(__file__).resolve().parents[1];D=P/"results"/"pacer_candidates_v01"/"final"
def write(frame,path):
 w=Chem.SDWriter(str(path))
 for r in frame.itertuples(index=False):
  m=Chem.MolFromSmiles(r.canonical_smiles);AllChem.Compute2DCoords(m);m.SetProp("_Name",r.candidate_id)
  for k,v in r._asdict().items():
   if k!="canonical_smiles" and pd.notna(v):m.SetProp(str(k),str(v))
  m.SetProp("SMILES",r.canonical_smiles);w.write(m)
 w.close()
def main():
 d=pd.read_csv(D/"final_candidate_hypotheses.csv");write(d,D/"final_candidate_hypotheses.sdf");write(d[d.portfolio_role=="local_exploitation"],D/"local_exploitation_12.sdf");print(f"exported {len(d)} molecules")
if __name__=="__main__":main()
