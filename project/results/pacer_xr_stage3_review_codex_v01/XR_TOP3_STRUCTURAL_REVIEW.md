# Pure six-channel PACER-XR Top3 structural review

## Measured result

PACERGEN01495, PACERGEN01580 and PACERGEN01593 are pure-six-channel ranks 1, 2 and 3; they are Cascade ranks 10, 11 and 12 because Cascade first orders its nine Glide-BEmin top-1% gate members. This is retrospective score/pose evidence, not PAM activity validation.

|Candidate|Cascade|Six calibrated percentiles: G-PDB/G-BEmin/G-BEavg; V-PDB/V-BEmin/V-BEavg|G-V|E-PDB|BEmin receptors|PDB engine overlap|Evidence|
|---|--:|---|--:|--:|---|---|---|
|PACERGEN01495|10|0.712/0.986/0.981; 0.885/0.997/0.999|-0.068|0.192|G cluster_00; V cluster_00|shared=R:ASN423,R:ASP432,R:GLN184,R:ILE93,R:LEU190,R:PHE186,R:SER436,R:TRP435,R:TYR92; Jaccard=0.56|4/4 representative contexts available|
|PACERGEN01580|11|0.799/0.923/0.989; 0.861/0.989/0.994|-0.044|0.144|G cluster_00; V cluster_00|shared=R:ASN423,R:ASP432,R:GLN184,R:ILE93,R:LEU190,R:PHE186,R:SER436,R:TRP435,R:TYR89,R:TYR92; Jaccard=0.62|4/4 representative contexts available|
|PACERGEN01593|12|0.843/0.982/0.998; 0.726/0.997/0.999|0.034|0.210|G cluster_01; V cluster_00|shared=R:ASN423,R:ASP432,R:GLN184,R:ILE93,R:LEU190,R:PHE186,R:SER191,R:SER436,R:TRP435,R:TYR92; Jaccard=0.77|4/4 representative contexts available|

## PACERGEN01495
- Minimum channel percentile: 0.712. G/V means are 0.893/0.960; calibrated PDB/ensemble means are 0.798/0.991.
- PDB contacts: Glide R:ASN423 2.731A/7; R:GLN184 2.894A/7; R:TYR92 3.196A/4; R:SER191 3.215A/5; R:ASP432 3.303A/5; R:GLY96 3.361A/7; R:LEU190 3.496A/6; R:TRP435 3.527A/12; R:PHE186 3.536A/7; R:SER436 3.717A/3; R:GLN427 3.793A/1; R:ILE93 3.842A/4; R:THR179 3.936A/2. Vina R:SER436 3.148A/5; R:THR433 3.218A/3; R:ASP432 3.357A/6; R:TRP435 3.361A/7; R:ILE93 3.408A/11; R:TYR92 3.477A/8; R:LEU190 3.534A/2; R:PHE186 3.565A/16; R:ASN423 3.591A/2; R:ILE35 3.712A/1; R:GLN184 3.754A/2; R:TYR89 3.834A/2.
- BEmin contacts: Glide :ILE691 2.803A/6; :TYR875 3.128A/5; :SER872 3.321A/3; :PHE784 3.428A/19; :TYR687 3.497A/2; :TRP871 3.509A/23; :ASN859 3.515A/4; :GLY694 3.553A/7; :ILE692 3.609A/3; :MET630 3.646A/4; :TYR690 3.698A/2; :ASP868 3.823A/3; :LEU788 3.926A/1. Vina :SER872 2.597A/6; :ASN859 2.9A/8; :LEU788 2.968A/4; :ILE691 3.0A/14; :ASP868 3.164A/15; :TYR690 3.386A/3; :TRP871 3.431A/21; :PHE784 3.454A/27; :TYR875 3.47A/5; :MET630 3.474A/3; :GLY694 3.586A/2; :ILE633 3.713A/2. shared=:ASN859,:ASP868,:GLY694,:ILE691,:LEU788,:MET630,:PHE784,:SER872,:TRP871,:TYR690,:TYR875; Jaccard=0.79.
- PDB overlap supports a similar observed 7TRS contact region across engines. No PDB/ensemble residue mapping or hydrogen-bond assignment is made.

## PACERGEN01580
- Minimum channel percentile: 0.799. G/V means are 0.904/0.948; calibrated PDB/ensemble means are 0.830/0.974.
- PDB contacts: Glide R:THR424 2.785A/3; R:ASP432 2.826A/7; R:ASN423 2.952A/7; R:TRP435 3.164A/6; R:TYR92 3.177A/7; R:SER191 3.179A/8; R:PHE186 3.49A/20; R:ILE93 3.563A/3; R:GLN184 3.585A/1; R:GLN427 3.636A/3; R:SER436 3.814A/4; R:LEU190 3.824A/1; R:TYR89 3.878A/2; R:GLY96 3.879A/1. Vina R:SER436 3.063A/5; R:THR433 3.09A/4; R:ASP432 3.209A/11; R:TYR92 3.284A/8; R:TRP435 3.308A/6; R:ILE93 3.485A/11; R:PHE186 3.558A/16; R:ASN423 3.567A/2; R:LEU190 3.588A/3; R:TYR89 3.65A/2; R:ILE35 3.665A/1; R:GLN184 3.944A/2.
- BEmin contacts: Glide :SER872 3.059A/4; :ASN859 3.173A/7; :THR869 3.315A/3; :LEU788 3.398A/2; :TYR690 3.403A/9; :PHE784 3.422A/27; :TYR687 3.452A/1; :TRP871 3.459A/15; :GLY694 3.62A/1; :MET630 3.632A/1; :ILE691 3.735A/6. Vina :TYR687 2.819A/6; :ASN859 3.069A/17; :CYX862 3.42A/2; :LEU788 3.443A/2; :ILE691 3.504A/3; :SER789 3.522A/5; :TYR690 3.548A/19; :TRP871 3.55A/10; :PHE784 3.566A/22; :SER864 3.619A/1; :SER872 3.777A/1; :THR860 3.847A/1. shared=:ASN859,:ILE691,:LEU788,:PHE784,:SER872,:TRP871,:TYR687,:TYR690; Jaccard=0.53.
- PDB overlap supports a similar observed 7TRS contact region across engines. No PDB/ensemble residue mapping or hydrogen-bond assignment is made.

## PACERGEN01593
- Minimum channel percentile: 0.726. G/V means are 0.941/0.907; calibrated PDB/ensemble means are 0.784/0.994.
- PDB contacts: Glide R:ASN423 2.848A/11; R:PHE186 3.112A/6; R:TRP435 3.199A/17; R:SER191 3.281A/11; R:ILE93 3.506A/3; R:LEU190 3.575A/6; R:ASP432 3.636A/2; R:GLN427 3.714A/4; R:GLY96 3.737A/1; R:SER436 3.771A/2; R:TYR92 3.868A/1; R:GLN184 3.934A/1. Vina R:TYR92 2.976A/6; R:ILE93 3.092A/8; R:GLN184 3.359A/4; R:ASN423 3.468A/3; R:ASP432 3.47A/6; R:THR433 3.485A/2; R:PHE186 3.558A/15; R:TRP435 3.627A/4; R:LEU190 3.672A/3; R:SER436 3.764A/2; R:SER191 3.869A/2.
- BEmin contacts: Glide :TRP876 3.359A/12; :SER872 3.431A/14; :ILE633 3.535A/1; :TRP871 3.602A/1; :ILE873 3.62A/5; :LEU640 3.67A/2; :TYR687 3.672A/1; :THR637 3.72A/1; :PHE784 3.721A/6; :LEU877 3.729A/3; :TYR690 3.755A/1; :ILE691 3.815A/2; :VAL636 3.824A/1. Vina :ASN859 2.909A/16; :TYR687 3.345A/5; :ILE691 3.356A/3; :LEU788 3.379A/4; :TYR690 3.481A/15; :ASP868 3.605A/2; :PHE784 3.617A/13; :SER789 3.688A/5; :TRP871 3.689A/15; :CYX862 3.835A/1; :GLN863 3.907A/1. not compared: different winning receptor conformations.
- PDB overlap supports a similar observed 7TRS contact region across engines. No PDB/ensemble residue mapping or hydrogen-bond assignment is made.

## Comparison with Cascade Top10

- The Top3 minima are 0.712, 0.799, and 0.726; their G-V differences are -0.068, -0.044, and +0.034; their E-PDB differences are 0.192, 0.144, and 0.210. Each has all four requested representative contact contexts.
- Four Cascade Top10 gate cases have pure-XR ranks 66-122 and E-PDB 0.564-0.684. The Top3 therefore have more balanced six-channel evidence under the existing transforms, without establishing which strategy predicts activity.
- BEmin winning receptor agrees across engines for 01495 and 01580 (cluster_00). It differs for 01593 (Glide cluster_01, Vina cluster_00), so ensemble cross-engine residue comparison is intentionally withheld for 01593.
