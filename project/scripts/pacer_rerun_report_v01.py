#!/usr/bin/env python3
"""Phase E, figures/report/validation from saved outputs only. Does not rerank.
Usage: MPLCONFIGDIR=/tmp/pacer_mpl .venv-pacer-admet-v01/bin/python project/scripts/pacer_rerun_report_v01.py [--output DIR]
"""
import argparse,json, math, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from rdkit import Chem
from rdkit.Chem import Draw
from pacer_rerun_medchem_v01 import ROOT,INPUT,DEFAULT,TOP3,KEY,sha,csv,dump,identity,encode

BLUE='#3477A5';GREEN='#268578';ORANGE='#D88A2B';GRAY='#B6BCC3';RED='#BB5864'
def markdown(d):return d.to_markdown(index=False,floatfmt='.3f')
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=DEFAULT);o=ap.parse_args().output
    read=lambda n:pd.read_csv(o/(n+'_v01.csv'))
    phy=read('ALL200_PHYSICOCHEMICAL');drug=read('ALL200_DRUGLIKENESS');st=read('ALL200_STEREOCHEMISTRY');al=read('ALL200_STRUCTURAL_ALERTS');rea=read('ALL200_REACTIVITY_ALERTS');syn=read('ALL200_SYNTHESIS_FEASIBILITY');retro=read('ALL200_RETROSYNTHESIS');cn=read('ALL200_CNS_COMPONENTS');a=read('ALL200_ADMET');t=read('ALL200_TOXICITY');allpred=pd.concat([a,t],ignore_index=True)
    metadata=json.loads((o/'MODEL_PROVENANCE_v01.json').read_text());iv=json.loads((o/'INPUT_VALIDATION_v01.json').read_text())
    if not (o/'MODEL_ENDPOINT_METADATA_v01.json').exists():dump(o,'MODEL_ENDPOINT_METADATA_v01.json',metadata['endpoint_metadata'])
    # Endpoint table retains raw outputs, usable values, status and model SD independently.
    top=phy[phy.candidate_id.isin(TOP3)].set_index('candidate_id').loc[TOP3].reset_index()
    for table in [drug,st,al,syn,cn,rea]:
        new=[c for c in table.columns if c not in top.columns or c in KEY];top=top.merge(table[new],on=KEY,validate='one_to_one')
    for ep,g in allpred.groupby('endpoint'):
        for col in ['prediction','raw_prediction','ensemble_sample_std','status']:
            if col in g:top[ep+'__'+col]=top.candidate_id.map(g.set_index('candidate_id')[col])
    sol=allpred[allpred.endpoint=='Solubility_AqSolDB'].merge(phy[KEY+['molecular_weight']],on=KEY,validate='one_to_one')
    sol['solubility_mol_per_L']=10**sol.prediction;sol['solubility_mg_per_L']=sol.solubility_mol_per_L*sol.molecular_weight*1000
    sol['conversion_note']='Algebraic transformation of model log10(mol/L), not independent measurement; pH, crystal form, temperature and salt not controlled.'
    csv(o,'ALL200_WATER_SOLUBILITY_v01.csv',sol)
    for col in ['solubility_mol_per_L','solubility_mg_per_L']:top[col]=top.candidate_id.map(sol.set_index('candidate_id')[col])
    top['A_physchem_summary']=['MW 516.034; MolLogP 5.262; TPSA 103.88; 6 aromatic rings; lower saturation.','MW 480.041; MolLogP 3.118; TPSA 74.72; FractionCSP3 0.417.','MW 494.068; MolLogP 3.426; TPSA 74.72; FractionCSP3 0.440.']
    top['B_druglikeness_summary']=['Lipinski 2 violations (MW/logP), Ghose 2, Muegge 1; Veber/Egan pass.','Lipinski/Veber/Egan/Muegge pass, Ghose 2 violations (MW just above 480 and MR).','Lipinski/Veber/Egan/Muegge pass, Ghose 2 violations (MW/MR).']
    top['C_ADMET_summary']='Exploratory 5-member ADMET-AI v2 outputs; per-endpoint status/SD retained; AD unknown. Solubility and permeability require assay confirmation.'
    top['D_CNS_BBB_risk']='BBB classification only; no brain exposure, Kp,uu, Pgp substrate, full CNS MPO, logD or pKa evidence.'
    top['E_toxicity_risk']='hERG, DILI and CYP3A4 inhibition model probabilities warrant experimental confirmation; no safe label.'
    top['F_synthesis_summary']='SA heuristic only. Route NOT_ASSESSED; no reaction step count, verified starting materials, yield or cost.'
    top['G_stereo_salt_summary']='No RDKit potential tetrahedral/double-bond stereocenters found; neutral single component. Ionizable motifs present, no pKa/speciation or best-salt evidence.'
    top['H_evidence_gaps']='Unknown AD and uncalibrated probabilities; no experimental ADMET; Caco2 units unresolved; microsomal/systemic PK gap; retrosynthesis asset/stock unavailable.'
    csv(o,'TOP3_COMPREHENSIVE_PROFILE_v01.csv',top)
    (o/'TOP3_STRUCTURES_v01.png').exists() and (_ for _ in ()).throw(FileExistsError('Top3 image exists'))
    Draw.MolsToGridImage([Chem.MolFromSmiles(top.set_index('candidate_id').loc[x,'canonical_smiles']) for x in TOP3],molsPerRow=3,subImgSize=(500,360),legends=TOP3).save(str(o/'TOP3_STRUCTURES_v01.png'))
    evidence=[]
    for r in phy.to_dict('records'):
        k=identity(r);cid=r['candidate_id']
        for ep,typ in [('physicochemical','COMPUTED_DESCRIPTOR'),('druglikeness_QED','EMPIRICAL_RULE_HEURISTIC'),('stereochemistry','COMPUTED_STRUCTURE'),('SA_score','COMPUTED_HEURISTIC')]:evidence.append({**k,'endpoint':ep,'status':'COMPUTED','evidence_grade':'A_REPRODUCIBLE_CALCULATION' if typ in ['COMPUTED_DESCRIPTOR','COMPUTED_STRUCTURE'] else 'B_HEURISTIC','evidence_type':typ,'reason':None})
        for ep,tab,col in [('PAINS',al,'PAINS_hit'),('Brenk',al,'Brenk_hit'),('bioactivation_motifs',rea,'prediction')]:
            hit=bool(tab.set_index('candidate_id').loc[cid,col]);evidence.append({**k,'endpoint':ep,'status':'STRUCTURAL_ALERT' if hit else 'SCREENED_NO_MATCH','evidence_grade':'B_STRUCTURE_SCREEN','evidence_type':'STRUCTURAL_ALERT','reason':'No-match does not establish safety'})
        for ep,status in [('retrosynthesis','NOT_ASSESSED'),('starting_material_availability','AVAILABILITY_NOT_VERIFIED'),('reaction_steps','NOT_ASSESSED'),('yield','NOT_ESTIMATED'),('cost','NOT_ESTIMATED')]:evidence.append({**k,'endpoint':ep,'status':status,'evidence_grade':'D_MISSING','evidence_type':'MISSING','reason':'No complete tool assets/routes or actual procurement/experimental evidence'})
    ev=pd.concat([pd.DataFrame(evidence),allpred[KEY+['endpoint','status','evidence_grade','status_reason']].rename(columns={'status_reason':'reason'}).assign(evidence_type='MODEL_OR_MISSING')],ignore_index=True);csv(o,'EVIDENCE_COMPLETENESS_v01.csv',ev)
    summary={ep:{'status_counts':g.status.value_counts().to_dict(),'min_prediction':float(g.prediction.min()) if g.prediction.notna().any() else None,'median_prediction':float(g.prediction.median()) if g.prediction.notna().any() else None,'max_prediction':float(g.prediction.max()) if g.prediction.notna().any() else None} for ep,g in allpred.groupby('endpoint')};dump(o,'ENDPOINT_SUMMARY_v01.json',summary)
    figdir=o/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':300})
    def save(fig,n,title,foot):
        fig.suptitle(title,fontsize=15,fontweight='bold',y=.99);fig.text(.02,.015,foot,fontsize=9,color='#444444');fig.tight_layout(rect=[0,.07,1,.94])
        for ext in ['png','pdf']:
            p=figdir/f'Figure_{n}_v01.{ext}';assert not p.exists();fig.savefig(p,bbox_inches='tight',facecolor='white')
        plt.close(fig)
    fig,axs=plt.subplots(2,3,figsize=(12,7.6))
    for ax,col,label in zip(axs.flat,['molecular_weight','MolLogP','TPSA','HBA','rotatable_bonds','FractionCSP3'],['MW (g/mol)','MolLogP','TPSA (Å²)','HBA','Rotatable bonds','Fraction Csp3']):
        ax.hist(phy[col],bins=16,color=BLUE,alpha=.85,edgecolor='white');ax.set_xlabel(label);ax.set_ylabel('Candidates')
    save(fig,1,'PACERGEN-200 | Physicochemical distributions','Blue: calculated RDKit descriptors (n=200). These are not ADMET experiments or activity measurements.')
    fig,axs=plt.subplots(1,2,figsize=(12,4.9));axs[0].hist(drug.QED,bins=18,color=BLUE,edgecolor='white');axs[0].set(xlabel='QED (weighted RDKit default)',ylabel='Candidates')
    rules=['Lipinski','Veber','Ghose','Egan','Muegge'];passed=[int((drug[x+'_status']=='PASS').sum()) for x in rules];failed=[200-x for x in passed]
    axs[1].barh(rules,passed,color=BLUE,label='Zero violations');axs[1].barh(rules,failed,left=passed,color=ORANGE,label='One or more violations');axs[1].set_xlabel('Candidates');axs[1].legend(loc='lower right',fontsize=8)
    for i,x in enumerate(passed):axs[1].text(2,i,str(x),va='center',color='white',fontsize=9)
    save(fig,2,'PACERGEN-200 | QED and empirical rule checks','Calculated heuristics; PASS means zero violations. Muegge uses MolLogP proxy; Egan uses common rectangular limits. No elimination.')
    eps=['Solubility_AqSolDB','BBB_Martins','Caco2_Wang','PPBR_AZ','Clearance_Microsome_AZ','hERG'];labels=['logS [log10(mol/L)]','BBB class-1 probability','Caco2 raw model scale*','PPBR (%)','Microsome CLint (µL/min/mg)','hERG blocking class-1 probability']
    fig,axs=plt.subplots(2,3,figsize=(12,8))
    for ax,ep,label in zip(axs.flat,eps,labels):
        g=allpred[allpred.endpoint==ep];v=g[g.status=='PREDICTED'].prediction;ax.hist(v,bins=17,color=GREEN,edgecolor='white');ax.set(xlabel=label,ylabel='Candidates');ax.set_title(f'{len(v)}/200 usable; AD UNKNOWN',fontsize=10)
        if ep=='PPBR_AZ':ax.text(.04,.85,'34 outputs >100%\nMODEL_FAILED',transform=ax.transAxes,color=RED,fontsize=9)
    save(fig,3,'PACERGEN-200 | Exploratory ADMET model outputs','Green: model predictions, no established per-molecule AD. *Caco2 units conflict; no conversion. Invalid outputs excluded and counted.')
    fig,axs=plt.subplots(2,3,figsize=(12,8));short=['02432','00123','00462']
    for ax,ep,label in zip(axs.flat,['Solubility_AqSolDB','BBB_Martins','Pgp_Broccatelli','hERG','AMES','DILI'],['logS [log10(mol/L)]','BBB class-1 probability','Pgp inhibitor class-1 probability','hERG blocker class-1 probability','Ames class-1 probability','DILI class-1 probability']):
        g=allpred[allpred.endpoint==ep].set_index('candidate_id').loc[TOP3];ax.errorbar(range(3),g.prediction,yerr=g.ensemble_sample_std,fmt='o',color=GREEN,capsize=5,markersize=7);ax.set_xticks(range(3),short);ax.set_ylabel(label);ax.set_title('AD UNKNOWN | model only',fontsize=9)
        if ep!='Solubility_AqSolDB':ax.set_ylim(-.05,1.22)
    save(fig,4,'Frozen Top3 | ADMET prediction comparison','Points: five-member means. Bars: sample SD (not confidence intervals). Inhibitor ≠ substrate; no calibrated risk threshold applied.')
    fig,axs=plt.subplots(1,2,figsize=(12,5.2));cats=['NOT_ASSESSED','NO_ROUTE_FOUND','ROUTE_FOUND'];counts=[int((retro.status==x).sum()) for x in cats];axs[0].bar(cats,counts,color=[GRAY,ORANGE,GREEN]);axs[0].tick_params(axis='x',labelrotation=12);axs[0].set_ylabel('Candidates')
    for i,v in enumerate(counts):axs[0].text(i,v+3,str(v),ha='center')
    axs[0].set_ylim(0,225);axs[1].hist(syn.SA_score,bins=15,color=BLUE,edgecolor='white');axs[1].set(xlabel='SA score (1 easier → 10 harder)',ylabel='Candidates');axs[1].set_title('Calculated complexity heuristic, not a route')
    save(fig,5,'PACERGEN-200 | Synthesis evidence and missing routes','Gray: unassessed (200). Policy/template/stock gate failed; searches executed=0. SA does not predict steps, yield or price.')
    matrix=[];mlabel=['Properties','Rules / QED','PAINS screen','Bioactivation\nmotifs','logS / BBB','hERG / DILI','PPBR','Pgp substrate','Routes','Costs']
    for cid in TOP3:matrix.append([0,0,2,2,1,1,1 if allpred[(allpred.candidate_id==cid)&(allpred.endpoint=='PPBR_AZ')].status.iloc[0]=='PREDICTED' else 4,3,3,3])
    fig,ax=plt.subplots(figsize=(13,4.6));ax.imshow(matrix,cmap=ListedColormap([BLUE,GREEN,ORANGE,GRAY,RED]),vmin=-.5,vmax=4.5,aspect='auto');ax.set_xticks(range(len(mlabel)),mlabel);ax.set_yticks(range(3),TOP3);ax.tick_params(axis='x',labelsize=9)
    for i in range(3):
        for j,v in enumerate(matrix[i]):ax.text(j,i,['Computed','Predicted','No match' if j==2 else 'Alert','Unassessed','Failed'][v],ha='center',va='center',fontsize=8,color='white' if v in [0,1,4] else '#333333')
    ax.legend(handles=[Patch(color=c,label=l) for c,l in [(BLUE,'Descriptor / heuristic'),(GREEN,'Prediction; AD unknown'),(ORANGE,'Structure screen'),(GRAY,'No assessment'),(RED,'Invalid model output')]],loc='upper center',bbox_to_anchor=(.5,-.15),ncol=3,fontsize=8)
    save(fig,6,'Frozen Top3 | Evidence types, coverage and gaps','Cell colors denote evidence type or missingness, not an overall drug-development score. No-match PAINS does not mean safe.')
    # Compact narrative tables from actual saved values.
    rule_table=pd.DataFrame({'规则':rules,'零违反PASS':passed,'FAIL':failed})
    ranges=phy[['molecular_weight','MolLogP','TPSA','rotatable_bonds','FractionCSP3']].agg(['min','median','max']).T.reset_index(names='指标')
    topvalues=top[['candidate_id','molecular_weight','MolLogP','TPSA','FractionCSP3','QED','Solubility_AqSolDB__prediction','BBB_Martins__prediction','hERG__prediction','AMES__prediction','DILI__prediction','SA_score']]
    counts=allpred.status.value_counts().to_dict();allpass=int((drug[[x+'_status' for x in rules]]=='PASS').all(axis=1).sum())
    def ep_summary(ep):
        g=allpred[allpred.endpoint==ep];valid=g[g.status=='PREDICTED'];return f"{len(valid)}/200 PREDICTED；"+(f"中位数 {valid.prediction.median():.3f}，范围 {valid.prediction.min():.3f}–{valid.prediction.max():.3f}" if len(valid) else '无可采用数值')
    report=f'''# PACER-200 药化性质、ADMET与合成可行性综合评价 v01

分析日期：2026-10-12（Asia/Shanghai）。这是独立事后评价，不是2026年10月原始筛选门槛。冻结Top3、PACER-XR、DrugCLIP、Stage4和PACER-FKG均不因本报告改变；没有执行分子生成、Docking或MD。全部数值来自本次保存的实际计算或推理，无本次实验数据。

## 1. 输入与环境（PHASE A）

唯一分子输入：`{INPUT.relative_to(ROOT)}`。SHA256：`{iv['sha256']}`。起始Git提交：`{iv['git_commit']}`；起始结果分支：`{iv['parent_result_branch']}`；分析分支：`{iv['git_branch']}`。初始工作区干净。

200行、200个唯一candidate_id、200个唯一原始canonical SMILES及200个RDKit规范异构SMILES；全部成功解析和价态消毒检查，candidate_id与SMILES一一对应；原始SMILES列与canonical列一致，规范化序列化差异0。Top3均从此文件读取；未使用或映射历史PACER编号。没有身份冲突。

隔离环境`.venv-pacer-admet-v01`采用Python 3.11.9，RDKit {metadata['RDKit_version']}；完整依赖见两个lock文件与ENVIRONMENT。原默认Python 3.14.7缺少RDKit/pandas/matplotlib；未改变项目既有环境。仓库检索未发现可直接复用的专用ADMET或逆合成权重/库存；原药化脚本涉及其他分析背景，因此本批独立实现，不借用历史分子。

## 2. 全200分子基础性质与类药性（PHASE B）

{markdown(ranges)}

{markdown(rule_table)}

五条规则全部零违反：{allpass}/200。Lipinski严格PASS指0违反，另保存“至多1违反”字段。规则每个具体阈值、违反条件与数量均保留；QED是原始连续分数，未赋予淘汰阈值或综合权重。

QED中位数{drug.QED.median():.4f}，均值{drug.QED.mean():.4f}，范围{drug.QED.min():.4f}–{drug.QED.max():.4f}。PAINS命中{int(al.PAINS_hit.sum())}/200；Brenk命中{int(al.Brenk_hit.sum())}/200，单独保存规则名称、catalog标识、输入分子0起始原子索引和片段SMILES。PAINS零命中不证明安全。

实现约定：RDKit MolLogP一致用于规则；Muegge采用MolLogP替代原始/其他实现中的脂溶性算法，Egan采用常见矩形阈值近似，不能称为原版判别器完全复现；Ghose原子数包括隐式H。QED为RDKit默认加权几何平均。规则是经验标准，通过不证明ADMET优良，违反也不等于不能成药。[规则实现背景](https://pmc.ncbi.nlm.nih.gov/articles/PMC5335600/)

## 3. ADMET模型、证据等级与不确定性（PHASE C）

实际运行ADMET-AI **2.0.1** / Chemprop {metadata['chemprop']}，分类及回归各5个本地权重，CPU推理覆盖200×41=8200个原始端点输出，逐成员预测均保存。权重、元数据、SA文件SHA256及端点训练来源、单位、报告误差全部记录于MODEL_PROVENANCE。包为MIT许可证；训练来源为TDC，底层数据集许可证各异。2024论文描述v1，v2重新训练后输出不同；本报告不将v1验证指标自动视为v2或PACERGEN验证。[官方版本说明](https://github.com/swansonk14/admet_ai)

证据等级：A为可复现描述符/结构计算（不是实验验证）；B为经验规则/结构警示；C为可追溯模型的探索性预测，适用域未验证；D为缺失、失效或不可靠。所有模型的逐分子适用域为**UNKNOWN**：权重包未提供端点训练结构或经验证的适用域估计器，不能声称IN_DOMAIN或任意设置OUT_OF_DOMAIN阈值。五成员样本标准差仅表示集成分歧，不能量化全部模型偏差，不能作为置信区间或保证低误差。分类保留P(class1)，没有应用安全阈值。

在52个记录端点（41个模型+11个明确缺项）中，每个分子都有记录；总体状态计数为 `{counts}`。PREDICTED仅表示数值推理成功且通过物理边界检查，不表示模型在本批化学空间已被验证。无实验参考数据被当成本次实验，也未计算综合成药概率。

### 水溶性

{ep_summary('Solubility_AqSolDB')}，单位log10(mol/L)。使用AqSolDB训练的回归模型；已另外保存代数转换的mol/L及mg/L，明确不是新实验。pH、盐型、晶型和温度未指定，不能视为某固定条件下的本征溶解度。[训练端点](https://tdcommons.ai/single_pred_tasks/adme/#solubility-aqsoldb)

### BBB、CNS与通透性

BBB：{ep_summary('BBB_Martins')}。这是穿透分类概率，不是logBB、脑组织浓度或Kp,uu。实际脑暴露200/200 NOT_ASSESSED。CNS MPO仅保存MW、MolLogP代理、TPSA、HBD四个单项分量及明确标记的0–4部分和；缺乏pKa与pH7.4 logD，完整0–6 CNS MPO未计算。[CNS MPO原始方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC3368654/)

Caco2：{ep_summary('Caco2_Wang')}。模型包声明`log(10^-6 cm/s)`，TDC榜单标注`cm/s`，两者不一致；本次保留原始模型尺度并明确标记单位冲突，**没有换算绝对通透率或作高/低通透判定**。[TDC端点单位](https://tdcommons.ai/benchmark/admet_group/01caco2/)

PAMPA为分类概率，200/200完成；它是人工膜体系，不与Caco2或BBB等价。P-gp抑制：{ep_summary('Pgp_Broccatelli')}。P-gp底物200/200 NOT_ASSESSED；不能从抑制概率推出外排底物性，CNS外排风险仍是关键证据缺口。

### CYP、微粒体与蛋白结合

五个CYP抑制端点及CYP2C9/2D6/3A4底物均完成200/200预测；底物与抑制端点独立保存。代谢位点200/200 NOT_ASSESSED，无原子位点或确定代谢途径结论。

人微粒体CLint：{ep_summary('Clearance_Microsome_AZ')}，µL/min/mg；随包报告R²≈0.277、MAE≈26.74，精度有限，适用域未知。肝细胞CLint也完成200/200，保留其原始µL/min/10^6 cells单位及混合物种来源说明。微粒体半衰期与稳定性分类均未评价，没有进行体内外缩放。

PPBR：{ep_summary('PPBR_AZ')}。另34个输出超过100%，保留raw_prediction并将prediction置空，状态MODEL_FAILED；不截断为100%。Fraction unbound没有独立可靠模型，未伪装为另一个预测结果。

### 清除率、半衰期与口服生物利用度

系统清除率0/200可靠预测；体外CLint没有被称作全身清除率。Half_Life_Obach虽实际推理200次，但随包元数据R²≈−2.386、MAE≈31.85 h，缺乏可靠支持，全部NOT_ASSESSED，原始模型值仅归档。VDss模型同样因负R²未采纳。Bioavailability_Ma完成200个分类概率，不能解释为口服F百分比；连续F为NOT_ASSESSED。未拼接独立模型计算任何PK参数。

### 毒性及反应性代谢物

hERG：{ep_summary('hERG')}；Ames：{ep_summary('AMES')}；DILI：{ep_summary('DILI')}。另保留ClinTox、Tox21等模型输出，不能将受体/应激端点活性等同于已确认系统毒性。DILI是特定数据集端点，未据此生成独立泛肝毒性结论。任何低概率都不是“安全”结论。

非穷尽性结构匹配发现{int((rea.status=='STRUCTURAL_ALERT').sum())}/200命中反应性或生物活化相关motif（具体规则、SMARTS、原子及片段单独保存）。所有此类结果均标为STRUCTURAL_ALERT；包含噻吩、呋喃、硝基芳环、酚及少数反应性官能团，不能推断这些分子必然形成反应性代谢物。环融合、取代位置、竞争代谢与解毒都会改变意义。专门反应性代谢物模型未运行。[结构警示研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC5871347/)

## 4. 逆合成、采购、收率与成本（PHASE D）

在另一隔离环境安装AiZynthFinder 4.4.1，修复PyArrow 26与NumPy 1.26不兼容问题（固定PyArrow 18.1.0），完成实际导入及对象初始化。正式搜索前检查策略模型、模板、库存及许可证。公开模型Zenodo记录为CC-BY-4.0；代码MIT。Zenodo首次预检超时，随后沙箱外复核模型URL返回HTTP200；官方ZINC库存标准GET仍返回403，本地未找到所需完整资产，因此没有继续下载策略权重并声称可搜索。没有绕过限制，也没有把安装包当作完成逆合成。[官方工具与资产](https://github.com/MolecularAI/aizynthfinder)

因此200个分子均为NOT_ASSESSED，真正运行路线搜索0个、获得路线0个；**NO_ROUTE_FOUND为0个**。这里的0路线是未评价结果，不能说明分子不可合成。没有路线中间体、起始原料、步骤或路线评分可提供。SA Score全部计算，只用于相对结构复杂度描述，不能替代路线可行性、实验收率或成本。

起始原料现实采购核实0个，AVAILABILITY_NOT_VERIFIED；无供应商、CAS、产品链接、库存或报价被虚构。反应步数0个可靠证据，收率和总成本均200/200 NOT_ESTIMATED；币种、日期、纯化/人工/耗材及收率考虑字段留空，因为尚无路线和真实报价基础。保护/脱保护、关键转化和特殊条件同样不作无证据推断。

## 5. 立体化学与盐型

RDKit识别的已指定四面体中心分子{int((st.specified_tetrahedral_centers>0).sum())}个、未指定四面体中心分子{int((st.unspecified_tetrahedral_centers>0).sum())}个、未指定双键立体分子{int((st.unspecified_double_bond_stereo>0).sum())}个。所有分子为中性单组分。按保留已指定构型、仅枚举未指定构型的方法，理论异构体计数与枚举上限均已保存。零RDKit中心不排除阻转异构、构象或实验稳定立体形式。

结构存在可能可电离的胺/芳香氮等motif；匹配不意味着这些位点都在生理pH质子化。未运行pKa或pH依赖质子化预测，未确定最佳盐型，原始SMILES从未被改变。

## 6. 冻结Top3完整比较

下列列名中的预测值均为模型原尺度；所有概率适用域未知。

{markdown(topvalues)}

'''
    for cid in TOP3:
        r=top.set_index('candidate_id').loc[cid];g=allpred[allpred.candidate_id==cid];maincols=g[g.endpoint.isin(['Solubility_AqSolDB','BBB_Martins','Caco2_Wang','Pgp_Broccatelli','CYP1A2_Veith','CYP2C9_Veith','CYP2C19_Veith','CYP2D6_Veith','CYP3A4_Veith','Clearance_Microsome_AZ','PPBR_AZ','Bioavailability_Ma','hERG','AMES','DILI'])][['endpoint','prediction','raw_prediction','ensemble_sample_std','unit','status']]
        report+=f"### {cid}\n\n"+ '\n\n'.join(f"{x}. {r[x+'_'+name]}" for x,name in [('A','physchem_summary'),('B','druglikeness_summary'),('C','ADMET_summary'),('D','CNS_BBB_risk'),('E','toxicity_risk'),('F','synthesis_summary'),('G','stereo_salt_summary'),('H','evidence_gaps')])+f"\n\n{markdown(maincols)}\n\n预测溶解浓度：{r.solubility_mol_per_L:.3g} mol/L，{r.solubility_mg_per_L:.3g} mg/L（仅模型logS的代数换算）。\n\n"
    report+='''PACERGEN02432的高MW/脂溶性、较高TPSA及较低预测logS构成更明显的药化暴露挑战；其hERG概率较另外两者低，但仍不能归为低风险。PACERGEN00123与PACERGEN00462具有更低TPSA、更低MolLogP及更高FractionCSP3，BBB模型倾向更高；两者hERG及Ames概率提示优先验证需求。00462的BBB/logS模型均值略优于00123，但五模型分歧存在，不能把微小差异当作实验优越性。三者的P-gp抑制、CYP3A4抑制和DILI输出均值得确认；P-gp底物性、脑游离暴露及真实合成路线仍未解决。这里不计算候选总分、不排序、不重新定义PAM活性。

## 7. 证据缺口、后续实验与比赛展示

首要证据缺口：模型逐分子适用域与概率校准、Caco2单位冲突、pKa/logD、P-gp底物、真实Kp,uu/脑暴露、代谢位点、可靠半衰期/全身清除率、连续口服F、专门反应性代谢物模型、完整计算路线及真实采购/收率/成本。PPBR超界尤其说明不能盲信输出。可考虑针对冻结Top3补充实际溶解度（明确pH/盐/晶型）、MDR1外排、hERG、Ames、微粒体稳定性/代谢物与蛋白结合试验，并开展药化人员路线评审；这些仅为建议，没有在本次执行实验。

适合比赛展示：输入身份审计；RDKit性质分布、QED及透明规则计数；Top3结构及实算描述符；带“模型预测、AD未知、非实验”标记的logS/BBB/CYP/hERG等；未评价项与失败项组成的证据矩阵。不得展示为“已证实脑暴露良好”“安全无毒”“已完成合成路线”“低成本可采购”或“10月原始筛选已通过的门槛”。Figure 1–6均从落盘CSV生成，PNG 300 dpi及矢量PDF；图中的灰色为未评价，不是PASS。

## 8. 复现与文件

从仓库根目录运行（环境依赖精确版本见lock；需可访问公开PyPI）：

```sh
/usr/local/bin/python3 -m venv .venv-pacer-admet-v01
.venv-pacer-admet-v01/bin/pip install -r project/results/pacer_rerun_admet_synthesis_v01/requirements_admet_lock_v01.txt
/usr/local/bin/python3 -m venv .venv-pacer-retro-v01
.venv-pacer-retro-v01/bin/pip install -r project/results/pacer_rerun_admet_synthesis_v01/requirements_retro_lock_v01.txt
MPLCONFIGDIR=/tmp/pacer_mpl .venv-pacer-admet-v01/bin/python project/scripts/pacer_rerun_medchem_v01.py --output /tmp/pacer_admet_reproduction
```

后续阶段以相同`--output`运行pacer_rerun_admet_v01.py、pacer_rerun_synthesis_v01.py（retro环境）、pacer_rerun_report_v01.py（admet环境）。D阶段需要把本次EXTERNAL_ASSET_PREFLIGHT_v01.json审计快照复制至新目录，复现的是本次资产不可用条件；若未来资产可用，应另建版本、实现并验证正式搜索，不可把旧失败状态当作新网络状态。报告与图表阶段只读取保存的计算结果，不再次推理。脚本以独占创建方式防止覆盖已有结果。

输出包括规定的13个CSV/JSON/报告文件、补充端点级原始模型均值/五成员输出、结构警示、反应性motif、CNS单项、溶解浓度换算、模型元数据、阶段完成记录、两个依赖lock、Top3结构图、6组PNG/PDF。ALL200_ADMET与TOXICITY为长表，每个candidate_id×endpoint一行，prediction空值代表缺失/失败，raw_prediction保留已运行但不可采用的模型原值；0从不替代失败。

Git仅创建指定分析分支；未commit、push或创建PR。所有输出原本被.gitignore排除，版本管理采用仅针对本分析目录的精确例外；冻结目录与既有排序未修改。最终Git状态、文件SHA256和验证见ANALYSIS_VALIDATION与DELIVERY_MANIFEST。
'''
    with (o/'FINAL_SCIENTIFIC_REPORT_v01.md').open('x') as f:f.write(report)
    dump(o,'PHASE_E_COMPLETE_v01.json',{'status':'COMPLETE','figure_pairs':6,'top3_profiles':3,'complete_endpoint_records':len(ev),'report':'FINAL_SCIENTIFIC_REPORT_v01.md'})
    print('Phase E complete: six PNG/PDF pairs, complete profiles/report/evidence table.')
if __name__=='__main__':main()
