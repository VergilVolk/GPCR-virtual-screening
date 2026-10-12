"""Render a Chinese scientific report from independently replayed CSVs."""
from pathlib import Path
import csv,json,hashlib
P=Path(__file__).resolve().parent
def read(n): return list(csv.DictReader((P/n).open(encoding='utf-8-sig')))
rows=read('STAGE4_PBCFIX01_REGION_COMPARISON_v01.csv')
cr=read('STAGE4_PBCFIX01_CANDIDATE_COMPARISON_v01.csv')
rel=read('CONTRAST_RELATIONS_v01.csv'); overlap=read('REGION_OVERLAP_DIFFUSION_v01.csv')
CS=['PACERGEN02432','PACERGEN00123','PACERGEN00462']; BS=['STATE_MOTION','SIGNED_DRIFT']; KS=['Delta_AGO','Delta_PAM','Delta_INT']
RS=list(dict.fromkeys(r['region'] for r in rows))
def f(x,n=3): return f'{float(x):.{n}f}'
def row(c,b,k,r): return next(x for x in rows if (x['candidate'],x['branch'],x['contrast'],x['region'])==(c,b,k,r))
def triple(x,field): return '/'.join(f(x[field.format(i=i)],3) for i in [1,2,3])
def pair(x): return '/'.join(f(x[a]) for a in ['R1_R2_cosine','R1_R3_cosine','R2_R3_cosine'])
def table(headers,rr): return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,x))+' |' for x in rr])
audit=json.loads((P/'ANALYSIS_VALIDATION_v01.json').read_text())
version=json.loads((P/'ANALYSIS_CODE_VERSION_v01.json').read_text())
lines=['# Stage4 PBCFIX01 PACER-FKG v02 独立科学分析 v01','',
'## 核心结论与结论等级','',
'**直接数值事实：** 三候选共同具有 STATE_MOTION 分支 PAM 接触区域的 ΔAGO、ΔPAM 重复方向一致性；该结论没有延伸到全受体、所有胞内微开关或两个分支。PAM 接触区域的 ΔINT 均出现至少一对负余弦。SIGNED_DRIFT 的 ΔPAM、ΔINT 没有呈现与上述 STATE_MOTION 接触信号相当的重复方向证据。幅度稳定和方向稳定必须分别判断。','',
'**探索性机制解释：** 冻结特征中存在候选相关、在无 probe 和有 probe 背景下均出现的 PAM 区域动态差异；它可能反映局部构象/波动分布的改变。STATE_MOTION 合并平均状态与运动强度，现有区域 RFF 向量不能分离两者贡献。没有形成三个重复都支持的、从接触区到胞内输出的完整一致链条。','',
'**待实验功能假说：** 三者可作为检验“候选自身效应、probe 背景调制、二者交互”的实验对象，但冻结 ΔAGO/ΔPAM 名称不等于 agonist/PAM 活性标签。本文不推断 PAM/ago-PAM 活性、效力、概率或因果机制，不给出综合药效排名。','',
'## 1. 输入身份、冻结协议及审计','',
'最终 JSON：`../pacer_stage4_pbcfix01_fkg_v02_v01/phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json`。其 SHA256 为 `ee3a6621b6df039c18cd762a5d2448cd60b1c69642e7131203df281a11558236`，与 Phase2b 最终冻结回执一致。全部18份 NPZ 的字节数和 SHA256 与 JSON 清单一致；候选×分支×R1/R2/R3 为完整 3×2×3 网格，路径、seed、cluster=1、A/P/C/CP 身份及上游 manifest 哈希均核对。未剔除、填补、调换任何重复。','',
f'独立复算了 {audit["n_arrays"]} 个有限的 (10,512) 数组、{audit["n_json_replica_region_summaries_replayed"]} 个 JSON 重复/区域摘要（包括两个辅助 probe 差分）及主要差分的 {audit["n_main_contrast_pair_cosines_replayed"]} 个重复对余弦。JSON 摘要复核容差 rtol=1e-10、atol=1e-12。上下文线性重构与冻结 contrast 的最大绝对误差 {audit["max_context_contrast_formula_absolute_error"]:.3g}；因分别扩散和 float32 存储产生舍入差异，按 rtol=2e-5、atol=2e-7 核对，不替换冻结数组。','',
'全量输入认证 `INPUT_AUTHENTICATION_v01.json` 状态为 `36_OF_36_ACCEPTED`，36条轨迹各100段，总计3600段。本文读取并记录该回执，不重新验证原始 DCD/checkpoint 内容。早期 `LOCAL_VALIDATION_v01.json` 和 `DEPLOYMENT_ACCEPTANCE_v01.json` 以及交接文档中的“未完成”属于先前部署/采样快照，不能代替最终输入认证与 Phase2b 回执，也不能将早期快照称为最终验证。本文验证范围是冻结结果身份与数值复算，不声称重新执行了原始 MD 到 Encoder 的最终流水线验证。','',
'冻结历史 manifest SHA256：`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`；图 SHA256：`25593c11dac0f70ad6811c8c6aa7a4f06c1755f14397f1162eaa1157026be122`。历史校准来自历史长 MD 的 R2；本次 Stage4 R2 与 R1、R3 一样仅用于评价。历史定义和数值资产只做读取/哈希核对。','',
'360 ns = 3候选×4上下文×3重复×10 ns；不是单条360 ns。每条1000原始帧，以 frames[::5] 取200帧；采样点10–9960 ps，间隔50 ps。每20帧形成1个名义1 ns块，块内首末跨度950 ps，每条10块。没有丢弃前段、重新划块或额外平衡期选择。10块彼此相关，仅作时间分布描述，独立重复数为3；全部区域亦非独立观测。','',
'## 2. 差分与统计量的解释','',
table(['名称','冻结公式','可作的描述'],[
['ΔAGO','C − A','无 probe 背景的候选相关动态差异'],['ΔPAM','CP − P','有 probe 背景的候选相关动态差异'],['ΔINT','CP − C − P + A = ΔPAM − ΔAGO','两个背景的候选差异之差']]),'',
'A=apo、P=probe_only、C=candidate_no_probe、CP=candidate_probe，全部取同一候选、同一重复匹配上下文。配对种子和时间块编号不代表两条随机轨迹的构象同步；ΔINT 的块幅度还包含上下文采样波动。ΔINT 非零不直接证明生物协同或非加和药效：RBF-RFF 本身非线性，其差分发生在冻结表示空间。','',
'STATE_MOTION = concat(mean(Z), sqrt(mean(diff(Z)^2)))；SIGNED_DRIFT = mean(diff(Z))，后者是块内50 ps采样步的平均有符号增量，经独立冻结归一化与 RFF 后的对比，不等于受体整体运动速度。两分支使用各自冻结尺度和随机特征，原始幅度不能跨分支比较，RFF 方向亦没有通用“激活方向”。正余弦表示同一差分在不同重复的特征方向对齐；负余弦表示方向冲突，不等于抑制/反向药理活性。','',
'令 v[r,b] 为一个重复的块向量，m[r]=mean_b v[r,b]。报告 ||m[r]||（时间平均向量幅度）、||v[r,b]|| 的 min/q10/median/q90/max，以及 R1–R2、R1–R3、R2–R3 的 cos(m[i],m[j])。幅度 CV 用三个 ||m[r]|| 的样本标准差(ddof=1)/均值。时间相干度 = ||mean(v)||/mean(||v||)，用于观察块间抵消；early/late 为前五块和后五块平均向量余弦，不能充当收敛检验。','',
'不做 p 值、置信区间、独立块 bootstrap、显著性宣称、数据驱动阈值、模型训练或综合权重。统计三个余弦是否全正仅作为几何符号事实（零点不是药效阈值）；接近零的正值不据此宣称强一致性。图中全部区域使用固定顺序。幅度图每个分支独立设置统一面板尺度。','',
'## 3. 区域依赖和特异性限制','',
'保留9个冻结区域：PAM consensus/union、cooperativity_mutagenesis、compound110_extension；正交 activation_core/contact_union；intracellular_microswitches；stable_core_control/distal_control。区域大小、每个残基的结构编号与 embedding_index 见 `REGION_MEMBERSHIP_v01.csv`。区域名称是结构先验，不能因名称将其动态变化认作该功能已发生。','',
table(['区域对','共同残基数','Jaccard','图扩散共享质量'],[[x['region_a']+' / '+x['region_b'],x['intersection_n'],f(x['jaccard']),f(x['diffusion_shared_probability_mass'])] for x in overlap if (x['region_a'],x['region_b']) in [
('pam_contact_consensus','pam_contact_union'),('pam_contact_consensus','cooperativity_mutagenesis'),
('orthosteric_activation_core','orthosteric_contact_union'),('pam_contact_consensus','orthosteric_activation_core'),
('orthosteric_activation_core','intracellular_microswitches'),('pam_contact_consensus','stable_core_control'),('pam_contact_consensus','distal_control')]]),'',
'图采用冻结 alpha=0.65、20步重启扩散，每步 H[t+1]=0.35I+0.65 T H[t]，H[0]=I。仅对结构算子计算区域平均影响权重 w，未对 MD/Encoder/RFF 数值重新运行扩散；这是评估共享信息的解析性诊断。共享质量=sum_j min(w_region1[j],w_region2[j])，范围0–1；它不是实测信号的独立性校正，也没有用于扣减、重加权或重新排名。','',
'PAM consensus 是 union 的子集；正交 activation core 也基本嵌在 contact union 内，cooperativity 与 PAM 区域重叠。因此多个此类区域一起重复，不是多份独立支持。即使直接无交集，图扩散仍共享信息；对照区域同样经历扩散。区域取均值降低残基数的简单求和效应，但区域大小、结构平滑程度与局部变异仍不同，区域范数比不能自动解释为选择性。','',
'特异性同时审阅两个对照：每个重复的 region/stable-control 与 region/distal-control 范数比均输出到区域 CSV，不跨分支、不先合并重复。远端对照并非静止：某些 ΔINT 的远端幅度超过 PAM 接触区域，限制“局部特异交互”的解释。`CROSS_REGION_DIRECTION_v01.csv` 另外给出同重复全部区域对的方向余弦，防止将平滑后的共同方向当作信号传播路径。','',
'## 4. 全部主要差分的候选横向比较','',
table(['候选','分支','差分','9区域中三对余弦均正的数量','最小余弦','三对最小余弦最高的区域'],[
[r['candidate'],r['branch'],r['contrast'],r['regions_all_three_pairs_positive'],f(r['all_pair_cosine_min']),r['most_directionally_consistent_region']] for r in cr]),'',
'该表展示空间分布，没有将重叠区域计数转化为证据强度分数。候选比较仅限相同分支、相同差分和区域；不按幅度、计数或任一余弦组合生成药效排序。各候选使用各自冻结 A/P 背景，不混用对照。','',
'### PACERGEN02432','',
'**支持（直接事实）：** STATE_MOTION 的 PAM consensus ΔAGO 余弦为0.562/0.424/0.744，ΔPAM 为0.742/0.729/0.583（顺序均R1–R2/R1–R3/R2–R3）。ΔPAM 在 consensus、union、cooperativity 的三对均正，orthosteric_activation_core 亦为0.225/0.349/0.274，但后者对齐较弱。ΔAGO 的 PAM consensus 幅度 CV=0.072左右，说明该局部差异在三个重复中幅度也较稳定；精确值见附表。','',
'**冲突（直接事实）：** PAM consensus ΔINT 为0.287/−0.182/−0.003，extension ΔINT 的R1–R3为−0.635。胞内 ΔPAM 为0.322/−0.158/−0.072，ΔINT 为0.426/−0.445/−0.052；局部差异未伴随可重复胞内方向。SIGNED_DRIFT ΔPAM 的 PAM consensus 为−0.191/0.121/−0.329。SIGNED_DRIFT 胞内 ΔINT 的三个幅度几乎相等（0.06578/0.06606/0.06620），方向却为0.306/−0.106/−0.043，是“幅度一致不等于方向一致”的具体例子。','',
'**探索性解释与边界：** 更符合可重复的 PAM 接触区状态/波动表示变化，并有较弱的正交区域背景依赖差异；不能证明协同交互、持续定向漂移或胞内激活。','',
'### PACERGEN00123','',
'**支持（直接事实）：** STATE_MOTION 的 PAM consensus ΔAGO=0.704/0.608/0.430，ΔPAM=0.624/0.315/0.692。cooperativity 的 ΔAGO=0.749/0.761/0.715，是本候选最一致的主要区域信号。PAM consensus ΔPAM 幅度 CV约0.094。胞内 ΔPAM 三对为0.122/0.189/0.413，虽全正但对齐较弱；胞内 ΔINT 的R1–R3仅0.003，不能用全正计数强化它。','',
'**冲突（直接事实）：** 正交 activation core 的 ΔAGO 三对均负（−0.443/−0.032/−0.046），ΔPAM 亦有R2–R3=−0.144。PAM consensus ΔINT=0.377/−0.407/−0.292，union R1–R3=−0.559。extension ΔINT=0.668/−0.689/−0.525，若只看R1–R2会遗漏严重冲突。distal ΔINT R2幅度0.1802，高于R1=0.06665及R3=0.09934，CV约0.507，方向又为−0.332/0.249/−0.637。SIGNED_DRIFT 的 PAM consensus ΔPAM=0.083/−0.085/0.035，缺乏重复对齐。','',
'**探索性解释与边界：** 接触/协同性结构先验区域的无 probe 差异较稳定；正交激活区和交互差分与之不一致。胞内弱正对齐不足以支持完整激活机制，更不能跳过远端R2差异。','',
'### PACERGEN00462','',
'**支持（直接事实）：** STATE_MOTION PAM consensus ΔAGO=0.820/0.640/0.814；cooperativity ΔAGO=0.868/0.749/0.819，其三对最小余弦高于其他两候选同区域同差分。ΔPAM PAM consensus=0.746/0.572/0.739，cooperativity=0.732/0.571/0.784。正交 activation core 的 ΔAGO=0.520/0.322/0.308，ΔINT=0.260/0.444/0.164，均全正但明显弱于局部接触区。SIGNED_DRIFT 正交 ΔAGO=0.134/0.177/0.237，提供弱方向一致事实，不能与 STATE_MOTION 幅度相加。','',
'**冲突（直接事实）：** PAM consensus ΔINT=0.052/−0.606/0.077、cooperativity R1–R3=−0.430。胞内 ΔAGO=−0.287/−0.141/0.197、ΔPAM=0.206/−0.258/−0.197、ΔINT=−0.278/−0.153/0.070，均混合。ΔPAM stable-control 的R3幅度0.06658约为R1=0.03290的2倍，CV约0.421；该重复的背景变化削弱选择性。ΔINT stable-control 三对为0.281/0.183/0.208，接近甚至超过部分正交交互对齐，进一步限制独特激活解释。SIGNED_DRIFT PAM consensus ΔPAM=0.080/0.205/−0.100。','',
'**探索性解释与边界：** 在同分支同区域比较中，它具有较整齐的接触/协同性区域 ΔAGO、ΔPAM 方向，并出现较弱的正交 ΔAGO/ΔINT 一致性；但 PAM 区域交互方向严重不稳定、胞内方向冲突、对照也变化，不能称其更可能有功能活性或已建立变构通路。','',
'## 5. 幅度、时间分布和交互几何复核','',
'附表列出全部三个重复的时间平均范数和三对余弦；CSV 进一步保留每个重复10块的分位数、极值、时间相干度和 early/late 余弦。范数是非负数，噪声或抵消也可产生非零范数。块范数 q10–q90 是观测分布，不是均值置信区间；三个重复的CV低也不能补偿向量方向相反。','',
'ΔINT=ΔPAM−ΔAGO 已逐块核对。`CONTRAST_RELATIONS_v01.csv` 输出同重复 AGO/PAM 夹角以及 ||INT||/(||AGO||+||PAM||)，用于区分“两个背景同向且近似”与“两个背景方向不同”。由于共享C/A/P/CP项，三差分不是独立证据；ΔINT 的更大幅度可能来自两个方向相减。现有10 ns重复和有限区域均值不能排除起始姿势、短时松弛、背景取样或随机特征表示的影响。','',
'## 6. 科研展示图与可追溯性','',
'每个图有同名PNG（220 dpi）和矢量PDF。英文图例避免跨平台缺字；图内数值对应CSV的原始精度字段，标注仅作四舍五入。全部来源见 `FIGURE_PROVENANCE_v01.csv`。','',
table(['图文件前缀','内容与用途'],[[b+'_replica_direction','3候选×3差分；全部9区域、全部3重复对，固定−1至1余弦色标'] for b in BS]+
[[b+'_magnitude_distribution','三个重复的块q10/median/q90与时间平均向量范数；仅同分支统一幅度轴'] for b in BS]+
[[b+'_region_specificity','各重复相对stable-control的范数比，并保留distal-control；1表示相等'] for b in BS]+
[['region_dependence','冻结区域残基Jaccard与结构扩散共享质量；评估信息依赖']]),'',
'`BLOCK_MAGNITUDES_v01.csv` 每行含候选/分支/差分/区域/重复/块时间及NPZ key，可追溯每个分位数。`INPUT_SHA256_v01.csv` 保存输入绝对路径、字节数、哈希及来源预期哈希；分析结束后已复核所读取输入哈希不变。`ANALYSIS_VALIDATION_v01.json` 给出检查计数及数值容差。','',
f'分析代码版本 v01，`analyze_v01.py` SHA256 `{version["script_sha256"]}`，仓库基线 HEAD `{version["git_head"]}`。代码与环境完整记录在 `ANALYSIS_CODE_VERSION_v01.json`，报告渲染脚本为 `write_report_v01.py`。运行命令见该记录，全部输出仅写入本独立目录；没有导入或执行上游模型流水线，没有改动历史结果，也没有提交GitHub。','',
'## 7. 最终回答：动态证据、冲突及实验优先事项','',
'**02432：** 最可重复的是 STATE_MOTION 的 PAM 接触区 ΔPAM，以及局部 ΔAGO；正交 ΔPAM 较弱但三对同向。主要冲突是 PAM ΔINT、胞内 ΔPAM/ΔINT 和 SIGNED_DRIFT 的缺乏一致性。优先检验候选在有/无 probe 背景下是否产生可重复功能变化，尤其检验“有背景差异但未形成可重复交互”是否只是动态表示层面的现象。','',
'**00123：** 最可重复的是 STATE_MOTION cooperativity ΔAGO，PAM consensus ΔAGO/ΔPAM 亦重复。主要冲突是正交 ΔAGO、PAM/extension ΔINT 与远端R2异常突出。优先区分候选自身作用、probe背景调制与非特异响应；胞内弱一致方向只作为次级探索，不应据此提前判定ago-PAM。','',
'**00462：** 最可重复的是 STATE_MOTION PAM/cooperativity 的 ΔAGO和ΔPAM；同区域ΔAGO方向一致性比其余两者更整齐。正交 ΔAGO/ΔINT 有较弱一致支持。最明显冲突是 PAM consensus ΔINT 的R1–R3=−0.606、胞内所有主要差分方向不一致，且stable-control也有变化。优先验证候选自身功能和probe背景调制是否存在，并检验局部接触变化能否对应可重复受体输出；这些优先事项不是功能排序。','',
'**实验设计建议（均为待验证假说）：** 对三个候选采用相同、预先定义的功能实验矩阵：vehicle、probe单独、候选单独、候选+probe；测候选浓度系列和probe完整浓度–反应曲线，分别评估基线效应、probe曲线位移和最大响应变化，不仅测单一浓度。保留跨日独立重复及同批对照，用受体相关功能读数，并加入受体缺失/阻断、细胞状态与检测干扰对照以区分非特异效应。实验是否支持PAM或ago-PAM应由预先约定的实验判据决定，不从本报告学习阈值或权重。','',
'若基础功能矩阵重复，再将PAM接触/协同性区域与正交/胞内区域的结构假说做定点扰动或其他正交验证；选点应遵循冻结结构先验并预先注册，结合表达/基线功能对照。不要从区域RFF向量推出某残基是因果位点，不用单次突变效应直接宣称变构传播路径。三个候选都值得在一致实验条件下检验；当前结果支持的是不同程度的局部动态差异重复性，而不是已验证的PAM/ago-PAM药理结论。','',
'## 附表：全部9区域、两个分支、三主要差分','',
'每行范数为R1/R2/R3的 ||mean(v)||，余弦顺序为R1–R2/R1–R3/R2–R3。表中数值三位小数；CSV保留全精度，接近零的符号以CSV为准。','']
def span(r,field):
    vals=[float(r[field.format(i=i)]) for i in [1,2,3]]
    return f'{min(vals):.3f}–{max(vals):.3f}'
extra=['### PAM consensus：幅度分布与对照审阅','',
       '以下每个范围是三个重复各自统计量的最小至最大值，不将30个相关块合并成独立样本。q10最小与q90最大给出三个重复描述区间的外包络，不是置信区间。', '',
       table(['候选','分支','差分','块q10最小–q90最大','时间相干度范围','前/后半段余弦范围','相对stable范数比范围','相对distal范数比范围'],[
           [c,b,k,f"{min(float(row(c,b,k,'pam_contact_consensus')[f'R{i}_block_norm_q10']) for i in [1,2,3]):.3f}–{max(float(row(c,b,k,'pam_contact_consensus')[f'R{i}_block_norm_q90']) for i in [1,2,3]):.3f}",
            span(row(c,b,k,'pam_contact_consensus'),'R{i}_temporal_coherence'),
            span(row(c,b,k,'pam_contact_consensus'),'R{i}_early_late_cosine'),
            span(row(c,b,k,'pam_contact_consensus'),'R{i}_ratio_to_stable_core_control'),
            span(row(c,b,k,'pam_contact_consensus'),'R{i}_ratio_to_distal_control')]
           for c in CS for b in BS for k in KS]),'',
       'STATE_MOTION 的局部接触信号相对stable控制增强，并不代表相对distal控制始终增强。SIGNED_DRIFT 的接触区范数可与stable控制相近、低于distal控制，且方向常混合；控制区背景变化不能忽略。时间相干度与前/后半段余弦仅描述同一10 ns轨迹内部变化，不解决独立重复冲突，也不能证明达到稳态。','']
idx=lines.index('## 6. 科研展示图与可追溯性')
lines[idx:idx]=extra
lines=[x.replace('CV=0.072左右','CV=0.071左右').replace('CV约0.421','CV约0.423') for x in lines]
for c in CS:
    for b in BS:
        lines.extend(['### '+c+' / '+b,'',table(['差分','区域','范数 R1/R2/R3','余弦 12/13/23','幅度CV'],[
            [r['contrast'],r['region'],triple(r,'R{i}_mean_vector_norm'),pair(r),f(r['replica_norm_cv'])]
            for r in rows if r['candidate']==c and r['branch']==b]),''])
(P/'STAGE4_PBCFIX01_SCIENTIFIC_ANALYSIS_v01.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('Report written')
