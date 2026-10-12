# PBCFIX01 下载与 PACER-FKG v02 部署验收 v01

验收日期：2026-10-11。当前状态：**REAL_INTEGRATION_PASS / FULL_ACCEPTANCE_PENDING**。
**本报告不标记 FULLY READY。** 真实两个分段和全部12个最终膜体系已经通过；
服务器快照尚只有12/36条轨迹完成，360 ns 全量数据、3600个分段及正式三阶段结果仍待验收。
没有修改服务器 MD 脚本、任务、checkpoint、科学参数或冻结资产，没有重新训练或提交 GitHub。

## 证据与验收状态

| 检查项 | 状态 | 实际证据及范围 |
|---|---|---|
| 真实部署根和目录 | 部分通过 | 已认证SSH读取 `/data/pacer_stage4_rerun_v01_deployment/pacer_stage4_rerun_v01_pbcfix01`，下载79文件/609277547字节；最终完整目录仍待36任务完成 |
| 36任务身份 | 通过 | 读取服务器真实矩阵；3候选×独立A/P/C/CP×R1/R2/R3，12独立体系，无共享A/P轨迹 |
| 服务器绝对路径 | 代码与回归通过 | 实际服务器矩阵目前使用相对路径；新增测试覆盖 output_path/checkpoint_path/system_path 的指定服务器绝对前缀；拒绝其他前缀、路径逃逸和旧路径 |
| 36轨迹/3600段 | 未完成 | 服务器采样快照：12 COMPLETE、2 RUNNING、22尚无progress，1353个提交回执；已核查这些回执的目录步数、连续顺序及身份，不代表已读取1353段DCD或checkpoint内容 |
| 真实最小分段输入 | 通过 | PACERGEN02432/A/R1 的00050000、00100000两个完成分段；各10帧、141470原子、10ps间隔；CSV、checkpoint、DCD与receipt SHA256一致 |
| 真实串联坐标/盒/时间 | 通过（两段） | 每一帧坐标SHA256及六维盒参数逐帧与独立DCD读取结果相同；时间10–200ps，跨段无丢帧、重复或重排；盒体积与CSV一致 |
| 全部12最终膜体系 | 通过 | 真实服务器minimized.pdb及build_audit；每个映射270残基、2139个ATOM14重原子；每个实际运行一次冻结Encoder，输出[1,270,256]且有限 |
| 两段Encoder真实重放 | 通过 | 20帧输出[20,270,256]；跨段重放索引0/9/10/19，relative RMSE约1.86e-7，冻结上限1e-4；不声称bitwise相同 |
| 冻结资产与依赖 | 通过 | 固定freeze SHA、归一化/RFF/图/区域逐文件认证；Geom2Vec安装代码与固定源码一致；checkpoint严格加载，161/161 keys，无缺失或额外key；CUDA实际推理成功 |
| 两段下游数值通路 | 诊断通过 | 实际20帧经stride5得到4采样点，执行冻结normalize/RFF/diffuse/regions；这是短窗诊断，不是正式20采样帧/1ns块，不产生Phase2a/2b完成回执，也未测试四情境真实对比 |
| 下载/Phase1恢复 | 通过（机制和回归） | rsync保留临时下载、最终SHA校验；新增OS互斥锁、孤立/损坏缓存隔离与下游回执失效机制；历史根拒绝写入，隔离文件不覆盖、不删除 |
| D盘/WSL/rsync本机可用性 | 未验证 | 本会话可访问文件系统只有C盘；真实样本落在独立C盘缓存。未实际向D盘执行rsync；下文提供WSL命令与Windows校验命令 |
| 正式360ns三阶段 | 未完成 | 需要36条完整数据下载后执行输入认证、Phase1、Phase2a、Phase2b与最终verify |

证据路径均相对于仓库根：

- `project/cache/pacer_stage4_pbcfix01_deployment_sample_v01/REMOTE_SAMPLE_DOWNLOAD_v01.json`：服务器目录、36任务快照、79个下载文件的源端SHA256及本地副本。
- `project/results/pacer_stage4_pbcfix01_fkg_v02_v01/deployment_acceptance/da46a8a65921/FINAL_12_TOPOLOGIES_v01.json`：12真实拓扑、原子数、冻结映射和Encoder结果。
- 同目录 `REAL_TWO_SEGMENTS_ENCODER_v01.json`：真实分段、边界重放和冻结数值通路。
- `input_audit/samples/da46a8a65921/`：逐帧坐标SHA256、盒信息、CSV体积及分段完整性。
- `deployment_acceptance/implementations/*/ASSETS_cuda_v01.json`：严格模型加载、依赖版本和固定资产SHA256；选择implementation与当前代码身份一致的记录。
- `deployment_acceptance/DEPLOYMENT_TEST_RESULTS_v01.txt`、`DEPLOYMENT_ACCEPTANCE_v01.json`：最终回归及范围总结。

验收记录使用实际 source-root 的哈希分目录，C盘样本与将来的D盘数据互不覆盖。
最终回归共50项，50通过、0跳过，包含已启用的真实两分段测试。
新下载缓存与SSH辅助依赖已加入精确范围的Git忽略规则，不自动提交原始数据或依赖包。
之前未下载真实数据时的 `local_validation/` 记录和中间检查结果仍保留，不能替代本报告的新证据。

## 下载后的精确目录契约

服务器源：
`/data/pacer_stage4_rerun_v01_deployment/pacer_stage4_rerun_v01_pbcfix01/`

Windows目标：`D:\data\pacer_stage4_rerun_v01_pbcfix01\`
二者下面的**相对目录原样保持**，不要多套一层同名目录：

```text
D:/data/pacer_stage4_rerun_v01_pbcfix01/
  SHA256SUMS_PBCFIX01.txt
  manifests/STAGE4_RERUN_JOB_MATRIX.csv
  manifests/STAGE4_RERUN_SYSTEM_MANIFEST.csv
  config/pacer_stage4_rerun_v01_pbcfix01.json
  scripts/...                         # 包清单全部成员，保留原始字节
  source.json
  DEPLOYMENT_GUIDE.md
  authorization.json
  evidence/source_snapshot.json
  evidence/...                        # 包清单成员
  logs/server_preflight_pacer_stage4_rerun_v01_pbcfix01.json
  systems/<system_id>/
    build_audit.json
    minimized.pdb
    system.xml
    state.xml
    ...                               # build_audit.outputs全部文件
  smoke/<system_id>/audit.json
  equilibration/<system_id>/replica_01/
    ready_pacer_stage4_rerun_v01_pbcfix01.json
    nvt/receipt.json, state.xml, state.csv, trajectory.dcd
    npt/...
    release_01/... through release_04/...
  production/<system_id>/replica_01/
    progress.json
    segments/00050000/trajectory.dcd, state.csv, checkpoint.chk, receipt.json
    segments/00100000/...
    ...
    segments/05000000/...
```

12个system各有replica_01、02、03。每段50000步、2fs，100ps；100个分段严格按照
50000、100000、…、5000000数值顺序读取。每段10帧，完整轨迹1000帧，10–10000ps。
1000帧的生产区间是10ns，首帧不是t=0。Phase2a保留索引0、5、…、995，时间10–9960ps。

真实DCD两个文件头都从约10ps重新开始，这与每段新建DCDReporter一致。
接口校验CSV绝对时间和回执步数后，只将读取时间设为`(global_frame+1)*10ps`。
不改变坐标或盒信息；不使用会丢弃重叠时间帧的`continuous=True`；不生成合并DCD。
PDB中部分溶剂resid字段产生解析警告；270残基使用唯一拓扑序列顺序及实际原子名验证，
12个真实体系均成功，没有使用缺失resid推测映射。

矩阵/ready中的路径若为服务器绝对路径，仅接受上述精确SERVER_ROOT前缀并映射到
`--source-root`；相对路径也支持。身份比对和SHA256始终使用原始矩阵与回执，**不要编辑下载文件来替换路径**。
`source.json`内旧轮source_root仅作原始证据，不会尝试打开旧服务器目录。

## 源端：生产全部完成后生成最终SHA256清单

以下命令在Linux服务器执行。只读生产树；清单写入生产目录外独立验收目录。
所有36条任务未完成时明确退出；不能把正在变化的progress作为最终下载快照。
已有v01清单时停止，保留既有证据；需要新快照请改验收目录版本。

```bash
set -euo pipefail
export SRC=/data/pacer_stage4_rerun_v01_deployment/pacer_stage4_rerun_v01_pbcfix01
export ACCEPT=/data/pacer_stage4_rerun_v01_deployment/fkg_download_acceptance_v01
mkdir -p "$ACCEPT"
python3 - <<'PY'
import csv, hashlib, json, os
from pathlib import Path
root=Path(os.environ['SRC']).resolve()
out=Path(os.environ['ACCEPT']).resolve()
assert out != root and root not in out.parents
ledger=out/'PBCFIX01_TRANSFER_SHA256_v01.txt'
listing=out/'PBCFIX01_TRANSFER_FILES_v01.txt'
assert not ledger.exists() and not listing.exists()
def relative(name):
    p=Path(name)
    if p.is_absolute(): p=p.relative_to(root)
    assert '..' not in p.parts
    return p
rows=list(csv.DictReader((root/'manifests/STAGE4_RERUN_JOB_MATRIX.csv').open()))
assert len(rows)==len({r['job_id'] for r in rows})==36
expected={}
for row in rows:
    base=root/relative(row['output_path'])
    progress=json.loads((base/'progress.json').read_text())
    assert progress['status']=='COMPLETE' and progress['completed_steps']==5000000 and progress['completed_ns']==10
    names={p.name for p in (base/'segments').iterdir() if p.is_dir() and p.name.isdigit()}
    assert names=={'%08d'%n for n in range(50000,5000001,50000)}
    for end in range(50000,5000001,50000):
        d=base/'segments'/('%08d'%end)
        receipt=json.loads((d/'receipt.json').read_text())
        assert receipt['start_step']==end-50000 and receipt['end_step']==end
        assert receipt['identity']==progress['identity']
        assert set(receipt['outputs'])=={'trajectory.dcd','state.csv','checkpoint.chk'}
        for name,digest in receipt['outputs'].items():
            p=(d/name).resolve()
            assert root in p.parents and p.is_file() and p.stat().st_size>0
            expected[p]=digest
files=sorted(p for p in root.rglob('*') if p.is_file()
             and not any('.interrupted.' in x or x=='__pycache__' for x in p.parts)
             and p.suffix not in ('.lock','.log','.tmp','.pyc'))
snapshot={p:(p.stat().st_size,p.stat().st_mtime_ns) for p in files}
tmp=ledger.with_name(ledger.name+'.partial')
with tmp.open('x',encoding='utf-8',newline='\n') as f:
    for p in files:
        assert not p.is_symlink() and root in p.resolve().parents
        h=hashlib.sha256()
        with p.open('rb') as source:
            for b in iter(lambda:source.read(8*1024*1024),b''): h.update(b)
        digest=h.hexdigest()
        if p.resolve() in expected: assert digest==expected[p.resolve()], 'receipt mismatch: '+str(p)
        name=p.relative_to(root).as_posix()
        assert '\n' not in name and '\r' not in name
        f.write(f'{digest}  {name}\n')
assert all((p.stat().st_size,p.stat().st_mtime_ns)==v for p,v in snapshot.items()), 'source changed'
os.link(tmp,ledger)
tmp.unlink()
with listing.open('x',encoding='utf-8',newline='\n') as f:
    for p in files: f.write(p.relative_to(root).as_posix()+'\n')
print('SOURCE_SHA256_PASS',len(rows),'jobs',len(expected)//3,'segments',len(files),'files')
PY
```

源端清单检查了全部3600个分段输出SHA256与生产receipt一致；DCD物理帧数、原子数、
时间轴、盒和冻结映射继续由目标端正式输入认证检查。
源端清单生成命令尚未对本轮全量执行，因为36条任务尚未完成。

## 正式rsync下载（Windows的WSL终端）

需要可访问D盘的WSL环境及rsync/ssh。密码只在SSH提示时输入，命令和清单不包含凭据。
`SRC/`末尾斜杠表示复制目录内容到目标根。清单明确选择文件；不下载活动日志、锁文件或失败归档。
不使用`--delete`。下载中断后原样重跑，rsync使用`.rsync-partial`恢复，并在最终文件发布前校验传输。
最终判据仍为SHA256和PACER输入认证。

```bash
set -euo pipefail
REMOTE=ubuntu@36.103.199.153
SRC=/data/pacer_stage4_rerun_v01_deployment/pacer_stage4_rerun_v01_pbcfix01
ACCEPT=/data/pacer_stage4_rerun_v01_deployment/fkg_download_acceptance_v01
DEST=/mnt/d/data/pacer_stage4_rerun_v01_pbcfix01
mkdir -p /mnt/d/data "$DEST"
rsync -rtvh -e ssh "$REMOTE:$ACCEPT/PBCFIX01_TRANSFER_SHA256_v01.txt" /mnt/d/data/
rsync -rtvh -e ssh "$REMOTE:$ACCEPT/PBCFIX01_TRANSFER_FILES_v01.txt" /mnt/d/data/
rsync -rtvh --info=progress2 --partial --partial-dir=.rsync-partial --delay-updates \
  --no-owner --no-group --omit-dir-times \
  --files-from=/mnt/d/data/PBCFIX01_TRANSFER_FILES_v01.txt \
  -e ssh "$REMOTE:$SRC/" "$DEST/"
cd "$DEST"
sha256sum -c /mnt/d/data/PBCFIX01_TRANSFER_SHA256_v01.txt
```

若只有少量数据可用于提前验证，可运行仓库提供的只读SFTP样本下载器。
它用交互式密码、已验证known_hosts，下载前两段、准备身份资产和12拓扑；
不会保存密码或修改服务器文件。传输辅助依赖与科学环境隔离：

```powershell
& C:/Users/barne/anaconda3/python.exe -m pip install --target project/cache/pacer_stage4_pbcfix01_fkg_v02_v01/ssh_runtime paramiko
& C:/Users/barne/anaconda3/python.exe -m project.pacer_fkg_v02.download_stage4_pbcfix01_sample_v01 --ssh-runtime project/cache/pacer_stage4_pbcfix01_fkg_v02_v01/ssh_runtime --destination project/cache/pacer_stage4_pbcfix01_deployment_sample_v01
```

目标已有内容不一致时样本下载器拒绝覆盖；请选择新样本目录保留旧快照。
正式360ns下载使用上面的rsync及最终清单，样本下载器不是全量下载工具。

## Windows目标端：可以直接复制的最终命令

从 `C:\projects\GPCR-virtual-screening` 仓库根执行。模型、源码与冻结目录使用已经实际验证的本地路径。
源端SHA清单放在campaign目录外的 `D:\data\`。

```powershell
Set-Location C:/projects/GPCR-virtual-screening
conda activate pacer_dc_geom2vec
$data = 'D:/data/pacer_stage4_rerun_v01_pbcfix01'
$env:PACER_STAGE4_FROZEN_ROOT = 'C:/projects/GPCR-virtual-screening-fkg-v02'

# 下载完整性、36任务与冻结资产
python -m project.pacer_fkg_v02.stage4_pbcfix01_deployment_v01 verify-sha256 --source-root $data --manifest D:/data/PBCFIX01_TRANSFER_SHA256_v01.txt
if ($LASTEXITCODE -ne 0) { throw 'Download SHA256 acceptance failed' }
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --contract --source-root $data
if ($LASTEXITCODE -ne 0) { throw '36-job contract failed' }
python -m project.pacer_fkg_v02.stage4_pbcfix01_deployment_v01 assets --device cuda
if ($LASTEXITCODE -ne 0) { throw 'Frozen assets or dependencies failed' }

# 在最终D盘副本上重做12拓扑、两段真实Encoder对接
python -m project.pacer_fkg_v02.stage4_pbcfix01_deployment_v01 topologies --source-root $data --device cuda
if ($LASTEXITCODE -ne 0) { throw 'Final topologies failed' }
python -m project.pacer_fkg_v02.stage4_pbcfix01_deployment_v01 sample-encoder --source-root $data --device cuda
if ($LASTEXITCODE -ne 0) { throw 'Real segment integration failed' }

# 正式36轨迹/3600段输入认证，以及全部冻结分析
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --all --source-root $data --device cuda
if ($LASTEXITCODE -ne 0) { throw 'Full PACER-FKG pipeline failed' }
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase2b_graph_region_v01 --verify --source-root $data
if ($LASTEXITCODE -ne 0) { throw 'Final frozen result verification failed' }
```

依赖检查已实际加载固定checkpoint
`b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`，
历史freeze为`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`。
本机torch CUDA可用；其他机器可以先用`--device cpu`检查资产，小规模测试后再选合适设备。
所有阶段的单独执行命令见 `PACER_STAGE4_PBCFIX01_FKG_V02_HANDOFF_v01.md`。

正式结果继续写入独立目录：
`project/results/pacer_stage4_pbcfix01_fkg_v02_v01/` 和
`project/cache/pacer_stage4_pbcfix01_fkg_v02_v01/`。
BS256缓存约9.27GiB、Phase2a约0.65GiB，建议额外预留至少12GiB分析空间，下载数据空间另计。
最终结果必须有36份Phase1、72份Phase2a、18份区域向量及通过的Phase2b回执。

## 可恢复运行及冻结成果保护

下载中断：重跑同一rsync命令，直至目标SHA检查通过；下载未完成不得跳过校验。
不要同时对同一数据目录下载和分析。输入认证检测读取期间源文件改变并拒绝发布有效回执。

正常Phase1失败：已验证的任务缓存可复用；本次尚未发布的临时缓存被清理，原样重跑`--all`。
强制终止或发布缓存后、发布manifest前崩溃：保留孤立缓存/`.partial.*`，运行明确的恢复命令：

```powershell
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --recover-phase1 --source-root $data
if ($LASTEXITCODE -ne 0) { throw 'Recovery failed; inspect preserved evidence' }
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --all --source-root $data --device cuda
```

恢复需要同一代码版本、同一完整源目录和有效正式输入回执。
它先验证全部路径，只移动本轮失败任务的缓存/manifest/partial到本轮`_recovery/<unique-id>/`，
并归档本轮Phase1完成回执和依赖它的Phase2a/2b结果及缓存；其他有效Phase1任务缓存保持可复用。
不会删除、覆盖或移动历史PACER-FKG目录、MD数据或模型。
恢复及分析共享OS文件锁，进程死亡自动释放，不靠手工删除锁文件绕过活跃任务。
代码版本、输入路径或冻结资产发生改变时，不会静默恢复旧科学结果；保留旧本轮输出后使用新的版本目录。

## FULLY READY的剩余条件

必须等待全部36任务完成，执行源端最终SHA、完整rsync、D盘目标校验和正式全量命令，
并确认所有真实回执及三阶段结果通过。当前已完成的是**真实最小对接与12拓扑验收**，
不能把两段短窗诊断、1353份回执目录检查或合成1000帧测试升级为3600段真实读取通过。
