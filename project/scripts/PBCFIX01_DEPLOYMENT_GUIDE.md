# Stage 4 Rerun v01 PBC fix 01

本包是独立增量包，不包含新膜构建、AM1-BCC、Glide 重采样，也不修改旧版49文件清单或任何旧体系。
原始分支：results/pacer-stage4-rerun-preparation-v01。服务器继续使用已经安装并验证的原始 OpenMM 环境，不重新求解或升级力场环境。

## 本地上传

```powershell
scp project/results/pacer_stage4_rerun_v01_pbcfix01/pacer_stage4_rerun_v01_pbcfix01_repair_bundle.tar.gz project/results/pacer_stage4_rerun_v01_pbcfix01/pacer_stage4_rerun_v01_pbcfix01_repair_bundle.tar.gz.sha256 ubuntu@36.103.199.153:/data/pacer_stage4_rerun_v01_deployment/
```

## Linux 验证顺序（只修复、最小化、500步 smoke）

先激活服务器上原先成功执行单 Apo 实验的环境，并确保旧版任务已停止，避免快照期间旧文件变化。不要使用旧版构建、preflight 或 scheduler 命令。

```bash
set -euo pipefail
cd /data/pacer_stage4_rerun_v01_deployment
sha256sum -c pacer_stage4_rerun_v01_pbcfix01_repair_bundle.tar.gz.sha256
test ! -e pacer_stage4_rerun_v01_pbcfix01
tar -xzf pacer_stage4_rerun_v01_pbcfix01_repair_bundle.tar.gz
cd pacer_stage4_rerun_v01_pbcfix01
sha256sum -c SHA256SUMS_PBCFIX01.txt
unset CUDA_VISIBLE_DEVICES
python scripts/test_pacer_stage4_pbcfix01.py -v
python scripts/pacer_stage4_pbcfix01.py snapshot
python scripts/pacer_stage4_pbcfix01.py recover-single-apo --device 0
python scripts/pacer_stage4_pbcfix01.py smoke-single-apo --device 0
# 必须实际看到 PBC_FIX_SINGLE_APO_SMOKE_PASS 后继续。
python scripts/pacer_stage4_pbcfix01.py recover-all --device 0
python scripts/pacer_stage4_pbcfix01.py smoke-all
python scripts/pacer_stage4_pbcfix01.py preflight
sha256sum -c SHA256SUMS_OUTPUT.txt
find systems smoke single_apo evidence -name SHA256SUMS_OUTPUT.txt -print
```

单 Apo 与全12体系 smoke 使用不同输出目录。每个 smoke 包含50/100/200/300 K各125步、0.5 fs步长，禁用 barostat，检查每段势能、动能、坐标、受力有限及 CUDA 二进制 checkpoint 回读。所有12体系轮流使用物理 GPU 0/1；任一失败立即停止。每个新体系从旧 input.pdb 与旧 system.xml 的原始参数出发，仅改变已识别的位置约束表达式；独立证明其他系统参数未变化。旧 state.xml、minimized.pdb、refined_state.xml 只记录 SHA256，绝不作为重建起点。

首次快照核对旧包完整性与原构建输出哈希，记录所有旧文件哈希，并复制失败日志。新 build_audit.json 提供逐体系源文件哈希、修复前后 OpenMM 约束能量、独立最近周期像距离计算、最小化指标、电荷与离子数。每个输出目录有独立 SHA256SUMS_OUTPUT.txt；请在该目录内运行 `sha256sum -c SHA256SUMS_OUTPUT.txt`。

不覆盖任何已有输出。中断目录与 failure JSON 必须保留；不要删除以绕过门禁。若需重新运行，完整保留失败部署并另行制作新版本包与清单，不修改本版本的输入哈希。

## 生产前必须人工审查

preflight PASS 仅表示12个 CUDA smoke 与磁盘预算通过。科学审查与用户授权始终为 PENDING。本包不生成 authorization.json，不执行36个平衡任务或生产。

审查逐体系 net_charge_e（现有0/+1/+2/+3）、NA/CL计数、原配体参数、电荷态及与离子策略的兼容性；本修复不补离子或改变电荷。检查最小化最大受力、所有 smoke 阶段指标，确认任务的科学可接受性。

logs/storage_budget.json 由真实 CUDA checkpoint 文件大小、每100 ps的10帧全体系DCD、CSV/receipt开销估算：36条轨迹 × 100段 = 3600 checkpoint，总360 ns；增加50%余量与10 GiB预留。读取新目录所在数据盘的实时 total/free，容量不足则预检失败。报告不包含其他项目未来增长，投产前须再次检查 `df -h /data`。运行器在 prepare/production 前都会重新检查可用容量。

在用户明确授权且科学审查通过后，操作者才可以创建 authorization.json，字段为 protocol_id（本版本）、preflight_sha256（本版本预检文件SHA256）、scientific_review_approved=true、user_explicit_authorization=true、reviewer、user_authorization_reference（实际授权记录）。不要提前创建或伪造批准。

后续单独授权阶段才允许调用 scripts/schedule_pacer_stage4_rerun_v01_pbcfix01.py 的 prepare/production，二者均要求 --authorize-production 标志及上述审批文件。此处有意不提供可直接执行的投产命令。

新运行器仅读取新 systems/.../state.xml，不读取旧 refined_state.xml；prepare、production、ready收据、任务身份、分段 checkpoint 全部绑定新版本、system、PDB、起始 state、原始输入快照、包清单、预检与人工审批哈希。36 job_id 与全部路径也使用新版本标识；旧 checkpoint 无法通过新身份校验。

## 验证状态

本地回归测试结果见修复包旁 LOCAL_TEST_RESULTS.txt。用户提供的单 Apo 实验数值保存在 evidence/server_experiment_report.json，明确标记为用户提供证据。此次新包的12体系 CUDA、真实磁盘容量和逐体系电荷结果必须由 Linux 执行上述命令后产生，不能把既有实验结果当作本包 PASS。
