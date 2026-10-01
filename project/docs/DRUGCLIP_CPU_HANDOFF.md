# DrugCLIP × PACER-DC CPU交付说明

## 这条分支交付什么

本分支把官方DrugCLIP作为PACER-DC的预训练结合表征层，并提供纯CPU可复现入口。它回答
“分子与M4别构口袋是否具有预训练模型意义上的结合相容性”，不把DrugCLIP相似度冒充PAM
功能效力。

本机已完成最小实测：官方权重成功加载；权重缺失键和异常键均为0；真实7TRS M4别构
口袋得到1×128向量，3个真实M4分子得到3×128向量，分数全部有限。529分子正式基准必须
由运行者按下述冻结流程完成后再报告。

## 固定版本与大文件

- DrugCLIP：`bowen-gao/DrugClip@7a3a3fa33673f8668c811790f2e4681c98af44ef`
- Uni-Core：`dptech-corp/Uni-Core@44f6386f4dcd7137fc1e5d5e768117d635d64a26`
- 官方checkpoint文件名：`checkpoint_best.pt`
- 文件大小：`1183713459` bytes
- SHA256：`dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e`

第三方源码、checkpoint、LMDB、embedding和结果都不进入Git。checkpoint从DrugCLIP官方
数据链接取得，或由队内点对点传输后放到：

```text
project/tools/DrugCLIP/artifacts/checkpoint_best.pt
```

## 设备要求

- Windows 10/11 + WSL2 Ubuntu；
- 纯CPU即可；建议16 GB以上内存；
- 首次安装约需数GB磁盘空间；
- 如网络较慢，优先由一名队员完成环境和权重准备，再共享校验一致的checkpoint，不要上传Git。

## 1. 安装CPU环境

进入WSL并切换到仓库根目录，然后执行（脚本会自动识别仓库位置）：

```bash
bash project/scripts/setup_drugclip_cpu_wsl.sh
```

脚本会固定第三方commit、创建`drugclip-cpu`环境、应用仅涉及设备位置的CPU兼容补丁，并以
`--disable-cuda-ext`安装Uni-Core。官方checkpoint内部仍使用旧注册名`binding_affinity`，
提取脚本会明确映射到公开代码中的`drugclip`，不会修改网络结构或权重。

## 2. 一条命令运行冻结baseline

```bash
bash project/scripts/run_drugclip_m4_cpu_baseline.sh
```

该入口依次构建529个分子、7TRQ/7TRP/7TRS三个固定口袋，提取embedding并生成
`project/results/drugclip_m4_benchmark_v01/baseline.json`。三个口袋使用同一组预声明残基，
禁止根据标签调整口袋或只报告最优结构。若本地没有三个PDB，入口会从RCSB下载并核验
冻结的SHA256；PDB本体不塞进Git。

验收条件：`molecule_shape=[529,128]`、`pocket_shape=[3,128]`、`finite=true`，且
`missing_checkpoint_keys=[]`和`unexpected_checkpoint_keys=[]`。

必须完整报告7TRQ、7TRP、7TRS、state-mean和state-max，不得根据测试标签只挑最好的口袋。
该结果是无监督、回顾性的结合相容性baseline，不是PAM功能验证。

## 3. 接入四上下文与Triplet

DrugCLIP产生每个上下文的`dc_p_*`和每个分子的`dc_m_*`；OneProt-MD产生`op_*`。
`build_drugclip_context_features.py`将它们组成：

```text
Fcontext = [OneProt-MD, DrugCLIP-pocket, pocket×molecule, cosine]
dPAM = F(candidate+ACh) - F(ACh-only)
dAGO = F(candidate-only) - F(apo)
```

随后由`train_dual_context_heads.py`训练轻量功能头。DrugCLIP高分、实验上non-PAM的分子只作为
PAM分支困难负样本；AGO分支保持独立负样本定义。完整创新必须在相同chemotype-held-out拆分
下优于DrugCLIP、OneProt-MD和无triplet融合三个baseline，才能晋级。
