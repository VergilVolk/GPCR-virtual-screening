# PACER-DC Geom2Vec 编码器小试与队友交接

## 结论先行

现有四上下文 MD **全部复用，不需要重跑**。本分支只替换：

`atom14 trajectory window -> OneProt embedding`

为：

`atom14 trajectory window -> frozen Geom2Vec ViSNet residue embedding`

本地 CPU 已完成真实 M4/compound110/replica 1 小试。当前只能确认工程可用和局部信号未塌缩，尚不能确认跨 replica 优于 OneProt；R2/R3 是必须完成的晋级审计。

## 已实测结果

- 官方 `visnet_l6_h64_rbf64_r75.pth`，SHA256：`b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`；
- 真实输入：100 帧、270 残基、2,139 个重原子；
- 输出：每帧 `270 x 128` 旋转/平移不变量残基特征；
- CPU纯PyTorch邻域后端：约 `3.0–4.8 s/frame`；
- identical repeat：相对误差 `0`；
- 旋转：相对误差 `0`；
- 平移：相对误差 `4.42e-7`；
- compound110 pocket 的 context separation / temporal variation = `1.18`；
- global = `0.68`，distal control = `0.68`。

后一点说明局部口袋信号没有像OneProt最终头那样被全局平均完全抹除，但这只是单replica描述性证据。

## 安装

权重已直接随本分支提交，队友不需要另行下载。先核验：

```bash
sha256sum project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth
```

应得到：`b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`。

```bash
conda env create -f project/environment_pacer_dc_geom2vec.yml
conda activate pacer_dc_geom2vec
```

若Linux/CUDA环境已有匹配版本的`pyg-lib`或`torch-cluster`，代码自动使用官方加速邻域构图；否则使用仓库内确定性PyTorch fallback。初次运行建议`--batch-size 1`。

## 队友现在运行什么

在保存R2/R3 atom14窗口的机器上：

```bash
python project/pacer_dc_training/run_geom2vec_four_context_batch.py \
  --input-root project/results/pacer_dc_four_context_v01/compound110 \
  --output-root project/results/pacer_dc_geom2vec_R2_R3_v01 \
  --checkpoint project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth \
  --replicas 2 3 \
  --stride 1 \
  --batch-size 1 \
  --device cuda
```

不要重跑MD，不要训练head，不要改region，不要读取blind candidates。先把`batch_audit.json`和全部`.geom2vec.npz`推回同一分支。

提取完成后运行同协议R2/R3审计：

```bash
python project/pacer_dc_training/audit_geom2vec_cross_replica.py \
  --embedding-root project/results/pacer_dc_geom2vec_R2_R3_v01 \
  --region-map project/results/pacer_dc_four_context_v01/compound110/G2_REGION_MAP_v02.json \
  --replicas 2 3 \
  --windows 0 1 2 3 4 \
  --output project/results/pacer_dc_geom2vec_R2_R3_v01/cross_replica_audit.json
```

## 晋级条件

使用与OneProt完全相同的R2/R3五窗口协议计算区域级：

`dAGO = z_C - z_0`

`dINT = z_CA - z_A - z_C + z_0`

只有同时满足以下条件才允许把Geom2Vec改为主编码器：

1. 生物学预定义区域的跨replica中位cosine高于OneProt对应层；
2. bootstrap 95% CI下界大于0；
3. context/replica方差比大于1；
4. 改善集中在ACh/PAM口袋或激活通路，不能在distal control同幅改善；
5. 至少在LY2119620纯PAM和compound110变构激动剂上验证`dAGO`方向差异。

## 主张边界

当前结果不证明Geom2Vec优于OneProt，不证明compound110是PAM，也不证明可以预测PAM效力。它只证明：现成MD可以无缝接入冻结Geom2Vec，编码器具备正确几何不变性，并保留了初步局部上下文差异。

## 提取完成后的新分析入口

不要再逐帧相减。R2/R3回传后，先用冻结的多结构图运行PACER-FKG：

```bash
python project/pacer_dc_training/pacer_factorial_kernel_graph.py \
  --embedding-root <单个candidate/replica/window的四上下文目录> \
  --graph project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json \
  --output <对应的PACER_FKG审计.json>
```

统计汇总必须以replica为独立单位做block bootstrap；不得把window或frame当作独立生物学样本。方法定义、基线和晋级条件见`docs/PACER_FACTORIAL_KERNEL_GRAPH_METHOD.md`。
