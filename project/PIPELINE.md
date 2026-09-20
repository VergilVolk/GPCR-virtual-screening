# CHRM4 PAM 项目主流程

## 核心问题

本项目不把“能够进入 M4 别构口袋”等同于“具有 PAM 功能”。需要判断的是：候选物是否在 ACh 存在时增强受体耦合，同时在没有 ACh 时不过度独立激活受体。

当前唯一主线为 **PACER-DC：有、无 ACh 对照下的 M4 PAM 动态评价**。生成、QSAR、对接和 ensemble docking 均为候选产生或基线，不再与主算法并列。

论文级主张、统一 baseline、统计检验和晋级标准冻结于
`docs/PACER_DC_PUBLICATION_PROTOCOL.md`。后续不得根据候选结果事后修改主要端点或成功阈值。

## 流程

1. **候选产生**：文献分子、生成分子和商业库共同提供候选。
2. **结合门控**：理化性质、适用域、单结构/ensemble docking 排除明显不合理分子。该阶段不预测 PAM 功能。
3. **四上下文 MD**：对同一候选建立 `candidate+ACh`、`candidate-no-ACh`，并使用匹配的 `ACh-only`、`apo` 对照；至少三个配对 replica。
4. **轨迹表征**：先计算可解释的距离、接触和 ACh 稳定性；再比较 PCA/tICA、VAMP/SPIB 等轨迹编码方法。
5. **双上下文差分**：分别计算 ACh 依赖的协同变化和无 ACh 时的内在激动风险。
6. **药理对比微调**：仅在具有多 chemotype、匹配上下文和明确功能标签后使用 triplet；hard negative 包括别构激动剂和可结合但无 PAM 功能的分子。
7. **严格验证**：以分子、化学系列或独立体系为切分单位；轨迹帧不得随机跨集合。与 2D QSAR、静态 docking、ensemble docking、手工 MD 指标和无 triplet 表征比较。
8. **候选输出**：只输出多维证据和不确定性。未经功能实验确认的分子称为“计算优先候选”，不称为已发现 PAM。

## PACER-DC 输出

- `CoupledShift`：`candidate+ACh` 相对 `ACh-only` 的耦合变化；
- `OrthostericStabilization`：候选物对 ACh 构象/接触稳定性的改善；
- `IntrinsicActivationRisk`：`candidate-no-ACh` 相对 `apo` 的激活样变化；
- `BindingCompatibility`：候选物在别构口袋中的稳定性与相互作用；
- `Uncertainty/AD`：replica 不确定性、化学适用域和缺失上下文。

## 当前闸门

- **结构筛选**：可运行，但只能作为结合门控。
- **PACER-DC 计分器**：已实现并有基础测试。
- **四上下文生产 MD**：尚未启动；主机缺少可重复的候选配体参数化环境。
- **Triplet/VAMP/SPIB 训练**：未达到数据闸门，禁止用轨迹帧数量冒充独立样本量。
- **功能性 PAM 结论**：必须经过前瞻功能实验；当前候选均为计算假设。

## 固定判废标准

1. 完整方法不能在留一体系/留一 chemotype 中超过简单 MD 差分基线；
2. 已知 PAM 与别构激动剂不能被无 ACh 上下文区分；
3. 结论依赖单个 replica、初始瞬态或任意聚类数；
4. 缺失任一必要上下文却仍输出 PAM 概率；
5. 将 docking、静态结构或一条轨迹的结果解释成功能效力。

自动闸门：

```powershell
py -3.13 project/scripts/build_pacer_dc_v2_manifest.py
py -3.13 project/scripts/audit_pacer_dc_phase1.py `
  --manifest project/config/pacer_dc_pilot_md_manifest_v2.csv `
  --outdir project/results/pacer_dc_phase1_audit_v2
```

结果写入 `project/results/pacer_dc_phase1_audit_v2/`。
