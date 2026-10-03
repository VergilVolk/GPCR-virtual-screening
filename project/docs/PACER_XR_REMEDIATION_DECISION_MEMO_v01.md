# PACER-XR 六路补救方案决策备忘（v01，2026-10-04）

依据：`INCIDENT_PACER_XR_PROTOCOL_DEVIATION_20261003.md`（事故冻结）、
`PACER_XR_CASCADE_METHOD.md`（六路协议）、本机引擎审计。

## 一、事实核（先纠三个错误认知）

1. **"前三是两模型共同最高分加权出来的"——错。** 最终三候选（PACER0010/
   0073/0027）= 2023 M4-safe DrugCLIP 排名（1/2/4）+ Murcko 骨架多样性，
   XR 从未参与。双模型交集 5 分子（old≤50 ∧ new≤25）是**早期链**，事故
   报告明确其"不能作为最终三候选晋级来源"。
2. **"重跑中间两模块还是这三个"——跑之前不可断言。** XR 六路分数对
   200 分子完全空白（连历史 24 候选的六路表都是空的，
   `awaiting_protocol_matched_six_channel_scores`）。PACER0014（rank 3，
   被骨架压缩挤掉）与 PACER0027 在 XR 融合下次序未知。
3. **"大不了凑一下（靠 LSTM/BRICS 随机漂移）"——禁止。** 这是事后候选
   操纵，违反事故报告 §5 自己冻结的规则；被评委复现即翻车。

## 二、六路缺口与可运行性

| 通道 | 需要什么 | 本机状态 | 判定 |
|---|---|---|---|
| Glide_PDB/BEmin/BEavg | Schrödinger 许可 + 协议匹配运行 | 无许可（历史即无） | **硬阻塞** |
| Vina_PDB | 官方晶体受体 + 协议匹配参数 | 分布表在库（M4R_PDB_Vina_scores.csv）；vina 在队友环境 | 可跑（需协议匹配重跑） |
| Vina_BEmin/BEavg | 官方系综协议匹配运行 | 分布表在库（115,771 行 ×2）；**既有 10-cluster 数据不算协议匹配**（box/协议不同，事故报告已定） | 可跑（需重跑） |

## 三、选项（按推荐排序）

### A. 拿到 Glide 许可 → 完整六路 post-hoc supplementary（最忠实）
- 学界/云按需许可；跑 6 路 → empirical percentile（官方分布已在库）→
  六路等权 + Glide-BEmin top-1% 级联 → 与既有三候选对比。
- 标签：**post hoc supplementary PACER-XR analysis**（事故报告 §5 措辞），
  明示候选已预先选定，不倒写为 prospective 证据。
- 风险：许可时间不可控；deadline 前赌不起。

### B. 无许可 → 冻结 XR v02：开源引擎替换 Glide 路（deadline 务实）
- 替换示例：GNINA CNN rescore / smina vinardo / QVina2 ×3 替 Glide ×3；
  Vina ×3 按协议匹配重跑。
- **必须**在官方 benchmark 分层子集上重验 v02 级联性能（标签在库），
  重冻结方法文档，所有材料标注"Glide→开源替换，v02"。
- 风险：v02 性能可能低于 0.7775 的原六路数字；讲稿里原数字只能标
  "原方法（含 Glide）在官方 benchmark 的成绩"，不得移花接木。

### C. 只跑三路 Vina 共识 post-hoc（最小计算量）
- 协议匹配跑 Vina ×3 → 三路 rank 共识 → 对照三候选。
- 标签："post hoc 三通道 Vina 共识"，**不得称六路 XR**。
- 六路方法创新在提交中定位为"benchmark 已验证的方法学"
  （AUC 0.7775 数字真实、可经库内脚本+已发布分数表复现），与本次
  prospective 链的执行差异如实写入 incident 报告引用。

### D. 不补算，叙事修正（零计算）
- 方法创新与执行链分开陈述 + 事故报告作为完整性证据。
- 评委复现风险：如果讲稿/PPT 声称候选"经六路 XR 筛出"则必翻车；
  只要不这么声称，D 是最低风险方案。

## 四、Stage4 投入保护

36 轨迹 ×10ns = 360ns 已跑在三候选上。若补救后 top-3 变化：
- 不重跑 Stage4：如实陈述"选择链按 DrugCLIP+系综 Vina 执行，XR 为
  post-hoc 补充分析"，360ns 投入保留；
- 重跑 Stage4（360ns，天级）：仅当团队坚持"XR 六路为唯一正式选择规则"
  且时间允许。

## 五、无后悔动作（立刻做）

1. PPT/讲稿数字链（附图注用）：
   17,146 BRICS joins → 5,000 生成 → 2,605 过生成规则 → Stage2 排除
   2,405（13 exact-known / 931 domain / 15 risk / 1,446 quota）→
   **200 入栈（local 160 + exploratory 40）** → 2,000 Vina 任务
   （200×10 构象）→ 200 结构门全 PASS → 140 骨架代表 → **3 入 MD**
   → 36 轨迹 360ns → PACER-FKG 冻结评估。
2. 冻结物不动：三候选 pose SHA、manifest、Stage4 receipts、事故报告。
3. 任何补算一律挂"post hoc supplementary"标签 + 独立时间戳 provenance。

## 六、建议

deadline 压力下：**先并行启动 C（三路 Vina，计算小时级）保底**；同时
追问学校 Glide 许可能否 72h 内到位，能则升级 A；否则按 C+D 组合交付，
B 仅在有余力重验 v02 时考虑。"凑一下"方案否决。
