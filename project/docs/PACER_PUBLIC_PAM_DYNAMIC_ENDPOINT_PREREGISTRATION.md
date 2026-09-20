# PACER Public PAM Dynamic Endpoint v0.1 预注册

日期：2026-09-01

## 数据

Wang et al. 2022官方Source Data，Supplementary Figure 8：M4R-iperoxo与
M4R-iperoxo-LY2119620各3条独立1 us MD，5000个iperoxo RMSD点/轨迹。

## 单位与假设

- 独立统计单位是trajectory replica（3 vs 3），不是30000个时间点；
- 5000点覆盖1 us，按均匀0.2 ns间隔解释；
- 主效应：`mean RMSD(ip-only) - mean RMSD(ip+LY)`，正值代表PAM稳定正构探针；
- 主分析使用完整production，去除首100 ns为敏感性分析；
- 50 ns block（250点）保留时间相关性。

## 推断

1. 报告每条replica均值/中位数/标准差；
2. 3 vs 3 replica mean的20种精确标签置换，单侧检验PAM RMSD更低；
3. 分层block bootstrap：先各条件重采样replica，再在replica内重采样50 ns blocks；
4. 不报告帧级AUC，不把N=6完全分离称为高性能分类。

## 支持条件

- 三个PAM replica的均值均低于三个no-PAM replica均值中的相应配对；
- 单侧精确置换p<=0.05；
- block-bootstrap 95% CI下界大于0；
- 去除首100 ns后效应方向保持。

通过只支持“LY2119620在该公开体系中动态稳定iperoxo”，不证明一般PAM效力预测。
