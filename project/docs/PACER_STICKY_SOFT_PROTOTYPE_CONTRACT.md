# PACER Sticky Soft Coupling Prototypes v0.2 开发合同

日期：2026-09-01

## 动机

v0.1有序soft prototype具有跨副本低JSD和ACh稳定性单调性，但瞬时membership时间
连续性不足。v0.2不调中心/K/带宽，而加入由训练轨迹自相关时间确定的sticky Markov
prior。

## 冻结方法

- 仍使用5个train-only分位数中心与相邻中心中位带宽；
- 每条训练replica计算coupling coordinate自相关首次降至1/e的lag（最多100 ns），
  取中位数tau；
- `p_stay=exp(-1/tau)`，其余概率只分配到相邻有序prototype；
- 高斯emission + forward-backward得到posterior membership；
- held-out数据不决定中心、带宽、tau或transition。

## Go

mean JSD<=0.30、posterior temporal continuity>=0.90、mean effective prototypes>=3、
mean max occupancy<=0.60、ACh稳定性单调性>=5/6折。功能预测仍No-Go。
