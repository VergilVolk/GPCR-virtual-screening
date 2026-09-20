# PACER Soft Coupling Prototypes v0.1 开发合同

日期：2026-09-01

## 定位

该模块是在看到hard-cluster不稳定与连续轴-ACh稳定性关联后定义的回顾性方法开发，
不是独立确认。它不再学习无序KMeans标签，而将冻结连续PAM-coupling coordinate
转换为可微、有序soft prototype membership。

## Leave-one-replica-out

- 每折只用五条训练replica的coupling coordinate确定5个中心：10/30/50/70/90%
  分位数；
- 带宽固定为相邻中心间距中位数；
- membership为各中心高斯距离的softmax；held-out replica只投影；
- 不使用ACh坐标确定中心或带宽。

## 指标

- held-out与train pooled prototype occupancy JSD；
- soft temporal continuity：相邻帧membership的`1 - total variation`；
- effective prototype count、max occupancy；
- 每个prototype的ACh RMSD加权均值与中心的Spearman单调性。

## 开发Go

mean JSD<=0.30，soft temporal continuity>=0.90，effective prototypes>=3，max occupancy
<=0.60；ACh稳定性单调性至少5/6折为正。即使通过也只晋级为轨迹表示/attention
head，不晋级为PAM效力预测器。
