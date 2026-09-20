# PACER Static-Anchored Trajectory Test v0.1

日期：2026-09-01（结果前冻结）

## 假设

若不同M4 PAM共享部分功能耦合构象，则两个独立PAM共晶状态7TRQ（VU0467154）与7TRP（LY2033298）的口袋蛋白几何应存在共同方向，并且第三种PAM MK-97 的GaMD轨迹应比ACh-only 7TRS更接近这两个PAM锚点。

## 固定分析

- 仅使用作者GaMD聚类口袋中的蛋白残基，不使用不同化学组成的配体原子；
- public GaMD residue编号按已验证偏移 `static M4 residue = public residue - 598` 映射；
- 只保留7TRQ/7TRP/7TRS/public topology中共同的重原子 `(residue, atom name)`；
- 每帧及三个静态结构均以共同口袋原子刚体对齐到public第一帧；
- 定义 `PAM-likeness = RMSD(frame,7TRS) - mean[RMSD(frame,7TRQ),RMSD(frame,7TRP)]`；正值表示同时更接近两个PAM锚点而非ACh-only；
- 不用静态锚点训练trajectory encoder，不调整cluster；只对现有public k=10 cluster作事后生物学注释。

## 支持条件

1. `RMSD(7TRQ,7TRP)`小于两条PAM–ACh-only距离；
2. 六条MK-97 replica的PAM-likeness中位数均>0；
3. 至少一个预存cluster在六条replica均有样本且PAM-likeness稳定为正。

该测试只有三个静态结构锚点，最多支持“共享几何方向假设”；不能证明功能微状态、因果协同性或候选PAM效力。

