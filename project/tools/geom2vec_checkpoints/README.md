# Geom2Vec frozen checkpoint

本目录包含 PACER-FKG 当前冻结使用的现成预训练权重：

- 文件：`visnet_l6_h64_rbf64_r75.pth`
- 大小：`1,815,684 bytes`
- SHA256：`b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`
- 上游项目：<https://github.com/dinner-group/geom2vec>
- 固定源码版本：`371d642ec1061664f16e49fcac702d07fc8d0b51`
- 原始权重路径：`checkpoints/visnet_l6_h64_rbf64_r75.pth`

该文件是上游提供的冻结 ViSNet 权重，不是本项目训练得到的。PACER-FKG 只在其输出上进行四上下文分布交互分析。不得把权重归属或预训练贡献写成本项目成果。

验证：

```bash
sha256sum project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth
```

Windows PowerShell：

```powershell
Get-FileHash project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth -Algorithm SHA256
```

Geom2Vec 默认代码许可为 MIT；其 VAMPNet 衍生文件另有 GPL-3.0 标记。当前PACER-FKG不复制上游VAMPNet代码。许可说明见`LICENSE.geom2vec`和上游完整`LICENSE`。
