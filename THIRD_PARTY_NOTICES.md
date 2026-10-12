# Third-party notices

本文件说明公开复现包依赖或引用的第三方组件。根目录 MIT 许可证只覆盖本项目自行编写、且有权许可的代码，不替代第三方条款。

## DrugCLIP

- 上游项目：DrugCLIP official implementation
- 本仓库位置：`project/tools/DrugCLIP-official/`
- 源代码许可证：Apache License 2.0
- 模型权重及模型输出：Creative Commons Attribution-NonCommercial 4.0 International（CC BY-NC 4.0）
- 本项目用途：受体—配体口袋兼容性迁移评分。该分数不等同于 PAM 功能预测。
- 许可证副本：[`LICENSES/DrugCLIP-LICENSE.txt`](LICENSES/DrugCLIP-LICENSE.txt)

## Geom2Vec 及 GraphVAMPnets 来源文件

- 本仓库位置：`submission/project/tools/geom2vec_checkpoints/`
- 大部分 Geom2Vec 代码：MIT License
- 明示源自 GraphVAMPnets 的文件：GNU GPL v3
- 许可证摘要副本：[`LICENSES/Geom2Vec-LICENSE.txt`](LICENSES/Geom2Vec-LICENSE.txt)

## Schrödinger Glide

Glide 是商业软件。公开包只保留申报结果中允许公开的派生数值，不提供 Glide 程序、许可证、商业数据库或用于替代授权的实现。冻结重放不会执行 Glide。

## 数据与模型责任

使用者应分别核对数据源、模型权重和软件依赖的原始条款。尤其是含 DrugCLIP 权重或输出的使用必须遵守非商业限制；本仓库不能向第三方授予超出上游许可的权利。
