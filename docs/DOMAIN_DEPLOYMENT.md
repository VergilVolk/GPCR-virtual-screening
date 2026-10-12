# 域名与部署方案

## 推荐名称

首选：`pacerm4.org`

- 短、可读，直接对应平台名称和靶点；
- `.org` 适合科研与开源项目，不暗示商业药物服务；
- 2026-10-12 通过 PIR 官方 RDAP 查询未发现已注册记录。可用性会随时变化，购买前必须再次核验。

备选顺序：`pacer-m4.org`、`m4pacer.org`、`pacerplatform.org`。不建议使用免费二级域名、动态 DNS 或带随机字符串的域名作为答辩主入口。

## 当前部署

网站由 `docs/` 生成，通过 `.github/workflows/pages.yml` 部署至 GitHub Pages。它是纯静态只读网站：

- 不上传分子或个人数据；
- 不在服务器侧运行 DrugCLIP、Glide 或分子动力学；
- 不需要维护本地电脑常开；
- 每次 `main` 分支更新后自动发布并获得 HTTPS。

## 购买后的绑定顺序

1. 在正规注册商购买域名并启用账户双因素认证、注册锁和自动续费。
2. 在 GitHub 账户中先验证域名所有权，再到仓库 Pages 设置添加自定义域名。
3. DNS 为 `www` 添加指向 `vergilvolk.github.io` 的 CNAME；根域按 GitHub 官方说明添加 A/AAAA 或 ALIAS/ANAME。
4. 等待解析后启用 Enforce HTTPS，核对证书和跳转。
5. 不使用通配符 DNS，避免子域接管风险。

域名购买和账户实名/付款属于仓库之外的操作；在团队确认名称和付款主体前，不应代为注册。
