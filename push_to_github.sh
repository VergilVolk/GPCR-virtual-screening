#!/bin/bash
# 推送到 GitHub (仓库已存在, remote 已配置)
# 用法: bash push_to_github.sh "commit message"
set -e
MSG="${1:-update: CHRM4 PAM project}"
git add -A
git commit -m "$MSG" || echo "nothing to commit"
git push -u origin main
echo "done"
