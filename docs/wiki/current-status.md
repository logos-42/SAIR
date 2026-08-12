---
title: SAIR (IGP24) 当前状态
source: session
created: 2026-08-12
tags: [status]
status: current
last_confirmed: 2026-08-12
audience: internal
stage: current
schema_version: 2
---

## 已支持

- **候选生成** (`gen_explore.py`)：`tri` / `subst` / `rand` / `tr` / `comp` / `cheb` / `cyc` 7 个函数族,均过滤为「首一、常数项非零、Q 上不可约」,输出带实根个数注释。
- **批次提交** (`submit.py`)：读取 `SAIR_API_KEY` 环境变量,支持 `--slot` / `--dry` / `--poll`。
- **数据固化**：LMFDB 基线、API 进度与剩余签名快照已入库。

## 最近变更

- 2026-08-12 双语 README(英文主 / 中文副)与 MIT 许可证。
- 2026-08-12 接入「维基 llm」知识系统(本 wiki)。
- 2026-08-12 首个 root commit,分支统一为 `master` 推送到公开仓。

## 待办 / 未支持

- 更多函数族(稀疏度更高、特殊群构造)。
- 提交结果的自动化回写与统计。
- raw 源批量的 `ingest_raw` 登记。

