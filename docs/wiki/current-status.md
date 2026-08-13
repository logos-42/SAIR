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

- 2026-08-13 分桶分析颠覆结论: 开放对 65% 的 2-部分 ∈ [2^10,2^14] (3,177 对), 73% 只要 3-部分 3^1 → 瓶颈是 2-部分不是 3-部分 (之前以为要 3^6)。最大桶 2^12·3^1 (695 对) / 2^14·3^1 (629 对)
- 2026-08-13 发现历史 bug: **PARI Pol() 收降幂, Vecrev 产升幂** — helix 生成器所有塔元素被"反转" (系数方向反), 导致群结构与设计不符。修正: Pol(bcoeffs[::-1]) / Polrev
- 2026-08-13 升维塔 B (gen_liftB.py): Q(ζ35)^+(√b) — 全实 abelian degree-12 基域 (2cos 型, 根∈(-2,2)) + √ 层共轭爆炸 → 闭包 2-部分 ≈ 2^14。b=低次多项式+常数偏移控制 r: 偏移 8-12→r=24, 0-3→r=16/20, -12→-8→r=0。系数 10^10-10^12 (vs helix2 的 10^76)。首批 176 个 (r: 0×51/16×11/20×7/24×82) 已提交 batch igp24_batch_9800bc...
- 2026-08-12 双语 README(英文主 / 中文副)与 MIT 许可证。
- 2026-08-12 接入「维基 llm」知识系统(本 wiki)。
- 2026-08-12 首个 root commit,分支统一为 `master` 推送到公开仓。

## 待办 / 未支持

- 更多函数族(稀疏度更高、特殊群构造)。
- 提交结果的自动化回写与统计。
- raw 源批量的 `ingest_raw` 登记。

