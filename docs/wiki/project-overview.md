---
title: SAIR (IGP24) 项目概览
source: session
created: 2026-08-12
tags: [overview]
status: current
last_confirmed: 2026-08-12
audience: internal
stage: current
schema_version: 2
---

## 一句话定义

SAIR 仓库是面向 SAIR Foundation 举办的 **IGP24 数论竞赛**的 Python 客户端,用于**生成并验证 24 次首一(monic)、Q 上不可约的整系数多项式**。

## 主线目标

- 通过多函数族穷举/搜索,生成满足 IGP24 约束的高价值候选多项式。
- 通过 SAIR 公开 API(`submit.py`)批量提交、轮询验证,追求获得分数。
- 把竞赛进度(标签、判别式、剩余签名)固化为可审计的数据与文档。

## 交付边界

- **Code**：`igp24/scripts/gen_explore.py`(候选生成)、`igp24/scripts/submit.py`(提交客户端)。
- **Data**：`igp24/data/`(LMFDB 基线、API 进度快照、候选/测试批次)。
- **Docs**：双语 README、MIT 许可证、本 wiki 知识系统。

## 目录结构

```
igp24/
├── scripts/
│   ├── gen_explore.py   # 候选多项式生成器(7 个函数族)
│   └── submit.py        # SAIR 公开 API 提交客户端
└── data/
    ├── lmfdb_baseline.csv   # LMFDB 参考判别式基线
    ├── labels_progress.json # API 按标签发现进度
    ├── remaining_pairs.json # 尚未发现的 label×r 组合
    ├── test_batch.txt       # 冒烟测试批次(文档示例)
    └── explore_batch*.txt   # 生成的候选批次
```

