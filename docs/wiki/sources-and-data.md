---
title: 资料与数据
source: session
created: 2026-08-12
tags: [data, raw]
status: current
last_confirmed: 2026-08-12
audience: internal
stage: current
schema_version: 2
---

## 已入库数据(随 Git 提交)

| 路径 | 说明 | 来源 |
|------|------|------|
| `igp24/data/lmfdb_baseline.csv` | 参考判别式(每标签 r、多项式判别式) | LMFDB |
| `igp24/data/labels_progress.json` | 按标签的发现进度/剩余签名 | SAIR 公开 API |
| `igp24/data/remaining_pairs.json` | 尚未发现的 `label×r` 组合 | SAIR 公开 API |
| `igp24/data/test_batch.txt` | 冒烟测试批次(文档示例) | 手写 |
| `igp24/data/explore_batch*.txt` | 生成的候选批次 | `gen_explore.py` |

## raw 根目录

原始资料默认放在本地 raw 根目录,不直接进 Git:

```text
../sair_raw/
```

GitHub 里只保留 manifest 与编译结果。新 raw 一多,直接跑：

```bash
python3 scripts/ingest_raw.py
python3 scripts/stale_report.py
python3 scripts/delta_compile.py --write-drafts
```

