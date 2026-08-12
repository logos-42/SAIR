# SAIR

> 由 SAIR Foundation 主办的 **IGP24** 竞赛工具集。

本仓库包含 IGP24 竞赛的 Python 客户端:一个生成**24 次、首一(monic)、不可约**多项式的候选生成器,以及一个通过 SAIR 公开 API 验证这些多项式的提交客户端。

**English documentation: see [`README.md`](./README.md) · [English README](README.md)**

---

## 什么是 IGP24?

IGP24 是关于 **Q** 上 24 次多项式的**逆伽罗瓦问题(Inverse Galois Problem)**竞赛。目标是将伽罗瓦群实现为 24 次首一整系数多项式的分裂域。多项式以 **25 个逗号分隔的系数**提交(升幂排列,即常数项在前、首项系数在最后且为 1),每行一个多项式,允许以 `#` 结尾的可选注释。

## 目录结构

```
igp24/
├── scripts/
│   ├── gen_explore.py   # 候选多项式生成器(第 1 批:探索)
│   └── submit.py        # SAIR 公开 API 提交客户端
└── data/
    ├── lmfdb_baseline.csv  # 取自 LMFDB 的参考基线
    ├── labels_progress.json# 来自 API 的按标签发现进度
    ├── remaining_pairs.json# 尚未找到的 label/r 组合
    ├── test_batch.txt      # 冒烟测试批次(文档示例)
    └── explore_batch1.txt  # 生成的候选批次
```

## 依赖

- Python 3.8+
- [`numpy`](https://numpy.org/) — 通过 `np.roots` 统计实根个数
- [`sympy`](https://www.sympy.org/) — 不可约判定、因式分解、结式(resultant)

安装:

```bash
pip install numpy sympy
```

## 使用方法

### 1. 生成候选

```bash
python3 igp24/scripts/gen_explore.py > igp24/data/explore_batch1.txt
```

生成器会搜索多个函数族,所有候选都经过过滤,满足**首一、常数项非零、且在 Q 上不可约**:

| 函数族   | 说明                                            |
|----------|-------------------------------------------------|
| `tri`    | 三项式 `x^24 + a·x^k + b`                       |
| `subst`  | `g(x^k)`,`g` 为随机的不可约 deg 24/k 多项式     |
| `rand`   | 随机稠密首一次 24,小系数                       |
| `tr`     | 全实搜索(deg 2,3,4,6,8),全实根判定            |
| `comp`   | 全实对的和域(compositum),通过结式计算          |
| `cheb`   | 切比雪夫 `T24(x) - c`                           |
| `cyc`    | `x^24 ± a` 循环(cyclic)候选                    |

每行输出末尾带注释,标明函数族和实根个数,例如 `... # tri a=1 k=4 b=2 r=2`。

### 2. 提交批次

```bash
export SAIR_API_KEY=<你的密钥>
python3 igp24/scripts/submit.py igp24/data/test_batch.txt           # 快速检查几个
python3 igp24/scripts/submit.py igp24/data/explore_batch1.txt --poll # 提交并等待验证
```

主要选项:

| 参数       | 说明                                           |
|------------|------------------------------------------------|
| `--slot`   | 可选的提交槽位编号                             |
| `--dry`    | 只打印载荷大小,不真正发送                      |
| `--poll`   | 持续轮询批次直到验证完成                       |

提交的多项式必须**恰好有 25 个系数**;少于 25 个的行会被跳过并给出警告。`submit.py` 从环境变量 `SAIR_API_KEY` 读取密钥,绝不写入文件。

## 数据说明

- `lmfdb_baseline.csv` 枚举参考判别式(含每个标签的 `r` 与多项式判别式),用作评分基线。
- `labels_progress.json` / `remaining_pairs.json` 是从公开 API 拉取的快照,用于跟踪哪些 `label × r` 签名仍未被发现。

## 许可证

基于 [MIT License](./LICENSE) 开源。
