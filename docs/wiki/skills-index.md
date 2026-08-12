---
title: 技能索引
source: session
created: 2026-08-12
tags: [skills, index]
status: current
last_confirmed: 2026-08-12
audience: internal
stage: current
schema_version: 2
---

## 项目特定技能

当前仓库**未内置项目本地 `skills/` 副本**,使用全局技能位置 `~/.config/opencode/skills/`。

| 技能 | 用途 / 触发条件 | 位置 |
|------|-----------------|------|
| 维基 llm | 为项目搭建 compile-first 知识系统、wiki-first 规则、知识图谱 | 全局 |

## 通用技能(全局)

全局技能位于 `~/.config/opencode/skills/`,由各模型(Claude Code / Codex / Cursor / Windsurf)按需加载,不在本仓库内。

## 技能位置说明

- **项目本地**：`skills/`(已 gitignore,不提交)
- **全局**：`~/.config/opencode/skills/`(不在此仓库内,无需 gitignore)
- **wiki 索引**：本页(`docs/wiki/skills-index.md`),随仓库提交
- **清单**：`skills-manifest.md`,随仓库提交

## 维护流程

详见 [`skills-manifest.md`](../../skills-manifest.md)。
