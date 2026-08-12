# skills-manifest.md

项目使用 `skills/` 目录存放本地技能副本。当前仓库**未内置本地 `skills/` 副本**,而是使用全局技能位置;后续若要为项目固化专属技能,按本清单维护。

## 技能目录结构(约定)

```
skills/
├── <skill-name>/          # 每个技能一个目录
│   └── SKILL.md           # 技能主文档
```

## 文件说明与 gitignore

| 路径 | 说明 | 是否 gitignore |
|------|------|----------------|
| `skills/` | 项目本地技能副本(全局位置在 `~/.config/opencode/skills/`) | **是**(gitignore) |
| `.opencode` | opencode 本地配置 | **是**(gitignore) |
| `docs/wiki/skills-index.md` | 技能索引(wiki 页面) | 否(提交) |
| `skills-manifest.md` | 本清单 | 否(提交) |

## gitignore 条目

```gitignore
# 技能目录(项目本地副本,全局位置在 ~/.config/opencode/skills/)
skills/

# opencode本地配置
.opencode
```

## 维护说明

- **新增技能**：复制/创建到 `skills/<name>/SKILL.md`,并更新 `docs/wiki/skills-index.md` 与 `docs/wiki/log.md`。
- **同步流程**：

```
项目本地 skills/  ←→  全局 ~/.config/opencode/skills/
       ↓
  docs/wiki/skills-index.md (索引)
  skills-manifest.md (清单)
  docs/wiki/log.md (变更日志)
```
