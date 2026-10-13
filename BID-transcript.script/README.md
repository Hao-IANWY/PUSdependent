# BID-transcript.script

BID-seq 转录组分析的脚本合集：Snakemake 工作流 + 流程内工具 + 独立下游工具。

## 目录结构

```
BID-transcript.script/
├── README.md                    ← 本文件（总览）
├── CONFIG_GUIDE.md              ← config.yaml 编写说明（单/双端、UMI、去重等）
├── config.example.yaml          ← 配置模板，新分析复制后修改
├── workflows/                   ← Snakemake 工作流（snakemake -s 指定）
├── bin/                         ← 流程内工具，即 config 里的 script_dir
└── downstream/                  ← 独立下游工具，不被任何 smk 调用，手动运行
```

各子目录的内容与用法见其自带文档：

- [workflows/README.md](workflows/README.md) — 两个 smk 的定位、运行方式、表达 count 输出
- [bin/README.md](bin/README.md) — 全部流程内工具清单与分组
- [downstream/README.md](downstream/README.md) — 下游工具的详细用法（dependent/temperature 位点判定、坐标换算等）

## 调用关系

- 工作流内所有工具均通过 config 的 `paths.script_dir` 引用（形如 `{SCRIPT_DIR}/samFilter`），
  **script_dir 必须指向 `bin/` 目录**；R 脚本由独立的 `paths.pseudo_callsites_r` 指定。
- 因此 `bin/` 是整体部署单元：同步到生产拷贝时整个 `bin/` 一起拷，不要拆开。
- 2 个 smk 不硬编码任何工具的相对位置，移动本仓库位置只需改 config 里的路径。

## 运行工作流

```bash
snakemake -s workflows/polyA_mapping.smk --configfile config.yaml --cores 120
```

config.yaml 从 `config.example.yaml` 复制修改，字段说明见 `CONFIG_GUIDE.md`。
