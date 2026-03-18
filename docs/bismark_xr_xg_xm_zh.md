# Bismark 中 XM / XR / XG 的定义（SE/PE）与 XM 重建

本文基于 `bismark` 主脚本逻辑，说明如何在 Python 中**根据 BAM 比对结果 + XR/XG + 参考序列重建 XM**，并与 Bismark 输出逐条比对。

## XR / XG 定义

- `XR`：read conversion state，取值 `CT` 或 `GA`
- `XG`：genome conversion state，取值 `CT` 或 `GA`

SE 的 4 种组合：

- `(XR=CT, XG=CT)` -> OT
- `(XR=CT, XG=GA)` -> OB
- `(XR=GA, XG=CT)` -> CTOT
- `(XR=GA, XG=GA)` -> CTOB

PE 时两个 mate 的 `XG` 一致，`XR` 在 read1/read2 上分别记录。

## XM 定义

`XM` 为 methylation call string，字符：

- `Z/z`：CpG（甲基化/未甲基化）
- `X/x`：CHG
- `H/h`：CHH
- `U/u`：未知上下文（N/X 或边界不足）
- `.`：非目标位点或不相关位点

核心判定：

- `XR=CT`：关注参考上的 `C`，读段 `C`=甲基化，`T`=未甲基化。
- `XR=GA`：关注参考上的 `G`（对链 C 位点），读段 `G`=甲基化，`A`=未甲基化。
- CIGAR 中 `I/S`（插入/软剪切）视为 `.`（不参与 C 上下文甲基化判断）。

## Python 实现

见 `tools/bismark_x_tags.py`：

- `qpos_to_rpos_from_cigar`：从 CIGAR 构建 query 位置到参考位置映射。
- `generate_xm_from_alignment`：按 Bismark CT/GA 规则逐位生成 XM。
- `validate_bam_with_reference`：读取 BAM + FASTA，重建 XM 并和 BAM 中 XM 精确比对。

## 使用方法

```bash
python tools/bismark_x_tags.py your.bam --reference genome.fa
```

通过测试脚本（单元测试 + 可选 BAM 对照）：

```bash
python tests/test_bismark_x_tags.py
python tests/test_bismark_x_tags.py --reference genome.fa --se-bam se.bam --pe-bam pe.bam
```
