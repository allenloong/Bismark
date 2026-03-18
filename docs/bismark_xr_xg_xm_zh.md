# Bismark 中 XM / XR / XG 的定义（SE/PE）

本文根据 `bismark` 主脚本实现整理，重点对应：

- `extract_corresponding_genomic_sequence_single_end`
- `single_end_SAM_output`
- `extract_corresponding_genomic_sequence_paired_end`
- `paired_end_SAM_output`

## 1) XR 与 XG 的定义

- `XR`：read conversion state（读段转换状态），取值 `CT` 或 `GA`。
- `XG`：genome conversion state（参考转换状态），取值 `CT` 或 `GA`。

### Single-end (SE)

SE 在内部根据 index=0/1/2/3 映射：

- index 0: `XR=CT, XG=CT`（OT）
- index 1: `XR=CT, XG=GA`（OB）
- index 2: `XR=GA, XG=CT`（CTOT）
- index 3: `XR=GA, XG=GA`（CTOB）

并与输出方向（SAM FLAG 0x10）对应：

- FLAG 非 reverse（`is_reverse=False`）只能是 `(CT,CT)` 或 `(GA,GA)`
- FLAG reverse（`is_reverse=True`）只能是 `(CT,GA)` 或 `(GA,CT)`

### Paired-end (PE)

PE 的 `XG` 对两个 mate 相同；`XR` 分别写在 read1/read2 上。

四种有效组合：

- index 0: `R1 XR=CT, R2 XR=GA, XG=CT`（OT）
- index 1: `R1 XR=GA, R2 XR=CT, XG=GA`（CTOB）
- index 2: `R1 XR=GA, R2 XR=CT, XG=CT`（CTOT）
- index 3: `R1 XR=CT, R2 XR=GA, XG=GA`（OB）

## 2) XM 的定义

`XM` 是 methylation call string（甲基化调用字符串），字符集：

- `Z/z`（CpG，甲基化/未甲基化）
- `X/x`（CHG）
- `H/h`（CHH）
- `U/u`（未知上下文）
- `.`（非C位点或不相关位点）

方向处理：

- 若该 read 在 Bismark 内部为 `strand '+'`，`XM` 直接写入。
- 若为 `strand '-'`，输出前会 `reverse`（仅反转，不互补）。

PE 时 read1/read2 分别按各自 strand 规则处理。

## 3) Python 实现

已在 `tools/bismark_x_tags.py` 实现：

- SE 规则检查（XR/XG 与方向一致性）
- PE 成对规则检查（read1/read2 的 XR 组合 + XG 一致性）
- XM 基础合法性检查（字符集、长度）

## 4) 测试脚本

`tests/test_bismark_x_tags.py` 提供：

- 单元测试（不依赖 BAM）
- 可选集成测试（读取你给的 Bismark SE/PE BAM）

示例：

```bash
python tests/test_bismark_x_tags.py --se-bam your_single_end.bam --pe-bam your_pair_end.bam
```

也可直接用主脚本：

```bash
python tools/bismark_x_tags.py your_single_end.bam --mode single
python tools/bismark_x_tags.py your_pair_end.bam --mode paired
```
