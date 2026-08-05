<!--
建立時間: 2026-07-23
目標: 獨立重算 HCC1395 strict k<=12 partition 的 block-level tree gate
處理範圍: HCC1395 chr1-22；production strict regions v2 + partition v2 + authoritative molecule sparse calls
關聯檔案: InterSubMod/research/20260723_k_gt12_chain_validity_audit/agents/audit_block_gate.py
-->

# HCC1395 block-level gate 獨立稽核

> **結論：現行 partition 會產生少數不連通 child blocks；`retained_molecule_weight=0` 不能當排除條件。** 服務 G1／G4／G5。

## 重算結果

| 檢核 | 數量 |
|---|---:|
| bounded blocks | 11,712 |
| `k>=2` | 11,700 |
| threshold-3 induced graph connected | 11,696 |
| threshold-3 induced graph disconnected | **4** |
| `k>=2` 且 retained segmentation weight=0 | **134** |
| 134 中經 molecule projection 後有 `MINREAD=3` local pattern | **75** |
| 134 中 pattern unsupported | **59** |
| 4 個 disconnected block 中 pattern unsupported | **4/4** |

134 個 zero-weight blocks 中另有 69 個至少具一個 multi-site supported pattern，53 個的 supported-pattern graph 可連通全部 block 位點。因此 zero weight 只表示沒有完整 constraint 被 partition 保留，不表示重新投影到 local block 後沒有可用 evidence。

四個 disconnected 例外：

- `chr6 U12b049833dea3ecf49966f5e:B0007`，k=8，分成 1／2／5 點。
- `chr6 U500892875397fd2b1986e107:B0001`，k=4，分成 3／1 點。
- `chr6 Uea95c3bd7e2a91b6fdd3fb32:B0007`，k=5，分成 3／2 點。
- `chr8 Uf209e406cd66a00dfb91598b:B0001`，k=2，分成 1／1 點。

## Adapter gate

現行 adapter 在 `exact_ps_partition_to_mlhp.py:470-480` 只排除：

1. `len(positions)<2`；
2. 沒有任何 pattern 達 `MINREAD`。

它**沒有**檢查 child block 的 threshold-3 induced graph 是否 connected，也沒有用 retained weight 作 gate。四個本次 disconnected blocks 恰好全部 pattern unsupported，因此在本批 HCC1395 仍會被第二道條件排除；這是資料上的巧合，不是 production contract。

既有 `exact_ps_mlhp_HCC1395_chr1_22.json` 綁定舊 `20260722...direct_big7_v2` 且 funnel 為 `exact_ps_membership_v1`，並非 `hcc1395_partition_v2`；目前沒有找到綁定 strict partition v2 的正式 adapter receipt，所以不可宣稱 strict downstream 已完成。

## 建議 gate

最小 production gate 應固定為：

```text
k_B >= 2
AND induced_threshold3_graph_connected
AND has_MINREAD3_supported_local_pattern
```

若不連通，先依 induced components 切開並重新投影 patterns；singleton component abstain。`retained_molecule_weight` 僅作 segmentation 診斷，不作 tree eligibility gate。

## 可重跑證據

輸入：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_partition_v2`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2`

命令：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  research/20260723_k_gt12_chain_validity_audit/agents/audit_block_gate.py \
  --partition-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_partition_v2 \
  --strict-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2 \
  --molecule-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2 \
  --min-read 3
```

實際輸出片段：`blocks_total=11712`、`k_ge_2=11700`、`connected=11696`、`disconnected=4`、`zero_weight_k_ge_2=134`、`zero_weight_pattern_supported=75`。

## Claim ceiling

本稽核證明 block 結構與 adapter 投影 gate 的工程行為；不證明任何 block 是真實 clone、唯一 topology 或 biological evolutionary unit。
