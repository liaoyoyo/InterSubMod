---
title: 5-version (baseline / V2b / V3F / V5 / V6) IGV snapshots — session note
date: 2026-05-12
source_session: InterSubMod/research/igv_sessions/v6_all_5versions.xml
base_session: InterSubMod/research/igv_sessions/v5_all_4versions.xml
---

# 5-Version IGV Snapshots — Session Note

## Source Session

- **產出 session**：`InterSubMod/research/igv_sessions/v6_all_5versions.xml`
- **基底 session**：`InterSubMod/research/igv_sessions/v5_all_4versions.xml`（保留不動，純 V5 並列版）
- **加入方式**：在 V5 Resource line 後加入 V6 Resource；在 V5 Panel block 後加入 V6 Panel block（不取代既有 panel）

## 5 Versions in Stack Order

| Panel | BAM path |
|-------|----------|
| baseline | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam` |
| V2b | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_tagged.bam` |
| V3F | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam` |
| V5 | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam` |
| **V6** (new) | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam` |
| tumor (untagged) | `/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam` |
| normal | `/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam` |

## 6 PNG 含義

### Set A — V5 鐵證個案（V5 翻 HP2 對齊 paired）

| 檔名 | 位點 | 預期觀察 |
|------|------|---------|
| `SP1_chr19_17565944_5versions.png` | chr19:17,565,444-17,566,444 | V5 與 baseline/V2b/V3F HP 方向相反；V6 在 germline-absent 是否 revert 回 baseline 方向 |
| `SP2_chr19_12452332_5versions.png` | chr19:12,451,832-12,452,832 | 同上 |
| `SP3_chr19_12467180_5versions.png` | chr19:12,466,680-12,467,180 | 同上 |

### Set B — V6 鐵證個案（chr19 germline-absent revert）

| 檔名 | 位點 | 預期觀察 |
|------|------|---------|
| `V6win_chr19_52081584_5versions.png` | chr19:52,081,084-52,082,084 | V6 在 germline 缺席區，HP 分布回到 baseline 4.19:1 偏 HP1 |
| `V6win_chr19_55347952_5versions.png` | chr19:55,347,452-55,348,452 | 同上 |
| `V6win_chr19_8349597_5versions.png` | chr19:8,349,097-8,350,097 | 同上 |

## 重生流程

```bash
# 1. 確認 V6 BAM 存在
ls -lh /big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam*

# 2. 跑 IGV headless batch
DISPLAY=:20 /big7_disk/liaoyoyo2001/IGV_Linux_2.19.3/igv.sh \
  --batch /tmp/v6_5versions_batch.txt

# batch 內容備份：見本 note 同層 _batch_template.txt（如後續需要保留）
```

Batch 關鍵 IGV 指令：
- `colorBy TAG HP` + `group TAG HP`
- `preference SAM.SHADE_ALIGNMENTS_BY HAPLOTYPE`
- `maxPanelHeight 400`
- 每位點 `goto <locus>` → `sort base` → `snapshot <name>.png`

## 驗證紀錄

- xmllint: exit 0
- V5 path count: 4（Resource + Coverage + Junctions + Alignment id）
- V6 path count: 4（同上）
- 6 PNG 全 > 100K，dim 1150 × {1887–2435}（高度反映 7 個 align panel）
- SP1 visual 確認：5 個 tagged panel + tumor + normal 並列，HP TAG color 套用正確

## PPT 使用建議

- **Slide 04a/04b（V5 鐵證）**：SP1/2/3 三張展示 V5 翻 HP2 對齊 paired；V6 在這些位點是否保持 V5 phasing 或 revert，是 V6 設計效果關鍵
- **Slide 16（V6 germline-absent revert）**：V6win 三張為 V6 鐵證 — V6 在 germline 缺席區回到 baseline 行為（4.19:1 偏 HP1），與 V5 的繼承 phasing 對比
