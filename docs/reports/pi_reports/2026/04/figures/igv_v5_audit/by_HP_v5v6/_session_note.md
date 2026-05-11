<!--
建立時間: 2026-05-12 01:35
最新更新: 2026-05-12
目標: 說明 by_HP_v5v6/ PNG 來源 session — V5 + V6 並列 audit snapshot
關聯檔案:
  - InterSubMod/research/igv_sessions/v5_v6_compare_with_paired.xml
  - InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6/_session_note.md
-->

# by_HP_v5v6 PNG 來源說明

## 目的

V5（pononly_v5_somatic_fallback）與 V6（v6_germline_absent_revert）**並列** audit snapshot，用於直接視覺比較 V6 germline-absent revert 對 V5 在三個 self-phasing 位點的 BAM tag 差異。

**前一版錯誤**：`v6_germline_absent_audit.xml` 用 `sed s|V5|V6|g` 整檔取代 V5 path 為 V6，V5 track 消失（無法比對）；該檔已 archive 為 `_archive/v6_germline_absent_audit_BROKEN_replaced_v5.xml`。新檔 `v5_v6_compare_with_paired.xml` 改為**加入** V6 panel 而非取代 V5。

## 來源 session

`InterSubMod/research/igv_sessions/v5_v6_compare_with_paired.xml`

含 7 BAM panels（panel order）：
1. PA_normal (paired germline)
2. PA_tumor (paired ground truth, 0.93)
3. BL_93 (longphase-to baseline, 0.93)
4. **V5_93** (pononly_v5_somatic_fallback)
5. **V6** (v6_germline_absent_revert) ← 新增
6. BL_06 (purity 0.6 baseline)
7. V5_06 (purity 0.6 V5)

加上既有 audit context（7 VCF + 6 BED：phase context / TP-FP / audit markers / LOH-GE）。

## 重生流程

```bash
IGV=/big7_disk/liaoyoyo2001/IGV_Linux_2.19.3/igv.sh
SESSION=/big7_disk/liaoyoyo2001/InterSubMod/research/igv_sessions/v5_v6_compare_with_paired.xml
OUT_DIR=/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v5v6

cat > /tmp/v5v6_audit_batch.txt <<EOF
load $SESSION
maxPanelHeight 400
snapshotDirectory $OUT_DIR
colorBy TAG HP
group TAG HP
preference SAM.SHADE_ALIGNMENTS_BY HAPLOTYPE

goto chr19:17565444-17566444
sort base
snapshot D_SP1_chr19_17565944_v5v6_audit.png

goto chr19:12451832-12452832
sort base
snapshot D_SP2_chr19_12452332_v5v6_audit.png

goto chr19:12466680-12467680
sort base
snapshot D_SP3_chr19_12467180_v5v6_audit.png

exit
EOF

DISPLAY=:20 $IGV --batch /tmp/v5v6_audit_batch.txt
```

## 檔案清單

| 檔案 | 尺寸 | 大小 | 位點 |
|---|---|---|---|
| `D_SP1_chr19_17565944_v5v6_audit.png` | 1150×4108 | 237 KB | chr19:17,565,944 |
| `D_SP2_chr19_12452332_v5v6_audit.png` | 1150×4098 | 254 KB | chr19:12,452,332 |
| `D_SP3_chr19_12467180_v5v6_audit.png` | 1150×3996 | 233 KB | chr19:12,467,180 |

## 與舊 by_HP_v6/*_audit.png 對照

| 項目 | by_HP_v6 (broken) | by_HP_v5v6 (corrected) |
|---|---|---|
| BAM panels | 6（V5 已被 V6 取代，缺 V5） | 7（V5 + V6 並列） |
| PNG 高度 | ~3523-3579 px | ~3996-4108 px |
| 高度差 | — | +473~530 px ≈ 1 panel |
| 適用情境 | 不可用（V5 缺席） | PPT 比對用 |

## Track 命名警告

BAM 檔名都是 `tumor_tagged.bam` — 視覺辨識依**載入順序**判斷（見 panel 順序 1-7）；PPT slide caption 必須明確標序。

## V6 在三 SP 的 HP 量化（已驗證，見 by_HP_v6/_session_note.md）

| 位點 | HP1 group | HP2 group | 偏向 |
|---|---|---|---|
| SP1 chr19:17,565,944 | 112 | 3 | HP1 |
| SP2 chr19:12,452,332 | 108 | 2 | HP1 |
| SP3 chr19:12,467,180 | 105 | 0 | HP1 |

V6 在三位點仍偏 HP1（與 baseline 同方向），確認 V6 germline-absent revert 並未修正 priority bug。
