<!--
建立時間: 2026-05-12 01:00
最新更新: 2026-05-12
目標: 說明 by_HP_v6/ PNG 來源 session 與重生流程
關聯檔案:
  - docs/reports/pi_reports/2026/04/figures/igv_v5_audit/sessions/v6_germline_absent_audit.xml
  - docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6/igv_batch_v6.txt
-->

# by_HP_v6 PNG 來源說明

## 目前狀態（2026-05-12）

⚠ **本目錄內三張 PNG 為 fast-batch 產出**（用 `igv_batch_v6.txt` 4-BAM only），**未啟用完整 audit session**（缺 phase VCF / TP-FP marker / LOH BED）。

正式 V6 audit session 為：
`InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/sessions/v6_germline_absent_audit.xml`

該 session 從 `v5_purity_compare_with_paired.xml` 複製改 path 而來，含完整 audit infrastructure：
- 6 BAM (normal + tumor untagged + baseline + V6 + V6 purity sim)
- 7 VCF (TO candidate + 5 phase context + 2 TP/FP)
- 6 BED (5 audit markers + LOH/GE)

## 重生流程（用 V6 audit session 重新 snapshot）

```bash
# Step 1: 從 audit session 重跑 IGV headless
IGV=/big7_disk/liaoyoyo2001/IGV_Linux_2.19.3/igv.sh
SESSION=/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/sessions/v6_germline_absent_audit.xml

# 寫 batch script 載 session + goto 三位點 + snapshot
cat > /tmp/v6_audit_batch.txt <<EOF
load $SESSION
maxPanelHeight 400
colorBy TAG HP
group TAG HP
goto chr19:17565444-17566444
sort base
snapshot D_SP1_chr19_17565944_v6_audit.png
goto chr19:12451832-12452832
sort base
snapshot D_SP2_chr19_12452332_v6_audit.png
goto chr19:12466680-12467680
sort base
snapshot D_SP3_chr19_12467180_v6_audit.png
exit
EOF

DISPLAY=:20 $IGV --batch /tmp/v6_audit_batch.txt
```

## 檔案分類

| 檔案 | 類型 | 來源 batch | 用途 |
|---|---|---|---|
| `D_SP{1,2,3}_chr19_*_v6.png` | fast-batch (4-BAM) | `igv_batch_v6.txt` | 快速產圖；缺 audit context |
| `D_SP{1,2,3}_chr19_*_v6_audit.png`（待產出） | full audit session | `v6_germline_absent_audit.xml` | 正式 PPT / 報告用 |
| `igv_batch_v6.txt` | batch script | — | fast-batch reproducibility |

## V6 HP1:HP2 量化數據（已驗證）

| 位點 | HP1 | HP1-1 (`11`) | HP2 | HP2-1 (`33`) | group HP1 : HP2 |
|---|---|---|---|---|---|
| SP1 chr19:17,565,944 | 36 | 76 | 1 | 2 | **112 : 3** |
| SP2 chr19:12,452,332 | 32 | 76 | 1 | 1 | **108 : 2** |
| SP3 chr19:12,467,180 | 30 | 75 | 0 | 0 | **105 : 0** |

→ V6 在三 SP 位點**仍偏 HP1**（與 baseline 同方向），未修 priority bug
→ 揭示 V6 在 germline-absent 區域的 Layer 1.5 inheritance caveat

## Track 命名警告

baseline 與 V6 在 IGV 中的 track label 都顯示為 `tumor_tagged.bam`（BAM 檔名相同）— 視覺辨識需依**載入順序**：
1. normal
2. tumor (untagged)
3. baseline `/longphase-to-mod/output/baseline/tumor_tagged.bam`
4. V6 `/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam`

PPT slide caption 需明確標序。
