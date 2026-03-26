# 20260316 Subclone Phase4 IGV 截圖索引

本資料夾整理 `20260315_subclone_phase4_casestudy_01.md` 的 9 個代表性 TP/FP 位點之 IGV session 與 PNG 截圖。

## 產出摘要

- 觀察視窗：所有位點統一使用 `SNV ±1000 bp`
- 位點清單：[case_regions.tsv](case_regions.tsv)
- 全批次 batch：[igv_batch.txt](igv_batch.txt)
- 單獨補跑 FP-B2：[igv_batch_FP_B2_only.txt](igv_batch_FP_B2_only.txt)
- manifest：[igv_snapshot_manifest.tsv](igv_snapshot_manifest.tsv)
- session 目錄：[sessions/](sessions/)
- snapshot 目錄：[snapshots/](snapshots/)

## 執行備註

- 本批最早一輪 PNG 是以 `/home/liaoyoyo2001/igv.sh -b ...` 直接連桌面 X session 輸出。
- 9 loci 全部已成功輸出 PNG。
- 初次整批執行在最後一張 `FP_B2_chr9_75383880` 遇到 Java heap 累積問題，因此以單獨 batch 補跑完成。
- 後續已驗證更穩定的新路徑：`user-space xvfb + slim template + per-locus runner`。
- 正式操作手冊請看：[20260316_IGV自動截圖操作手冊_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260316_IGV自動截圖操作手冊_01.md)

## 後續較穩定路徑

- slim template：[template_hcc1395_subclone_slim.xml](template_hcc1395_subclone_slim.xml)
  - 將原始 template 從 `23 resources / 11 panels` 縮到 `9 resources / 3 panels`
  - 移除遠端 `Refseq Select` track、Dorado 對照、重複 BAM 與不必要的輔助 VCF
- slim session 輸出目錄：[/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions_slim](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions_slim)
- per-locus runner：[run_igv_snapshots_per_manifest.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_igv_snapshots_per_manifest.sh)
  - 每個 locus 啟動一個新的 IGV process，比單一 Java process 連續 `load session` 更穩
  - 目前腳本預設已改為 `--use-xvfb always`
- user-space xvfb installer：[install_xvfb_user_space.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/install_xvfb_user_space.sh)
  - 不需要 `sudo apt install`
  - 可把 `xvfb-run` 與 `Xvfb` 裝到使用者可寫目錄
- `xvfb-run` 的角色主要是提供獨立 display，不是解決 IGV 的 `loadSession` UI null-pointer；真正降低 OOM 風險的是 slim template 與 per-locus runner。

## xvfb 安裝狀態

- `xauth` 已安裝。
- `xvfb` 已確認可用 user-space 方式安裝，不必依賴 `sudo apt install`。
- 在代理環境中，`Xvfb` 需於 sandbox 外啟動；session / manifest 生成可留在 sandbox 內。
- 標準建議流程是：先安裝 user-space xvfb，再用 per-locus runner + slim template。

一次性安裝：

```bash
scripts/analysis/install_xvfb_user_space.sh \
  --target /home/liaoyoyo2001/.local/opt/xvfb-user
```

批次截圖：

```bash
scripts/analysis/run_igv_snapshots_per_manifest.sh \
  --manifest /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions_slim/igv_snapshot_manifest.tsv \
  --use-xvfb always \
  --xvfb-run /home/liaoyoyo2001/.local/opt/xvfb-user/usr/bin/xvfb-run
```

## Case 索引

| Case | 類型 | 變異 | PNG | Session | 輔助圖 | 觀察重點 |
|---|---|---|---|---|---|---|
| TP_01 | TP | `chr6:145444893 G>A` | [PNG](snapshots/TP_01_chr6_145444893.png) | [XML](sessions/TP_01_chr6_145444893.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/chr6_145444893_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/chr6_145444893_distance_heatmap.png) / [methylartist](../../20260315_subclone_phase4_casestudy_assets/chr6_145444893_methylartist.svg) | 看 `HP2` 與 `HP2-1` 是否呈 LOH-like 偏態，並確認 REF 是否仍保留約半數支持。 |
| TP_02 | TP | `chr4:70548355 G>A` | [PNG](snapshots/TP_02_chr4_70548355.png) | [XML](sessions/TP_02_chr4_70548355.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/chr4_70548355_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/chr4_70548355_distance_heatmap.png) / [methylartist](../../20260315_subclone_phase4_casestudy_assets/chr4_70548355_methylartist.svg) | 看大量 `HP0` 是否真的是 phase 資訊不足，而 ALT read 是否只在其中一個大群中形成較小次群。 |
| TP_03 | TP | `chr5:153209947 C>A` | [PNG](snapshots/TP_03_chr5_153209947.png) | [XML](sessions/TP_03_chr5_153209947.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/chr5_153209947_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/chr5_153209947_distance_heatmap.png) / [methylartist](../../20260315_subclone_phase4_casestudy_assets/chr5_153209947_methylartist.svg) | 看低甲基背景下 ALT 群是否仍有穩定局部 pattern，以及 ALT 是否集中在某個次群而非均勻散布。 |
| TP_04 | TP | `chr16:35118902 G>A` | [PNG](snapshots/TP_04_chr16_35118902.png) | [XML](sessions/TP_04_chr16_35118902.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/chr16_35118902_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/chr16_35118902_distance_heatmap.png) / [methylartist](../../20260315_subclone_phase4_casestudy_assets/chr16_35118902_methylartist.svg) | 看 ALT 是否對應較低甲基 read 群，並注意 `HP3` / `HP0` 是否聚在複雜區塊而非隨機散布。 |
| TP_05 | TP | `chr7:109185781 G>T` | [PNG](snapshots/TP_05_chr7_109185781.png) | [XML](sessions/TP_05_chr7_109185781.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/chr7_109185781_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/chr7_109185781_distance_heatmap.png) / [methylartist](../../20260315_subclone_phase4_casestudy_assets/chr7_109185781_methylartist.svg) | 看高讀數下 `HP` 與 `Allele` 是否共線，ALT read 內部是否還能再分兩層。 |
| FP_A1 | FP | `chr8:93565727 C>T` | [PNG](snapshots/FP_A1_chr8_93565727.png) | [XML](sessions/FP_A1_chr8_93565727.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/fp_chr8_93565727_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/fp_chr8_93565727_distance_heatmap.png) | 看 2 條 ALT 是否完全混在大群中，以及高 `HP0` 區域是否伴隨多重比對/雜訊跡象。 |
| FP_A2 | FP | `chr9:137953060 T>C` | [PNG](snapshots/FP_A2_chr9_137953060.png) | [XML](sessions/FP_A2_chr9_137953060.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/fp_chr9_137953060_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/fp_chr9_137953060_distance_heatmap.png) | 看 ALT 是否其實與 REF 沒有可見差異，並觀察 158 CpG 高維資料下的「視覺無效應但統計顯著」案例。 |
| FP_B1 | FP | `chr7:52087777 A>T` | [PNG](snapshots/FP_B1_chr7_52087777.png) | [XML](sessions/FP_B1_chr7_52087777.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/fp_chr7_52087777_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/fp_chr7_52087777_distance_heatmap.png) / [methylartist](../../20260315_subclone_phase4_casestudy_assets/fp_chr7_52087777_methylartist.svg) | 先看 `chr7:52087776` 附近是否有 SEQC2 HighConf INDEL adjacency，再看高 CramersV 是否其實由 HP 分層主導。 |
| FP_B2 | FP | `chr9:75383880 T>A` | [PNG](snapshots/FP_B2_chr9_75383880.png) | [XML](sessions/FP_B2_chr9_75383880.xml) | [cluster](../../20260315_subclone_phase4_casestudy_assets/fp_chr9_75383880_cluster_heatmap.png) / [distance](../../20260315_subclone_phase4_casestudy_assets/fp_chr9_75383880_distance_heatmap.png) / [methylartist](../../20260315_subclone_phase4_casestudy_assets/fp_chr9_75383880_methylartist.svg) | 重點看 `75383879-75383880` 是否呈現一致的 adjacent mismatch / MNP-like pattern，而非獨立 SNV。 |

## 建議觀察順序

1. 先看 `TP_04` 與 `TP_01`，建立「看起來像真實 allele/subclone」的參考。
2. 再看 `FP_A1` 與 `FP_A2`，對照低 VAF 與高維假顯著。
3. 最後看 `FP_B1` 與 `FP_B2`，聚焦 HP confounding、INDEL adjacency 與 MNP 機制。

## 若要重跑

整批重跑：

```bash
scripts/analysis/igv_snapshot_from_template.sh \
  --template-xml /big8_disk/liaoyoyo2001/IGV_session/template.xml \
  --regions /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/case_regions.tsv \
  --output-dir /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions \
  --padding 1000
```

直接用 IGV batch 產生 PNG：

```bash
/home/liaoyoyo2001/igv.sh -b \
  /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/igv_batch.txt
```
