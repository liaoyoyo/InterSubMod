<!--
建立時間: 2026-07-15T04:25:00+08:00
目標: 保存單一 FP focal-ALT 多群假說的四路獨立 agent 稽核與交叉共識
處理範圍: provenance / statistics / topology / legacy compatibility
關聯檔案: InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/report_dataset.json
-->

# Multi-agent consensus: single-FP ALT-read methylation multigroup

## 審查設計

四位 agent 使用不同責任邊界獨立審查，主 agent 最後以共同輸入重新核算交集。這是 adversarial review，不以多數決取代原始數據。

| Reviewer | 審查邊界 | 主要輸入 | 核心質疑 |
|---|---|---|---|
| Kant | provenance / LongPhase-S tags | latest manifest、VCF、sidecar、BAM receipts | 是否真為同一次最新 PASS；BAM 是否覆寫；HP 詞義是否誤用 |
| Bernoulli | statistics / null / controls | FP/TP site results、null traces、REF control | 分群穩定是否受 null、coverage、shared reads、threshold 影響；是否 FP-specific |
| Linnaeus | topology / biological claim | layered topology、HP families、case rows | 單一 sSNV 能否識別 linear/branching；same-pipeline tree 是否正交 |
| Copernicus | historical compatibility | legacy/current reads、methylation、distance、endpoints | 舊版比例與圖片能否直接合併；舊 null 是否有 fail-open |

## Reviewer findings

### Kant: provenance 完整，但 HP tag 不能升級為 clone label

- 最新輸入由 `20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json` 鎖定，manifest SHA-256 為 `16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45`。
- 7 dataset 的分析輸入都是同一次 LongPhase-S recalibrated `FILTER=PASS`，truth label 只在下游 post hoc split。
- 最新 HP/PS 來自 `20260711_longphase_s_raw_all_production_sidecars_v2` exact coordinate sidecar；新 BAM 全落在 observation workspace，未覆寫原始或 canonical BAM。
- `885,047/885,047` unique output alignments 具有 MM+ML；sidecar missing 為 `0`。
- `phase_anchored_robust_epigenetic_candidate` 舊欄名會誤導。`1-1/2-1` 只表示 LongPhase-S 第一 somatic haplotype family，不能稱為已由 phase 確認的 clone。

### Bernoulli: 多群可重現，但不是 FP-specific，嚴格存活率很低

- forced silhouette 與 null-gated stable 是兩個非巢狀 endpoint：交集 `364`、forced-only `527`、stable-only `219`、Jaccard `0.328`；不可畫成 `891 -> 583` 的順序淘汰。
- stable `583/4,967 = 11.74%`；其中 `582/583` 通過 modal assignment pairwise minimum ARI `>=0.8`。因此 high-threshold cluster membership 本身穩定，但生物詮釋仍未成立。
- 發現率隨 ALT depth 從 `6.02%`（6-9 reads）升至 `19.23%`（>=40 reads），site-level prevalence 受 power 強烈影響。
- 1:1 matched TP 對照：高門檻 FP `103/7,745 = 1.330%`、TP `98/7,745 = 1.265%`，paired `p=0.772`；沒有 FP specificity。
- 6 biological samples 等權高門檻 risk difference 為 `-0.111` percentage points，方向 2 正、1 平、3 負。
- row-circular null、REF background 與 assignment gate 全部取交集只剩 `10` sites / `9` unique readsets；此結果未做 genome-wide FDR calibration，不能當 PPV。
- `583` stable sites 共享 reads 後只形成 `380` connected components；一般 Wilson CI 忽略此依存，偏樂觀。

### Linnaeus: 單一位點無法識別 linear evolution

- 一個 focal sSNV 最多形成 `ROOT -> ALT` 的 trivial edge，無第二 mutation state 就沒有 ancestor/descendant ordering，不能區分 linear 與 branching cellular ancestry。
- 103 個高門檻候選分布於 linear regional context `22`、branching/recurrence `25`、no tree `28`、incomplete `22`、trivial `6`。
- 相對其他 FP，高門檻候選未富集 linear regional context：Fisher OR `1.217`、`p=0.441`。
- focal single-site linear topology identifiable `0`；layered tree orthogonal subclone confirmation `0`。
- 10 個 strict follow-up sites 仍含 linear `3`、branching `1`、incomplete `3`、no tree `3`；strict methylation pattern 也不指定 topology。
- 尚缺 copy number、purity-adjusted CCF、多 sSNV read co-occurrence、single-cell 或獨立 lineage barcode。

### Copernicus: 舊圖可參考現象，舊率不可直接比較

- 最新與 legacy 共同 FP `3,302`；ALT readset `3,302/3,302` 相同，methylation matrix `3,302/3,302` 相同。
- Bernoulli distance matrix 僅 `2,081/3,302 = 63.0%` 完全相同，反映 distance/missingness contract 改變。
- 共同位點 stable endpoint：latest `372`、legacy `670`、Jaccard `0.548`；高門檻 Jaccard `0.331`。
- 舊 null trace 曾容許 `n_valid_null=0` 的節點通過；legacy stable 中 `374/687` sites 至少有此類 pass，屬 fail-open 風險。
- 因此 prior InterSubMod 圖可作同一 read-level pattern 的歷史視覺證據；legacy `21.55%` 不可當 current prevalence 或與 current `11.74%` 合併。

## 共識與 claim tier

| Claim | 共識 | 證據層級 | 可用敘述 |
|---|---|---|---|
| focal ALT reads 存在甲基化多群 | 支持 | L2：本 pipeline 直接量測 + null/seed sensitivity | stable focal-ALT methylation multigroup |
| 少數位點值得後續驗證 | 支持 | L2：strict sensitivity，不含 FDR/PPV | strict epigenetic follow-up candidate |
| truth FP 多群比 TP 特異 | 不支持 | L3 反證：matched TP 高門檻差異不顯著 | 不可宣稱 FP specificity |
| 多群就是 high-probability subclone | 不支持 | L4 假說；缺正交 cellular identity | hypothesis-generating only |
| 單一 sSNV 支持 linear evolution | 不可識別 | L5：結構上缺第二 mutation state | one-site topology unidentifiable |

四位 reviewer 無實質結論分歧：都同意保留「甲基異質性候選」價值，但拒絕直接升級為 subclone 或 linear ancestry。

## 外部審查完成

外部 Claude Code 已用 read-only 設定完成終版稽核，15 項獨立重算全部通過；判定為 `VALID WITH CORRECTIONS`，且未發現會翻轉主結論的數據錯誤。三項呈現修正（both-evaluable stable 對稱揭露、legacy 分子/分母、ALT readset 與 `reads.tsv` byte identity 分離）均已納入主報告。完整輸出見 `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_external_claude_code_review.md`。
