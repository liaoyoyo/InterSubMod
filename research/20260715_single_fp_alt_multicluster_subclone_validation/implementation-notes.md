<!--
建立時間: 2026-07-15T01:46:40+08:00
目標: 保存單一 FP focal-ALT 多群全量驗證的進行中設計決定、偏離、折衷與未決事項
處理範圍: 7 datasets / 6 biological samples；最新 LongPhase-S PASS + 既有 InterSubMod complete matrices
關聯檔案: InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/pre-decision-audit.md
-->

# Implementation Notes: FP focal-ALT multicluster validation

## 設計決定

- 主單位固定為 `(dataset, CHROM, POS, REF, ALT)`；HCC1395 與 DORADO 分列 dataset、合併 biological sensitivity 另列。
- `focal ALT` 僅接受 `reads.tsv:alt_support=ALT` 或由 raw BAM 在焦點 base 直接重算；`HP=1-1/2-1` 只作 carrier/HP annotation，不可代替 focal ALT。
- 舊 silhouette `best_k>=2` 僅作 legacy sensitivity；主 verdict 使用 phylo-v4.1 的 per-CpG permutation-null、minimum group size 3 與 modal stability。
- 結果分五層：輸入 FP → 有 ALT read → 可評估（>=6 ALT reads）→ null-validated multigroup → 排除主要混雜後 residual candidate；只有再具獨立 genetic evidence 才可升 subclone-supported。
- 「linear」只在多 sSNV mutation-state 共現可比較 ancestor/descendant 時討論；單一 sSNV 不作 linear/branching 判別。

## 偏離

- 最新 LongPhase-S production 不保存完整 tagged BAM，而保存 truth-unaware recalibrated PASS VCF + exact HP/PS sidecar；因此不能把 canonical 舊 tagged BAM 宣稱為最新 tag payload。
- 最新 layered v5 的 L3 methylation 為 `not_evaluated`；本輪需 bounded materialization 或 raw-tag direct extraction 才能產生最新 focal-ALT methyl evidence。

## 折衷

- 先以全 3,458 canonical LongPhase-S FP complete matrices做可重現歷史基線；最新 PASS primary subset/重算另列，兩套分母不混合。
- 對 FP<100 的 dataset 仍完整報數，但比例以 Wilson 95% CI 與低 N 警示呈現，不做單樣本 prevalence 強解讀。
- CN/LOH 只在有 measured source 的 dataset/case 做 annotation；不可用缺失樣本假設 neutral。

## 結案狀態

- 本輪分析、數據驗算、外部審查、standalone HTML、desktop/mobile browser QA、SoT 索引與 evidence ledger 終版 entry 均已完成；cycle 已關閉。

## 2026-07-15 已解決事項

- 最新 LongPhase-S recalibrated `FILTER=PASS` 經 post hoc truth split：FP `7,745`、TP `335,296`；與 legacy canonical FP 重疊 `3,302`，新增或先前未 materialize `4,443`。
- 7 個最新-tag bounded BAM 已完成：`885,047` 個 unique alignments、`885,047` 個 MM+ML 完整、`86,291` 個 duplicate identity 依預先規格收斂；輸出均在 observation workspace，未覆寫 canonical 或原始 BAM。
- 最新 FP InterSubMod per-site reads/methylation/distance `7,745/7,745` 完整；summary omission `119`，其中 `46` 可評估，但主分析直接讀 per-site files，不依賴 summary coverage。
- 主分析 frozen 結果：可評估 `4,967`；stable `583`；residual `453`；HP-tag-covered high-threshold `103`；orthogonal subclone confirmation `0`。
- 1:1 matched TP `7,745` 對已完成。高門檻 endpoint：FP `103/7,745`、TP `98/7,745`，paired `p=0.772`；無 FP specificity。
- 6 biological samples 等權後 stable risk difference 僅 `+0.113` percentage points；移除 HCC1395 + Dorado 後為 `-0.153` percentage points、`p=0.780`。
- REF 背景與 row-circular null 等 strict gates 交集後剩 `10` sites / `9` unique ALT readsets，即 `0.201%` evaluable FP 或 `0.129%` 全部 FP；仍只命名為 strict epigenetic follow-up candidates。
- 所有 focal ALT reads 共 `113,845`，LongPhase-S 明確第二 somatic haplotype tag (`1-2/2-2`) 為 `0`；HP tag 不構成獨立 clone 證據。
- topology：所有 endpoint 的 focal single-site linear identifiability 與 layered-tree subclone confirmation 都是 `0`；E4 linear regional-context enrichment Fisher OR `1.217`、`p=0.441`。
- legacy compatibility：共同 `3,302` 位點 ALT readset 與 methylation matrix 完全相同，Bernoulli distance matrix 僅 `63.0%` 相同；舊圖可作歷史敏感度，舊比例不可作新版估計。

## 最終敘述鎖定（外部 review 後）

- 可支持：部分 truth-labeled FP focal-ALT reads 呈 null-gated、跨 seed 穩定的甲基化異質性，少數通過更嚴格背景與 null 檢核，值得正交追蹤。
- 不可支持：把單一 FP 的 ALT-read 多群直接稱為高機率 subclone；把單一 sSNV 或 same-pipeline regional tree context稱為 linear evolution。
- 對外術語：`stable focal-ALT methylation multigroup`、`strict epigenetic follow-up candidate`；禁止使用 `confirmed subclone`、`linear subclone`。
- 外部 Claude Code 結論：`VALID WITH CORRECTIONS`；15 項重算通過，未發現會翻轉 verdict 的數據錯誤，三項呈現修正均已併入終稿。
