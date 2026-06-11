<!--
建立時間: 2026-06-09
狀態: proposal (位點甲基分群標籤 catalog — schema 提案，待用戶逐欄確認)
報告類型: paper_focus_catalog_schema_proposal
受眾: 廖子游（逐欄確認 schema 是否合理/有無缺漏）
provenance_note: 欄位定義為設計提案（🔵）；範例位點數字沿用已驗證集合（🟢/🟡）；閾值多為待校準（標 ⚠TBD）。
-->
<!-- provenance-verified: schema 為方法設計提案；範例 chr17/BRCA2 數字沿用 wf_a8ccbb34-3f7 對賬集合（🟡 待原檔）。 -->

# 位點甲基分群標籤 catalog — schema 提案（DEC-3，請逐欄確認）

> **L0 一眼結論**：catalog = **一位點一列**，把每個位點貼上「**分群明不明顯 / 是什麼性質（cis/copy/germline/artifact）**」的標籤 → 直接變論文 Results R6 + 回答你「哪些位點甲基分群明顯」。本檔是 schema **提案**，**請你逐欄確認定義/閾值是否合理、有無缺漏**，確認後才動手組表（CODE-catalog）。
>
> **L1 重點邏輯**：
> ① **核心欄 = `clustering_reliability`**（reliable/latent/none）—— 這欄直接答「哪些位點甲基分群明顯」。
> ② **每個位點最後收斂成一個 TAG**（7 類，見下）—— 讓零散觀察變成可發表的分類。
> ③ **多數閾值還沒校準（⚠TBD）+ 有已知前置依賴（5mC/5hmC dup-bug、cis 未接 SEQC2 CN）** —— 這些就是「缺漏」，要你看過再定。
> ④ **守 characterization≠filter**：catalog 是描述目錄，TAG-D（FP-prone）只標註不當 filter。

---

## L2 — 欄位提案（每欄：定義 / 閾值 / 來源 / 為何要 / 缺漏自評）

| # | 欄位 | 值 | 定義 | 閾值 | 來源 | 為何要這欄 | 缺漏自評 |
|---|------|----|------|------|------|-----------|---------|
| 1 | `locus_id` | chrom:pos (+gene) | 位點座標 | — | VCF | 主鍵 | gene annotation 來源待定 |
| 2 | `somatic_status` | somatic / germline-het / FP | 變異性質 | — | ClairS-TO VCF + SEQC2 truth | 區分 somatic vs germline（ALLELE-axis confound 根源）| 只 HCC1395 有 truth；single-pipeline |
| 3 | **`clustering_reliability`** ⭐核心 | reliable / latent / none | 甲基 read 分群可靠度 | reliable=CramersV Cochran 最小期望格≥5 且顯著；**latent=CV gated=0 但 PERMANOVA p<⚠TBD**；none=皆不顯著 | ISM RegionProcessor CramersV + PERMANOVA | **直接答「哪些位點甲基分群明顯」** | ⚠ PERMANOVA p 閾值 + minN 待定；latent(487) minN<15 FP-leaning，**只描述不可 filter** |
| 4 | `dbeta_magnitude` | large / medium / small | per-position allele 甲基平均差 | ⚠TBD：large≥0.2 / med 0.1–0.2 / small<0.1（待校準）| ISM-3 合併後的 Δβ 欄（**依賴 ISM-3**）| 抓「分群不明顯但平均差大」+ 與 clustering 交叉 | ⚠ 閾值需校準；高-Δβ regime FP-enriched（regression-to-extreme）→ 必配 clustering 解讀 |
| 5 | `clustering_origin` | somatic-linked / germline-allelic / drift | 分群來源 | blind-ARI vs germline-het null（imprinting 正控 GNAS/RB1=1.000）| ARI ruler | **D6：分出來是 somatic 還是 germline 背景** | somatic-specific 必跑 germline-het null（B2 證多為 germline-allelic）|
| 6 | `cis_status` | cis / copy / drift / untestable | cis 判定 | within d > HP-axis d 且 perm p<⚠TBD → cis；d_copy 主導 → copy（=subclone/copy-partition，HP1-1=subclone tag 非 CN-dosage）；pure-ALT/CGI-desert → untestable | normal-anchored cis-test + subclone/copy-partition（**依賴 D3**）| D3：真 cis vs subclone/copy-confounded | ⚠ 未 anchor SEQC2 CN；6 個 untestable 無解（除非 COLO829）|
| 7 | `axis` | HP-axis / ALLELE-axis | 哪個軸（somatic-controlled vs confounded）| — | ISM | ALLELE-axis 被 germline 基線 confound；LOH 位點只有 ALLELE-axis | LOH 位點全 ALLELE-axis → 須 germline-het null |
| 8 | `mod_type` | 5mC / 5hmC / mixed | 修飾類型 | — | BAM MM/ML | **避免 5mC+5hmC 雙列 collapse artifact**（BRCA2 buggy 教訓）| ⚠ MSA Level1 dup-bug 須先修否則此欄不可信 |
| 9 | `recurrence` | single(HCC1395) / cross(N樣本) | 跨樣本復現 | — | cross_sample | phenomenon vs private | 0/38 underpowered（同-locus 復發稀有非 private）|
| 10 | `provenance_tier` | 🟢 / 🟡 | 數字是否原檔對賬 | — | 02 文件庫 | 引用前知道能不能信 | 🟡 須投稿前補 |

---

## L2 — 輸出標籤分類（每位點收斂成 1 個 TAG = 你要的「移至某種標籤分類」）

| TAG | 條件（欄位組合）| 意義 | 範例 |
|-----|----------------|------|------|
| **TAG-A 乾淨 cis somatic-ASM** | reliable + cis + somatic-linked + somatic | 寶物（罕見，真 ASM exemplar）| chr17/TBC1D16（🟡 perm p=0.001 待對賬）|
| **TAG-B subclone/copy-confounded** | cis_status=copy | 排除於 cis 主張 | BRCA2（🟡 subclone/copy 主導，HP1-1=subclone tag 非 copy、非 CN-dosage；focal cis d_within=−0.023 邊際；% 不 robust）|
| **TAG-C germline-allelic 背景** | clustering_origin=germline-allelic | 大多數 read-clustering 屬此（非 somatic）| 多數位點 |
| **TAG-D high-Δβ FP-prone artifact** | large-Δβ + (none/latent) + LOH + 低覆蓋 | regression-to-extreme，**只標註不當 filter** | strong-ASM FP（🟡）|
| **TAG-E latent 真分群** | latent（CV gated 但 PERMANOVA 顯著）| 487 類，納回提升 characterization 分辨 | minN<15 FP-leaning |
| **TAG-F untestable** | cis_status=untestable | 6 個最強訊號（pure-ALT/CGI-desert），需新法/COLO829 | 6 strongest |
| **TAG-G 無訊號** | none + small-Δβ | 背景 | 多數 |

> 這 7 個 TAG 就是 Results R6 的主表骨架：論文可報「N 個位點分 7 類，其中 TAG-A（乾淨 cis）僅 X 個、TAG-C/D（germline/artifact）佔多數」→ 誠實 characterization。

---

## 🔴 L1 — 缺漏自評（你要確認的「是否合理/有無缺漏」）

**已知缺漏/前置依賴（組表前要先處理）**：
1. ⚠ **閾值未校準**：`clustering_reliability` 的 PERMANOVA p + minN、`dbeta_magnitude` 的 large/med/small cutoff、`cis_status` perm p —— 須先用 imprinting 正控/已知例校準。
2. ⚠ **`dbeta_magnitude` 依賴 ISM-3**（合併 Δβ）才有料 → 排序：ISM-3 先。
3. ⚠ **`cis_status` 依賴 D3 + 未 anchor SEQC2 CN** → copy 判定還沒接真值。
4. ⚠ **`mod_type` 依賴 MSA Level1 dup-bug 修復**（否則 5mC+5hmC collapse 重演 BRCA2 artifact）。
5. ⚠ **single-pipeline**：所有 `somatic_status` 來自 ClairS-TO → 封頂 ⭐3。

**可能還缺的欄（請你補）**：
- `n_CpG` / `coverage` / `region_size`（影響 reliability，建議加）
- `LOH_status`（in/out LOH，影響 axis + artifact 判讀，建議加）
- `gene_annotation` / `is_TSG_promoter`（生物意義，建議加）
- （若有）`RNA_expression`（cis-silencing 佐證，optional）

---

## L1 — 請你確認 3 件（確認後我就動手組 catalog）
1. **10 欄 + 7 TAG 合理嗎？**（特別核心欄 `clustering_reliability` 的 reliable/latent/none 三態定義）
2. **要補的欄**（n_CpG/coverage/LOH_status/gene/expression）哪些要加？
3. **前置順序**接受嗎：ISM-3（Δβ 欄）→ 閾值校準 → 組表？（或你要先用現有欄跑骨架，Δβ/cis 之後補）

---

## 🔒 schema 定稿（2026-06-09 用戶確認）

- ✅ 10 欄 + 7 TAG 照原提案；核心欄 `clustering_reliability` 三態（reliable/latent/none）OK。
- ✅ **補欄**：`n_CpG` / `coverage` / `region_size`（影響 reliability）· `LOH_status`（in/out LOH）· `gene` / `is_TSG_promoter`（生物意義）。
- ✅ **新增（用戶要求「明確數值 + 觀察分佈圖」）**：見下 B + C —— catalog **不只分類**，每分類欄**配實際數值** + 一組**觀察分佈圖**（見樹也見林）。
- ✅ **組表順序 = 先跑骨架再補**：先用現有欄（clustering_reliability / origin / somatic_status / axis / LOH_status / recurrence / coverage / n_CpG / gene）跑**骨架表** → 之後補 `dbeta`（依 ISM-3）+ `cis_status`（依 D3）。

### B. 數值欄（每分類欄配實際數值 — 不只 bucket，可回溯）

| 分類欄 | 配的數值欄 |
|--------|-----------|
| clustering_reliability | `cramersV_value` · `permanova_p` · `min_expected_cell` |
| dbeta_magnitude | `dbeta_max`（actual 值）|
| clustering_origin | `ari_value` · `ari_null` · `ari_p` |
| cis_status | `within_d` · `hp_axis_d` · `cis_perm_p` |

→ 每個 TAG 都能回溯到原始數值（標 🟢🟡 tier），不是黑箱標籤。

### C. 觀察分佈圖（forest view → 存 `04_figures/`）

| 圖 | 內容 |
|----|------|
| P1 | Δβ 分佈（all loci + large/med/small cutoff 線）|
| P2 | CramersV 分佈（reliable 閾值 + latent 區標色）|
| P3 | ARI 分佈（somatic TP vs germline-het null vs imprinting 正控 0.758/1.000）|
| P4 | within_d vs HP_axis_d scatter（cis vs copy 對角線分界）|
| P5 | **per-TAG 計數長條（7 TAG 各幾個 = Results R6 主圖）** |
| P6 | coverage / n_CpG vs reliability（檢查 reliability 是否被覆蓋驅動）|

> 最終 catalog ≈ **16 欄（9 分類 + 7 數值）+ 6 張分佈圖**；骨架先跑（不含 dbeta/cis），P5 per-TAG 長條 = R6 主圖。
> 互動逐欄確認版：`05_html_staging/catalog_schema確認.standalone.html`。
