<!--
建立時間: 2026-06-27
類型: L8 甲基輔助有效性審查 — per-region 漏斗 + power dose-response + cis-control 可評估性（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3）
data_sources: data/sm_methyl_sufficiency_audit.json, data/sm_methyl_reextract_ALL_perregion.tsv, data/sm_methyl_reextract_ALL.json
-->

# L8 — 甲基資訊對 clone/subclone 準確度的輔助有效性審查

> **問題（使用者）**：甲基的「分群結構 / 差異狀況 / 差異對齊位點」是否能**提升 clone/subclone 準確度**？作為**標記/輔助是否足夠有效**？
> **方法**：把 L4 的終點數字（6.6% corroborate）拆成**per-region 有效性漏斗** + **power dose-response** + **cis-control 可評估性**三軸。資料 = genome-wide BAM MM/ML 重抽（740 區，每區落 17 欄）。
> 圖：`figures/08_methyl_sufficiency.png`。資料：`data/sm_methyl_sufficiency_audit.json` + `data/sm_methyl_reextract_ALL_perregion.tsv`（740 行）。
> 對抗稽核：3 獨立 verifier 從 raw TSV 重算（matches=true ×3）+ evaluator 綜合（見 §6）。

## §0 一句話裁決（頂層）

🔴 **BOUNDED_AUXILIARY — 甲基不足以作為獨立或主要的 subclone 標記/偵測器。** 甲基差異**確實存在且強烈受覆蓋度閘控（power-gated）**：高功率區（popB_n≥20）corroboration 達 **54.9%**；但 (1) ONT 覆蓋度使**僅 51 區達此功率**，(2) 甲基**從未獨立偵測任何新 partition**（0；corroborate≠detect），(3) **cis-control 修正後 14/49 可評估**（str/int bug 修正，見 §4），**其中 11 解析為 germline 等位特異甲基化（HP1≠HP2，非 subclone）、僅 3 為 subclone/somatic-cis 候選** → where testable，甲基更像 **allele-specific methylation 偵測器**而非 **subclone 偵測器**。CN 分層（使用者校正）：**neutral 8/8 可評估、7 germline-cis**；**LOH 35/41 單倍型不可評估**（germline ASM 本就不可能）。→ 對應既有 **auxiliary-not-driver / ⭐3** 上限。**［L1 源碼/數據重現］**

## §1 三軸方法

| 軸 | 問的問題 | 量化 |
|---|---|---|
| **充分度漏斗** | 有多少位點「有足夠資料」能跑甲基分群？ | target→tested→powered→corroborated |
| **power dose-response** | 6.6% 偏低是「甲基沒訊號」還是「覆蓋不足」？ | corroboration rate vs popB_n bin |
| **cis-control 可評估性** | corroborated 真的是 subclone-specific 還是 cis-ASM？ | hp_control_eval（HP1 vs HP2 各≥3 reads）|

> **「corroborated」定義**：兩大基因型 population（如 RR vs AA）間，per-CpG MWU + BH-FDR，≥1 個 CpG 達 q<0.05 且 |Δβ|≥0.2。此為 **genotype-anchored**（非循環），非無監督 clustering。

## §2 充分度漏斗（verified）

| 階段 | 數量 | 比例 | 意義 |
|---|---|---|---|
| structured CN-clean target | **754** | — | 有 ≥2 read-population 的結構區 |
| − 不足 population（<2 geno-pop ≥3 reads）| −14 | | |
| **tested** | **740** | 100% | 進入甲基檢定 |
| 其中 **0 可測 CpG**（無 CpG 兩側各≥3 read）| 193 | 26.1% | 雖「tested」但實質無法測甲基 |
| **powered**（CpG≥10 且 popB_n≥5）| 459 | 62.0% | — |
| underpowered | 281 | 38.0% | corroborated=0（全部）|
| **corroborated**（≥1 sig CpG）| **49** | **6.6%** | — |

🔑 **覆蓋是瓶頸**：median 每區 reads=24、median 較小 population=**6 reads**（35.7% 的區 <5）；median 可測 CpG=185 但 **29.3% 的區 <10 CpG**。**［L1］**

## §3 🔑 Power dose-response — 訊號「存在但 power-gated」（對抗稽核修正）

> ⚠ **本節修正了初版誤判**：初版說「低 corroboration = 無訊號非缺檢定力」。對抗稽核（power-vs-signal verifier）揭露 popB_n≥5 門檻**過寬**遮蔽了 power 主效應 → 修正如下。

| popB_n（較小 population read 深度）| corroboration rate | n |
|---|---|---|
| 5–7 | **0.6%** | 1/156 |
| 8–11 | **1.5%** | 2/131 |
| 12–19 | **13.0%** | 18/138 |
| **≥20** | **54.9%** | **28/51** |
| strict high-power（CpG≥100 且 popB≥15）| 33.3% | 43/129 |

→ 🔑 **甲基差異訊號確實存在，且隨覆蓋度單調陡升**：高功率（popB_n≥20）達 **54.9%（overall 6.6% 的 8.3×）**。**6.6% 偏低主因是覆蓋不足**——**58.8% 的「powered」其實落在 popB_n[5,12) 的功率飢餓帶（rate~1%）**，且**僅 51 區達 popB_n≥20**。**［L1］**

> **意涵（雙面）**：① 正面 — 若覆蓋夠深，甲基**能**區分基因型 population（非雜訊）；② 限制 — HCC1395 ONT 覆蓋使**絕大多數結構區達不到該功率**，故「實際可用」的比例仍小。

## §4 🔴 cis-control 可評估性（str/int bug 修正 + CN 分層；最關鍵修正）

> **🔴 2026-06-29 BUG 修正（推翻先前 structural-zero 敘述）**：先前本節記「hp_control_eval = 0/740（structural zero、測試從未執行）」**是 str/int 型別 bug 假象**。腳本以 `reads[rn][2] in (1, 11)`（整數）比對 HP，但 pysam 對此 BAM 的 `HP:Z:` **字串** tag 回 str `'1'`/`'1-1'`，`'1' in (1,11)` 恆 False → **HP 從未被讀到**。BAM 實有大量 HP reads（chr17 抽樣 `HP:Z:1`×70,709 + `HP:Z:2`×112,070）。已修為字串比對（`sm_methyl_reextract.py:140-143`，兼容 longphase-to 整數）+ 全基因組重跑；corroboration 完美重現（754→740→49，證實核心 genotype-anchored 不受影響）。**［L1 pysam 型別實測 + 全基因組重算］**

**修正後真實 cis-control**（49 corroborated，`sm_methyl_reextract_ALL_perregion.tsv` 重算）：

| 指標 | bug 版（假象）| 修正後真值 |
|---|---|---|
| hp_control_eval=1（HP1{1,1-1} vs HP2{2,2-1} 各≥3 reads）| 0/49 | **14/49** |
| cis_explained（sig CpG 在 HP1 vs HP2 也分離 = germline 等位特異）| 0 | **11** |
| hp 可評估但 cis 未解釋（→ subclone/somatic-cis 候選）| — | **3** |
| hp 不可評估（單倍型/真 LOH）| 49 | 35 |

**按 CN 分層**（CN 來自 SEQC2 bed，不受 bug 影響）：

| CN | corroborated | hp 可評估 | cis_explained=germline-cis | 候選 |
|---|---|---|---|---|
| **neutral** | 8 | **8/8** | **7** | 1 |
| **LOH** | 41 | 6/41 | 4 | 2 |

→ 🔑 **修正後重大結論**：
1. **HP 排解 cis-ASM 在此資料「可行」**（14/49 corroborated 可評估；先前「0 不可評估」是 bug）。
2. **neutral 區（germline ASM 可能）：8/8 可評估、7/8 解析為 germline-cis** → 該處甲基 corroboration **主要是 germline 等位特異甲基化、非 subclone**（**首次有正向 cis-control 證據**，取代先前的 structural-zero 推託）。
3. **LOH 區（germline ASM 結構不可能）：35/41 正確為單倍型不可評估**（與 LOH 一致）；6/41「可評估」+4「cis_explained」屬 **annotation/tagging 邊界效應**（真 LOH 不應有乾淨雙-HP germline read → segment-level CN 與 read-level 區不完全對齊 / `2-1` somatic-anchored read 混入），其 HP1≠HP2 **非 germline ASM**（LOH 不可能），仍需排 somatic-cis。
4. **僅 3 個區**（hp 可評估但 cis 未解釋）= 真正抗拒 germline-cis 解釋的 **subclone/somatic-cis 候選**。

> **修正方向（對裁決的影響）**：先前「cis-control structural-zero → 撤回 subclone-specific」改為 **「cis-control 可評估 14/49；11 解析為 germline-cis（甲基非 subclone）、3 候選」**。這對 **BOUNDED_AUXILIARY 裁決是更強的反向證據**：where testable，甲基 corroboration 多為 germline 等位特異甲基化而非 subclone — 即甲基更像「allele-specific methylation 的偵測器」而非「subclone 偵測器」。CN 分層（41 LOH / 8 neutral，使用者 LOH 校正）仍有效；唯 LOH 的「hp_control_eval=0」原句改為「35/41 單倍型不可評估、6 為邊界效應」。
> 與既有 memory 一致：`within_clean≠subclone`（Jaccard 0.123）、double-dip + somatic-cis confound 仍須排（即使 11 germline-cis 已解析，3 候選仍非「已證 subclone」）。

## §5 邊際效用 — 0 新 partition + 反向趨勢

- **new_partitions_detected_by_methylation = 0**：甲基從未獨立偵測任何基因型遺漏的 subclone。
  - ⚠ 部分屬 **by-construction**（corroboration 設計只測「對齊」不搜尋「新切割」）；偵測端的 null 另由專案獨立的 **tumor-only 無監督分群 NEGATIVE** 佐證（見 memory）。
- **by tree shape**（corroboration rate）：full_tree **2.6%（1/39）** < sibling 4.8% < linear_nested 6.8% < co_linked 9.7%。→ 🔴 **結構最複雜、subclone 最多的 full_tree 反而最低**，與「好的 subclone 標記」應有的方向相反。
- corroborated 中 **18.4%（9/49）僅靠單一 CpG**（最弱證據）。**［L1］**

## §6 對抗稽核裁決 + 對既有報告的修正

**3 獨立 verifier（從 raw TSV 重算）+ evaluator 綜合**：
- 全部 raw 計數 **matches=true ×3**（459/49/0.1068/410/281/0 + by_shape + dose-response 逐位元重現）。
- **裁決收斂：BOUNDED_AUXILIARY（不足以作獨立/主要標記）**，支撐改由「**覆蓋限制 + 0 偵測 + cis-control 可評估後多為 germline-cis（11/14）**」三線，非初版的「無訊號」。**（2026-06-29 更新：原「cis-control 不可評估」是 str/int bug，已修正為 14/49 可評估、11 germline-cis；裁決方向不變且證據更強。）**

🔴 **必修正既有報告（L4 / L5 / HTML）的「全 subclone-specific（0% cis）」**（2026-06-29 cis-control bug 修正後）：
| 處 | 原 | 改 |
|---|---|---|
| `05_methylation_corroboration.md` §2/§4 | 全 49 subclone-specific（0% cis）| cis-control **14/49 可評估，11 解析為 germline-cis（甲基非 subclone）**，僅 3 候選；subclone-specific **多被否決** |
| `06_integrated_narrative.md` §2/§3.5 L4 行 | 全 subclone-specific（0% cis）| 同上 + 補 power-gated（高功率 55%）|
| `index.standalone.html` L4 box | 全部 subclone-specific（cis 0%，proxy）| cis-control 14/49 可評估、11 germline-cis → 甲基多為 allele-specific 非 subclone |

**底線結論不變**：甲基 = auxiliary / ⭐3 上限（存活）。崩解的只是「**獨立 subclone 表觀支持 / 0% cis**」這一層。

## §7 給使用者的誠實回答

> **「甲基資訊/分群結構/差異對齊位點，能否提升 clone/subclone 準確度？作為標記/輔助是否足夠有效？」**

1. **不能當偵測器或主要標記**：0 個新 subclone 由甲基獨立發現；full_tree（subclone 最多）corroboration 最低。
2. **作為事後佐證（corroborate）有條件有效但有界**：覆蓋夠深（popB_n≥20）時 55% 區的甲基差異對齊基因型切割 → 高覆蓋區可提供正交信心；但這類區**僅 51 個**（ONT 覆蓋限制）。
3. **🔴 cis-control 修正後可評估（str/int bug，2026-06-29）→ 佐證多為 germline-cis 而非 subclone**：cis-control **14/49 可評估**（先前「0/740 不可評估」是型別 bug 假象），其中 **11 解析為 germline 等位特異甲基化（HP1≠HP2）、僅 3 候選**。即 corroborated 的 all-or-nothing Δβ **where testable 多被證為 germline cis-ASM 而非 subclone** → **更坐實「甲基非 subclone 標記」**。要下「甲基支持 subclone」之語，那 3 候選 + LOH 單倍型區仍須 **matched-normal 的 germline cis-control**（B 軸 deferred）。
4. **準確度提升**：對「基因型骨幹已確立」的 clone/subclone 結構，甲基**不增量**（不改任何 call）；少數高覆蓋區的弱佐證，cis-control 顯示多為 germline allele-specific。

**結論：甲基是「characterize 有界輔助」，不是「提升重建準確度的標記」。** 與論文 ⭐3 / auxiliary-not-driver 定位一致；2026-06-29 cis-control bug 修正後，更有正向證據（11/14 germline-cis）顯示甲基訊號多為 **allele-specific methylation 非 subclone**，裁決方向不變且證據更強。
