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

🔴 **BOUNDED_AUXILIARY — 甲基不足以作為獨立或主要的 subclone 標記/偵測器。** 甲基差異**確實存在且強烈受覆蓋度閘控（power-gated）**：高功率區（popB_n≥20）corroboration 達 **54.9%**；但 (1) ONT 覆蓋度使**僅 51 區達此功率**，(2) 甲基**從未獨立偵測任何新 partition**（0；corroborate≠detect），(3) **cis-control 在 0/740 區可評估** → 這些 corroboration **無法與 germline cis-ASM 分離**。→ 對應既有 **auxiliary-not-driver / ⭐3** 上限。**［L1 源碼/數據重現］**

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

## §4 🔴 cis-control 不可評估 → 「subclone-specific」UNVERIFIED（最關鍵修正）

**hp_control_eval = 0 / 740**（含 **0 / 49 corroborated**）。

cis-ASM 控制（在 sig CpG 上比 HP1{1,1-1} vs HP2{2,2-1}）需「兩 germline 單倍型各 ≥3 reads」——此門檻在**全 740 區從未滿足過一次**。原因：corroborated 區以 **LOH 為主**（單一單倍型），且 genotype-population 本身被單倍型分離（somatic sSNV phased 到單一 HP）→ 結構上無法在 tumor 內做 HP1-vs-HP2 對照。**［L1 源碼 `sm_methyl_reextract.py:143` + 740 行 TSV 重算］**

🔴 **故 `cis_explained=0` 是「測試迴圈從未執行的 structural zero」，非「檢查後排除 cis」。** 推論：
- **「全 49 subclone-specific（0% germline cis）」UNVERIFIED** —— 精度為零（n_evaluable=0，true cis fraction CI=[0,1] 無資訊）。
- corroborated 的 **max_dbeta median=0.974、全 ≥0.857（all-or-nothing）= 正是 cis-ASM 特徵**（allele-linked methylation 跟隨基因型）。因 somatic sSNV phased 到單倍型，**任何 germline-haplotype-linked 甲基差異都會自動對齊 genotype 切割** → 在單樣本資料中**無法分離 subclone vs cis-ASM**。
- 與既有 memory 一致：`within_clean≠subclone`（Jaccard 0.123）、cis-ASM/double-dip confound。

## §5 邊際效用 — 0 新 partition + 反向趨勢

- **new_partitions_detected_by_methylation = 0**：甲基從未獨立偵測任何基因型遺漏的 subclone。
  - ⚠ 部分屬 **by-construction**（corroboration 設計只測「對齊」不搜尋「新切割」）；偵測端的 null 另由專案獨立的 **tumor-only 無監督分群 NEGATIVE** 佐證（見 memory）。
- **by tree shape**（corroboration rate）：full_tree **2.6%（1/39）** < sibling 4.8% < linear_nested 6.8% < co_linked 9.7%。→ 🔴 **結構最複雜、subclone 最多的 full_tree 反而最低**，與「好的 subclone 標記」應有的方向相反。
- corroborated 中 **18.4%（9/49）僅靠單一 CpG**（最弱證據）。**［L1］**

## §6 對抗稽核裁決 + 對既有報告的修正

**3 獨立 verifier（從 raw TSV 重算）+ evaluator 綜合**：
- 全部 raw 計數 **matches=true ×3**（459/49/0.1068/410/281/0 + by_shape + dose-response 逐位元重現）。
- **裁決收斂：BOUNDED_AUXILIARY（不足以作獨立/主要標記）**，但支撐改由「**覆蓋限制 + 0 偵測 + cis-control 不可評估**」三線，非初版的「無訊號」。

🔴 **必修正既有報告（L4 / L5 / HTML）的「全 subclone-specific（0% cis）」**：
| 處 | 原 | 改 |
|---|---|---|
| `05_methylation_corroboration.md` §2/§4 | 全 49 subclone-specific（0% cis）| cis-control **0/49 不可評估** → subclone-specific **UNVERIFIED**；corroboration 無法與 germline cis-ASM 分離 |
| `06_integrated_narrative.md` §2/§3.5 L4 行 | 全 subclone-specific（0% cis）| 同上 + 補 power-gated（高功率 55%）|
| `index.standalone.html` L4 box | 全部 subclone-specific（cis 0%，proxy）| cis-control 不可評估（0/740）→ 未驗證 |

**底線結論不變**：甲基 = auxiliary / ⭐3 上限（存活）。崩解的只是「**獨立 subclone 表觀支持 / 0% cis**」這一層。

## §7 給使用者的誠實回答

> **「甲基資訊/分群結構/差異對齊位點，能否提升 clone/subclone 準確度？作為標記/輔助是否足夠有效？」**

1. **不能當偵測器或主要標記**：0 個新 subclone 由甲基獨立發現；full_tree（subclone 最多）corroboration 最低。
2. **作為事後佐證（corroborate）有條件有效但有界**：覆蓋夠深（popB_n≥20）時 55% 區的甲基差異對齊基因型切割 → 高覆蓋區可提供正交信心；但這類區**僅 51 個**（ONT 覆蓋限制）。
3. **🔴 目前這份佐證的「獨立性」未成立**：cis-control 在本資料 0/740 不可評估，corroborated 的 all-or-nothing Δβ 是 cis-ASM 特徵 → **可能近乎全是 cis-ASM 而非 subclone**。要下「甲基支持 subclone」之語，**必須先補 matched-normal 的 germline cis-control**（B 軸 deferred）。
4. **準確度提升**：對「基因型骨幹已確立」的 clone/subclone 結構，甲基**不增量**（不改任何 call）；只在少數高覆蓋區提供未經 cis 校正的弱佐證。

**結論：甲基是「characterize 有界輔助」，不是「提升重建準確度的標記」。** 與論文 ⭐3 / auxiliary-not-driver 定位一致；唯一須修正的是把過度的「subclone-specific（0% cis）」降級為「cis-control 未評估、未驗證」。
