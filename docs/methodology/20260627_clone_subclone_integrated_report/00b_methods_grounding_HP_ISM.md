<!--
建立時間: 2026-06-27
類型: 方法學基礎 — longphase-S HP + ISM 甲基 定義/工作流/衝突/可信度/錯誤模式（clone/subclone 整合報告）
狀態: 基礎參考（源碼/KB grounded；L1-L5 皆引用此）
authority: src/core/ReadParser.cpp, src/core/DistanceMatrix.cpp, docs/method_comparison/.../01_ism_method_spec_from_source.md, docs/methodology/20260625_subclone_verification_nact_cistest_P1_3_01.md, docs/methodology/20260616_readset_provenance_observation_01.md
-->

# 方法學基礎 — HP（longphase-S）與甲基（ISM）的定義、衝突、可信度、錯誤模式

> 本檔在整合敘述前**先釘死兩個輔助資訊源的定義與可信度先後**，所有層引用。源碼/KB grounded（行號見 authority）。

## §1 longphase-S HP tag — 定義（source: `ReadParser.cpp:123-162`）

| code | string | 意義 |
|---|---|---|
| 1 / 2 | "1"/"2" | germline HP1/HP2（該染色體**無** somatic）— **驅動 phasing** |
| 11 / 21 | "1-1"/"2-1" | somatic 突變在 HP1/HP2 germline 背景上（**被動 annotate**，germline 主導）|
| 33 | "3" | somatic 方向不明 / 低信心 / 未 phase（**排除於方向投票**）|
| 0 | — | unphased read |

🔑 **HP1-1 ≠「HP1 的拷貝 1」**；= 「這條 read 帶著一個 phased 到 HP1 germline 單倍型的 somatic」。**HP1 ∪ HP1-1 = 同一 germline 單倍型**。

**工作流**：germline het VCF + somatic VCF + BAM → 2-pass phasing（**PON-only germline 主導、somatic 排除於 graph build → 避免 self-phasing 循環**）→ getVote 逐 read（germline-first 投票 → 無 germline 時 somatic fallback Layer1.5 → 仍不明=HP3；信心閾 0.6）→ 寫 HP tag。當前 production = **v5**。

## §2 ISM 甲基 pipeline — 定義（source: `01_ism_method_spec_from_source.md`, `DistanceMatrix.cpp:254-301`）

6 階段：① read×read 甲基距離（**BERNOULLI** 預設：δ=p(1−q)+(1−p)q 期望不一致，weight=2|p−0.5| 信心，Dist=Σw·δ/Σw；C_min CpG 重疊）→ ② UPGMA 分群（silhouette k=2..10）→ ③ PERMANOVA（pseudo-F, 999 perm）→ ④ per-CpG Fisher+BH-FDR → ⑤ Cramér's V + Cochran 閘 → ⑥ **NACT**（normal-anchored cis-test，扣 germline cis）。輸出 `methylation.csv`(read×CpG, 窗 ±5kb 中位 76 CpG) / `distance/BERNOULLI/matrix.csv` / `phylo_groups.tsv`。

## §3 🔴 核心可信度先後（precedence）— 三層

| 層 | 資訊源 | 循環性 | 可信度 | 角色 |
|---|---|---|---|---|
| **1（最高）** | **sSNV read-level 連鎖** | **非循環**（同分子物理共現）| 直接觀測 | 克隆共現的 ground truth |
| **2** | **HP tag（genetic phasing）** | 非循環（genetic）但**有 phasing 誤差** | per-read ~85–90%（**clean PS 上的跨版本 concordance，非本資料集直接量測**）| sibling vs allelic 鑑別 |
| **3（最低）** | **甲基（epigenetic）** | **循環**（若用甲基分群定 subclone）| 可在同克隆內變動 | corroborate genotype-anchored 群（**非 detect**）|

**衝突解決鐵則**：**genetic（sSNV 連鎖 > HP）> 甲基**。甲基與 genetic 衝突 → **genetic 勝**（甲基異質 = 同克隆內可塑性，非 subclone）。
**實證衝突**（`nact_cistest_P1_3 §`）：chr5:88881646 全 read hp=2-1+ALT（基因型 100% 同質）卻被甲基分 2:1 → **double-dip 假象**；HP-gating 擋掉。

## §4 為什麼甲基單獨不可信（cis-ASM/double-dip 循環，source: NACT P1.3）
1,139 甲基候選 → NACT 102「survivor」→ genotype-axis 重測**只 9/30 存活** → 這 9 = **somatic-allele-specific 甲基狀態（characterization）非獨立 subclone**。**結論：單樣本甲基無法非循環證 subclone**；genotype-同質群內甲基子結構**永遠 double-dip**。→ 甲基**只能 corroborate 已由 sSNV/HP 定好的群**。✅ 本報告 L4 用的正是 **genotype-anchored**（非 clustering），故非循環。

## §5 錯誤模式 + 本分析如何處理

| 錯誤模式 | 來源 | 影響本分析 | 處理 |
|---|---|---|---|
| **problem PS block**（**34 塊**，germline concordance 51-69%；「96%」= V5-vs-前版**改變的 sites 比例**，非「誤差集中率」— 修對抗稽核）| longphase TO 低-AF 難 phase | 落在這些塊的 region，HP 判 sibling/allelic **不可信** | 🔴 L3 PS 抽取後**標記**這些 region；乾淨結論排除之 |
| HP3 ungateable（v5 ~8%）| 無 germline 錨 | 1,317 linked-somatic 無法 HP-gate | L1 已標：這群不判 sibling/allelic |
| tumor-only axis-swap | 無 paired normal phasing | HP1/HP2 軸可能整體翻轉 | 不影響 same/diff-HP（相對），影響絕對 HP1/2 標籤 |
| self-phasing 循環 | somatic 進 phasing graph | 已 fix（PON-only v2b+）| 用 v5 tag，已修 |
| 甲基 double-dip | 甲基分群定 subclone | 會把同克隆甲基異質誤判 subclone | L4 用 genotype-anchored + NACT cis-control |
| sSNV mapping 偽影 | CN-gain segdup | 灌水 dense cluster（骨幹對抗稽核 F2）| 乾淨集 = LOH/neutral |

## §6 對整合的意涵（L5 將遵循）
1. **subclone 判定鏈**：sSNV 連鎖（共現）→ HP（同 germline 單倍型=克隆 vs 異=allelic）→ 甲基（corroborate）。**順序不可顛倒**。
2. **HP 可信度有上界**（~85-90% + problem PS blocks）→ 最乾淨的 sibling 結論應**排除 problem-PS-block region**（L3 標記後）。
3. **甲基永遠是佐證**，與 genetic 衝突時 defer；單獨不可定 subclone。
4. 三者**同源於同一批 ONT read**（同 BAM），故是「同分子上的 genetic + epigenetic 雙觀測」，可信度先後如 §3。
