<!--
建立時間: 2026-06-27
類型: L4 ISM 甲基 corroboration — 單變量（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3；🔴 既有甲基輸出覆蓋不足，genome-wide 待重抽）
data_sources: data/sm_methyl_corroboration.json, data/sm_methyl_reextract_merged.json
-->

# L4 — 甲基 corroboration（單變量；既有輸出不足 → 已從 BAM 重抽補完 genome-wide）

> 問題：① 既有 ISM 甲基輸出能否直接用於全基因組結構區的甲基 corroboration？② 甲基對遺傳定義群的輔助程度？
> 🔴 紅線：甲基 = **corroborate 非 detect**；只在 **genotype-anchored**（用 sSNV 定好的群，非 clustering）才非循環。
> 圖：`figures/05_methylation_corroboration.png`。資料：`data/sm_methyl_corroboration.json`。

## §1 決定性發現 — 既有甲基輸出覆蓋不足（不能直接全基因組用）

- 既有 ISM 甲基窗 = **1,139 個 cis-candidate anchor 窗**（records_v6 候選位置，**±5kb 寬**）。
- 這組位置與我的 TP∪FP 連鎖 sSNV **是不同集合** → 窗不對齊我的區域。
- 實測：CN-clean 結構區（span≤8kb）**323 個中，只有 9 個**有 covering candidate 窗（其餘 314 缺）。

→ **既有甲基輸出無法直接做全基因組區域甲基 corroboration**。🔴 **誠實分母（修對抗稽核）**：9/323 = 2.8% 是「CN-clean + span≤8kb」這個最有利子集的覆蓋率；對**全部 4,678 個有結構區域**，覆蓋率僅 **9/4,678 = 0.19%**。**全基因組需從 BAM（MM/ML tag）重抽甲基、對齊連鎖區域** = pipeline 整合步驟（後補；對應你「缺的後續整合成流程」）。

## §2 genome-wide 重抽結果（補完缺口；754 區，data_source: sm_methyl_reextract_merged.json）

**方法**：直接從 tumor BAM 的 **MM/ML（5mC）重抽**甲基、對齊每個連鎖區域（**非** ISM anchor 窗）→ 驗證 **corr=1.000 vs ISM methylation.csv**。單次 fetch 每 read 同時取甲基+基因型+HP。genotype-anchored：兩大基因型 population 間 per-CpG MWU+BH-FDR（sig: q<0.05 且 |Δβ|≥0.2）。cis-ASM control：sig CpG 在 HP1 vs HP2 是否也區分（=germline cis）。

| 指標 | 值 |
|---|---|
| 目標區（CN-clean 結構區 ≥2 population）| 754 |
| 可測（≥2 population 有甲基）| 740 |
| **甲基 corroborate（≥1 sig CpG）** | **49（6.6%）** |
| **cis-ASM 控制（HP1 vs HP2）** | 🔴 **0 / 49 可評估**（hp_control_eval=0/740）→ cis vs subclone **UNVERIFIED** |
| median sig CpG / 測試區 | **0** |

→ 🔑 **genome-wide-clean：6.6%（49/740）區域甲基能區分遺傳定義的 group**（median sig CpG=0）。**⚠ 詳細有效性審查見 `08_methylation_sufficiency_audit.md`**，兩項對抗稽核修正：
- **power dose-response**：6.6% 偏低主因是**覆蓋不足非無訊號** —— 高功率區（popB_n≥20）corroboration 達 **54.9%**，但僅 51 區達此功率（58.8% 的「powered」落在功率飢餓帶）。
- 🔴 **「全 subclone-specific（0% cis）」已撤回**：cis-ASM 控制（HP1 vs HP2 各≥3 reads）在 **0/740 區可評估**（corroborated 區多為 LOH 單倍型）→ `cis_explained=0` 是 structural zero（從未檢測）非排除；corroborated 的 all-or-nothing Δβ（median 0.974）= **cis-ASM 特徵** → 無法與 germline cis-ASM 分離。
→ **甲基 = 有界弱 corroborator**：少數高覆蓋區提供（未經 cis 校正的）佐證，多數區不區分 —— 符合「characterize 非 detect」。

## §3 chr17 canonical（已驗證，established）
chr17:48360161（有窗）先前完整分析：甲基 BERNOULLI/UPGMA **只切得出「有無 somatic」（L0/γ vs L1+L2），切不出 L1↔L2 fine subclone**；per-CpG 可歸因（23 CpG→α 祖先 ASM / 6→L1-L2）。→ 甲基 = **corroborate（刻畫已驗群）非 detect**，且細分 subclone 需遺傳錨。

## §4 data-supported 結論（誠實）
1. **既有 ISM 甲基輸出 = 稀疏 anchor 窗，覆蓋 0.19%（9/4678）結構區** → 故**已直接從 BAM 重抽**（驗證 corr=1.000），不靠既有輸出。
2. **genome-wide-clean 結果（740 區）：甲基只在 6.6%（49）區域 corroborate 遺傳 subclone**（power-gated：高功率區 54.9%，但僅 51 區達此功率）；🔴 **cis-control 0/740 不可評估 → 不能稱「subclone-specific」**（可能為 cis-ASM）；多數（93%）區域甲基不區分。→ 甲基 = **有界弱 corroborator**，不下「強佐證 / 獨立 subclone 表觀證據」定論（見 `08`）。
3. 🔴 甲基 corroborate 非 detect；細分 subclone 必須靠遺傳（sSNV 連鎖 + HP）。

## §5 後補流程（✅ 已執行 — `sm_methyl_reextract.py`）
- 從 tumor BAM 解析 MM/ML（5mC）tag → 對**每個連鎖區域**重抽 read×CpG 甲基矩陣（對齊 region，非 ISM anchor 窗）= **已完成**（§2）。剩餘：normal 端 NACT cis-control 用 HP-axis proxy（已含），完整 normal-baseline residual = 進一步精煉。
- 全基因組 genotype-anchored per_cpg_diff + NACT cis-control（normal）→ 真 genome-wide 甲基輔助程度。
- 這是論文「methylation profiles = characterize」一章的完整數據來源。
