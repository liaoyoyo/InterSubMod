<!--
建立時間: 2026-06-26
類型: methodology audit — 全基因組 sSNV 單分子連鎖 pipeline 對抗稽核 + 主 agent 獨立複核
狀態: NEEDS_WORK（fresh-context evaluator verdict + 主 agent 數據確認）
data_sources: sm_linkage_genomewide.json, sm_completeness_ledger.json, sm_verify_null.json, sm_seqc2_overlay.json
-->

# 全基因組 sSNV 單分子連鎖 pipeline — 對抗稽核（NEEDS_WORK）

> Fresh-context `evaluator` agent 獨立稽核 5 腳本 + 5 JSON，verdict **NEEDS_WORK**。主 agent 對其最尖的 F1/F2/F3
> **獨立複核（讀回 landed JSON）確認屬實**。本檔保存裁決 + 確認數據 + 修正清單。所有數字落檔可重算。

## Verdict: NEEDS_WORK — 5 findings

| # | Finding | Severity | 主 agent 複核 |
|---|---------|----------|--------------|
| F1 | `tally_powered_somatic` 是 **per-pair 非 per-sSNV** → dense 群以 O(n²) 灌水 | major | ✅ 確認：distinct linked sSNV=**21,554** vs pairs=**53,094**；max **90** links/sSNV，226 個 sSNV ≥80 links，top 1% 佔 21% link-degree。10,843 nested 等是 pair-events 非 distinct 變異數 |
| F2 | somatic gate（normal-VAF<0.05）admits **mapping/segdup 偽影**（非 somatic 證明） | critical | ✅ 確認：FP-somatic-linked 含經典偽影簇 |
| F3 | 「3,204 γ 類 = SEQC2 漏掉的 clonally-informative somatic」**OVERCLAIM**；3204≈3200「收斂」近**套套邏輯**（同 census 子集，只差 SEQC2-region 限制） | critical | ✅ 確認：見下偽影簇；convergence 非獨立 |
| F4 | (a) Fisher exact ≠「exact allele-label permutation null」（Fisher 固定兩 margin）；承諾的 MC shuffle 未交付。(b) **negative control diffHP 50% 顯著** → Fisher-sig **零 subclone-鑑別力**，鑑別全靠 HP-gate | major | ✅ 確認：diffHP sig/ns=4511/4465；headline「92% sig」誤導 |
| F5 | 「no-miss」僅 **Tier-R**（≤50kb same-read）；same-PS（Tier-PS）deferred → isolated_singleton 是上界非真值 | minor | ✅ 確認：已標註但須一致 |

**Evaluator 確認正確處**：chr17 nested 方向映射正確（α→β2）；ledger sum-check 真實（8320+5458+21554=35332）；component scope 在程式碼層 red-line-compliant（per-chrom ≤50kb 局部區域非 genome-wide tree）；ledger per-sSNV 分桶無 double-count。

## 主 agent 獨立複核數據（讀回 landed JSON）

### F1 — per-pair 灌水
- distinct **LINKED sSNV = 21,554**（乾淨 per-sSNV 分母）；**n_pairs = 53,094**。
- n_realized_links per sSNV：max **90**、median(linked)=**2**、mean=**4.1** → 多數 linked sSNV partner 少，少數 dense 群驅動 pair 暴增。
- ≥20 links: 808 個；≥50: 426；≥80: 226。top 1% linked-sSNV 佔 **21%** 總 link-degree。
- **結論**：headline 用 pair-count（10,843 nested / 10,641 co_linked / 3,949 sibling）誤導；須改報 **distinct sSNV 參與數** + 揭露 dense-cluster 灌水。

### F2/F3 — 偽影簇（artifact signature）
FP-somatic-linked = 3,204；其中 **54 個 ≤5kb-簇 ≥5 sSNV = 816 sSNV（25%+）** 帶經典偽影簽名（uniform-low-VAF + tight-clustering + 已知 segdup 區）：

| 簇 | 範圍 | sSNV 數 | VAF 範圍 | 訊號 |
|---|---|---|---|---|
| chr8:81.39–81.47M | 76kb | **97** | 0.10–0.26 | segdup-like |
| chr8:82.93–82.98M | 50kb | 82 | 0.11–0.28 | segdup-like |
| chr8:82.12–82.18M | 60kb | 65 | 0.09–0.27 | segdup-like |
| chr9:41.78–41.79M | 12kb | 40 | 0.10–0.17 | uniform-VAF |
| chr9:41.80–41.80M | 4.5kb | 38 | 0.11–0.16 | uniform-VAF |
| chr14:16.09–16.10M | ~3kb | 6+ | **0.027–0.042**（極均勻~3%） | mapping-artifact，coread~400 |

→ 這些**不是**獨立 clonally-informative somatic，是 cross-mapped reads / segdup 偽影（同一批 mis-mapped read 同時帶多個「ALT」→ 也過 linkage+Fisher）。normal-VAF<0.05 擋不掉（偽影多映到 ref 故 normal-ALT 低）。

## 修正清單（在任何 HTML/doc/論文宣稱前必做）

1. **F1**：報 **distinct sSNV 參與數**（per-sSNV）為 headline，pair-count 為輔並揭露 dense-cluster 灌水；考慮 collapse dense cluster 後再計。
2. **F2**：somatic gate 加 **mappability/segdup mask + CN/LOH context**（SAVANA per-segment CN，見 memory `project_ont_cnv_sv_subclone_verification_feasibility`）；語言降為「normal-VAF-filtered **candidate**」非「somatic-confirmed」。
3. **F3**：γ 類重framing 為「FP-labeled ONT 呼叫、somatic-like normal、單分子連鎖 —— **候選，需正交驗證（CN/mappability/multi-platform）**」；明列 segdup-artifact 為**首要替代解釋**；移除「兩獨立推導收斂」（套套邏輯）。
4. **F4**：改 method string（Fisher = conditional 2×2 test 非 free permutation）；交付 MC allele-shuffle 對照或移除該句；**明寫 Fisher-sig 必要非充分、不分 subclone/allelic（HP-gate 才分）**；承認 diffHP negative-control 50% 顯著。
5. **F5**：ledger 全標「**Tier-R completeness（same-PS deferred）**」。
6. **可重算性 gap**：產小 summary JSON（非 15M-token 單行 dump）讓獨立 evaluator 能 recompute tally/component aggregate。

## 不變的 solid 結論（adversarial 後仍站得住）
- 方法可行且**高效**（group-based single-pass，~40min 全基因組）。
- chr17:48360161 worked example **被 pipeline 自動完整重建**（α→β2 nested + γ∥α sibling，VAF 0.82>0.48>0.18，0 違反）。
- **完整性帳本（Tier-R）真實**：35,332 sSNV 每個歸一桶，sum-check ✓；21,554 distinct linked。
- null **sparse 2.2%≈機率底線** vs 結構類高顯著 → power-gate 不從噪聲製造結構（此為 null 唯一真正有效部分）。
- ⭐3 單樣本 = 分子證據非 single-cell confirmation；6,146 components = 局部 haplotype-連鎖區域非 subclone。
