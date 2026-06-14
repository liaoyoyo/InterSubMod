<!--
建立時間: 2026-06-12
狀態: handoff worklist（多 agent 稽核校正 → 待修正交辦單；本 session 不擅自改 Hard-Gate/他 session 檔）
報告類型: paper_focus_correction_worklist
受眾: 廖子游 · evidence_ledger owner · CURRENT_FOCUS / 整合篇章 / 正式學術版 owner session
data_sources: docs/paper_focus/00_共識證據台帳_20260612_01.md,docs/paper_focus/00_論文章節材料對應總表_20260612_01.md,research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv
provenance_note: 校正項出自 2026-06-12 paper-claim-provenance-audit（wf_ba28f1e3-a85）§3 矛盾清單；same-hap 一項經本 session 親驗 primary TSV。本單只列「需改但本 session 不擅自動」的檔；master 表自身 4 處已校正不在此。
-->
<!-- provenance-verified: 校正來源為共識台帳 §3 + 親驗 obs18 TSV；build_commit 069cadb。 -->

# 🔴 待修正交辦單 — Hard-Gate / 他 session SoT（2026-06-12 稽核校正）

> **L0**：稽核找到 **13 處跨文件矛盾（🔴 6）**。其中**本 session 自有檔（master 對應表）的 4 處已校正**（見該檔 §0）。**本單列「需改但屬 Hard-Gate 或他 session owner、本 session 不擅自動」的項**，交你/owner session 決定何時改。
> **L1**：① evidence_ledger 覆寫 = Hard-Gate（永遠需你確認）；② CURRENT_FOCUS/整合篇章/正式學術版 = 他 session owner（並行互撞風險，不擅改）；③ 其餘 paper_focus 檔本 session 可改，列「可即改（待你 go）」。每項附**驗證來源**。

---

## A. 🔴 Hard-Gate — `research/autoresearch/evidence_ledger.jsonl`（必經你確認）

| 位置 | 現值（stale）| 應改為 | 驗證來源 | 嚴重 |
|---|---|---|---|---|
| ledger:95 附近（same-hap overstated_flag）| same-hap「3/6 ≥93%（overstated，非 6/6）」 | **撤回 overstated_flag**；恢復 **6/6 ≥93% 並標 metric = same-hap occupancy（phasing 占用率，非 same_HP1 TP-rate）** | 親算 `research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv`：occupancy=0.932/0.990/0.988/0.965（4/6 直驗，與稽核吻合）| 🔴 |

> ⚠ evidence_ledger 為 append-only SoT，覆寫屬 Hard-Gate。建議**不改舊行，append 一條 correction entry** 註明 ledger:95 度量誤標 + 正確 6/6 occupancy，保留歷史軌跡。

---

## B. 🔴 Live 主軸 — `docs/CURRENT_FOCUS.md`（他 session owner / 你）

| 位置 | 現值 | 應改為 | 驗證來源 | 嚴重 |
|---|---|---|---|---|
| :137（同-hap 引用處）| same-hap「3/6 ≥93%」/ 引用 bug | 6/6 occupancy ≥93%（標度量）；修引用行號 | 同 A | 🔴 |

---

## C. 🔴 他 session SoT — `docs/concepts/2026/06/20260611_..._整合篇章_01.md`

| 位置 | 現值 | 應改為 | 驗證來源 | 嚴重 |
|---|---|---|---|---|
| §2.6（same-hap）| 「same-hap 中位 ~92% 但僅 3/6 ≥93%（0.840/0.939/0.759/0.429/0.932/0.920，非 6/6）」 | **occupancy 6/6 ≥93%**；那組 6 值是 same_HP1 TP-rate（準確率）被誤標成占用率 | 同 A | 🔴 |
| §2.5 / §5（chr17）| chr17/TBC1D16「唯一**乾淨** cis」 | 「**最強 cis 候選**（copy 軸已控）」；揭露殘餘 allele-axis d=0.183>d_within、未過 Cochran gate、n=30 Bonferroni-marginal | 稽核矛盾 #4 | 🔴 |
| §reconcile（subclone）| （已 flag G-B，但 Methods-Neg 不可寫成已 NEGATIVE）| subclone-甲基 somatic-specificity **undetermined**（押 G-B null）；勿寫 READY-AS-NEGATIVE | 稽核矛盾 #5（紅隊 REFUTED ×5）| 🔴 |

---

## D. 🔴 他 session SoT — `docs/paper_focus/02_paper_framework/論文架構_正式學術版_Slide2Thesis格式.md`

| 位置 | 現值 | 應改為 | 驗證來源 | 嚴重 |
|---|---|---|---|---|
| title / abstract_en | "only **1/332,705** loci is a clean somatic cis" | "1 of **816 testable** loci (~0.12%)" + 顯式分兩分母（catalog 332,705 vs cis-tested 816）| 稽核矛盾 #2 | 🔴 |
| §2.3（catalog floor）| 「乾淨 somatic cis 全基因組僅 **1/332,705**」 | 同上 1/816 | #2 | 🔴 |
| §2.5（chr17）| 「clean somatic cis exemplar」 | 「最強 cis 候選」+ caveats | #4 | 🔴 |
| §2.6（same-hap）| 「same-hap 中位 ~92%、3/6 ≥93%」 | 6/6 occupancy | A | 🔴 |
| §2.3 / §4.3（reliable 計數）| 12,868 vs 11,689 vs 12,876 並列當「reliable 總數」 | 加口徑註：12,868=TAG-C / 11,689=reliable(TP+FP-only) / 12,876=reliable(all-classes) | #3 | 🟡 |
| §2.3（none 計數）| 跨口徑補數 292,762 | **291,518（TAG-G）**；禁跨口徑補數 | #6 | 🟡 |

---

## E. 🟡 可即改（paper_focus，本 session 可改，待你 go）

| 檔案 | 校正 | 驗證來源 | 嚴重 |
|---|---|---|---|
| `01_focus_notes/02_文件庫...md` | same-hap「只 3/6 ≥93%」記述本身 stale → occupancy 6/6 | A | 🔴 |
| `02_paper_framework/位點甲基分群catalog_結果_R6.md` | floor 1/816；reliable 口徑註；none=291,518；chr18:11741161 改 genuine T3（非 cis-mechanical）| #2/#3/#6/#13 | 🟡 |
| `02_paper_framework/pre_registration_與fallback...md` | floor 1/816；reliable 口徑；chr17 軟化 | #2/#3/#4 | 🟡 |
| `01_focus_notes/03_研究方向卡...md` + `04_對應表` | chr17 軟化；schema 16→**32 欄**；刪「condition_fp_consensus.json 檔不存在」stale（檔實存）| #4/#7/#8 | 🟡 |

---

## 統計

- 矛盾總數：**13**（🔴 6 / 🟡 7）。本 session 已校正（master 表）：4。本單待修：A-E 共涵蓋全部 13。
- **0 個 F（捏造/孤兒）**：所有校正都是「口徑/範圍/誤標」修正，無數字憑空捏造。
- 完整逐條對抗證據 → `docs/paper_focus/_claim_audit_20260612/`（audited / challenge / contradictions JSON）。

## Provenance footer
- build_commit `069cadb`；稽核 run `wf_ba28f1e3-a85`（51 agents / 606 claim）。
- 共識台帳：`docs/paper_focus/00_共識證據台帳_20260612_01.md`（§3 矛盾清單原文）。
- same-hap 親驗：`research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv`。
