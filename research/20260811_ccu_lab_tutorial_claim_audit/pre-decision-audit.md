<!--
建立時間: 2026-08-11 23:55
目標: 在完整稽核 lab-tutorial 前，先確認問題可否證、證據是否足以支撐完整判定
處理範圍: 網站 25 個 live 頁面、InterSubMod 現行 authority、KB 與既有 NEGATIVE 結論
關聯檔案: InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/00_INDEX.md
-->

# Pre-Decision Audit：CCU lab-tutorial claim validation

## §0 Cynefin front-gate

- **Domain**：Complicated。
- **理由**：逐頁 claim extraction 與 evidence mapping 可重複；真正複雜處在不同 analysis unit／denominator 的詮釋，需專家判讀但不需先做新生物實驗。

## §1 Observation completeness

| 層級 | 現有證據 | 狀態 |
|---|---|---|
| L1 | 網站 source commit、live HTTP、authority artifacts、JSON/TSV counts、tool runtime KB | 完整 |
| L2 | 7 technical datasets／6 biological IDs、chr1–22；hash 與比例獨立重算 | 完整 |
| L3 | validated reports、audit cards、CURRENT_FOCUS 衝突修正 | 可用但不得凌駕 L1 |
| L4 | SR2/SR4/SR5/SR6 specs 與 toy derivations | 僅判定為 proposal |
| L5 | 首頁／教學類比／未附 provenance 的數字 | 需降權 |

## §2 Credibility score

| 維度 | 分數（/20） | 理由 |
|---|---:|---|
| 理論一致性 | 18 | molecule→cell、DNA fraction→purity 的界線可由明確模型檢查 |
| 直接觀察 | 20 | 網站 source 與現行 artifact 皆可讀取 |
| 機制可驗證性 | 18 | 重要主張可對照 tag type、CN/LOH 狀態與 denominator |
| 反例風險 | 10 | 已存在多項 NEGATIVE／confound 結論，依規自動降 10 |
| 資源可行性 | 20 | 唯讀掃描與重算，不需長計算 |
| **總分** | **86/100** | **GO** |

**Falsifier observable**：若所有疑似問題都能在同一 scope、同一分母、同一工具版本下被 authority artifact 支持，且沒有站內矛盾，則「網站需重大修正」假設被否證，結論降為 minor errata。

## §3 Assumption map

| 重要性／確定性 | 假設 | 處置 |
|---|---|---|
| 高／高 | read ≠ cell；f ≠ p；exact topology 尚無 cellular truth | 當 P0 ceiling |
| 高／中 | 網站工具 benchmark 與 repo run 可直接比較 | 僅在 dataset/version/scope 相同時比較 |
| 高／低 | 所有未附來源數字都錯 | 不採用；先重算再判 |
| 低／高 | print-all 是其他頁複本 | 檢查可達性，科學 claim 不重複計數 |

## §4 Quick pilot

已以首頁、M09、M12、M13、SR2b/SR2c 做快速 probe；發現至少三種可機械確認的問題：

1. 分子／細胞層級混用；
2. tumour DNA fraction／cellular purity 站內矛盾；
3. 巢狀分母與 reproducibility／accuracy 混用。

Probe 已超過 GO threshold，進入全站稽核。

## §5 Gap diagnosis

| 缺口 | 影響 | 優先級 | 本輪處置 |
|---|---|---:|---|
| 若干 tool 內部參數只見 KB，未釘正式 release | 中 | P1 | 標示 version-sensitive，不做無條件錯誤判定 |
| SR6 在 source 存在但 live 404 | 中 | P1 | 分離內容科學性與部署可達性 |
| 部分 paper benchmark 缺 CI／原始分母 | 中 | P1 | 限定 paper-scope，不外推 |
| 無 single-cell／multi-region truth | 高 | P0 | cellular clone／lineage 一律不得升格 |

## §6 Evidence conflict scan

| 既有結論 | Tier | 與網站關係 | Source |
|---|---:|---|---|
| 甲基不能獨立區分 cell groups／決定樹邊 | L1–L2 | 直接衝突 | InterSubMod/docs/CURRENT_FOCUS.md；20260801 authority |
| TO 多類特徵 AUC 低，安全 FP removal=0% | L2–L3 | 限縮 M7/M12/M13 的一般化敘事 | research landscape／ledger |
| exact topology 是 model-conditional graph shape | L1 | 限縮 SR1/SR2/SR2c | authority manifest |
| paired-pure ΔF1 約 +0.0112、具樣本異質性 | L1–L3 | 支持小幅 PoC，反駁普遍增益 | machine TSV／C11 |
| MEMORY.md | — | repo 未找到 | 以 ledger／validated NEGATIVE 取代 |

## §7 Decision path

- **Verdict：GO（86/100）**
- **允許範圍**：完成唯讀全站 claim audit、數字重算與本地報告。
- **禁止升格**：不因網站或工具作者的語句而突破現行 cellular／lineage claim ceiling。
- **紅隊檢核**：
  1. 教學簡化不自動等於錯誤；
  2. 外部 paper benchmark 不用本 repo 的 negative filter result 直接推翻；
  3. 數字正確但 denominator 不同時判為需補限定，不判造假。

