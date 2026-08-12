<!--
建立時間: 2026-08-11 23:55
目標: 完整稽核 CCU Bioinformatics Lab lab-tutorial 公開網站的科學敘述、數字口徑與現行 InterSubMod 證據上限
處理範圍: live 25 個 HTTP 200 頁面、1 個 source-only/HTTP 404 頁面、網站 Git commit 46b6f5b3016c187ad742fecbfa813f835b09e605、InterSubMod 20260801 authority bundle 與截至 2026-08-11 的研究狀態
關聯檔案: InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/pre-decision-audit.md; InterSubMod/docs/reports/validated/2026/08/20260811_CCU_lab_tutorial網站科學敘述完整稽核_01.md
branch: feat/lineage-tag-methylation-axes
commit: 73afaeac8e61c767241fa59c1ca6043a1c95290c
worktree: dirty; 本任務不覆寫既有未追蹤 crosscheck HTML
-->

# CCU lab-tutorial 科學敘述完整稽核

## 任務契約

- **Task Type**：B — Comprehensive validation。
- **服務目標**：G2、G3、G4、G5。
- **敘述框架**：SCQA 主線；逐項採 Claim → Evidence → Verdict。
- **主要問題**：網站敘述哪些可保留、哪些需限縮、哪些被現行數據反駁、哪些缺 provenance？
- **固定網站版本**：46b6f5b3016c187ad742fecbfa813f835b09e605（2026-08-11）。
- **完整範圍**：25 個 live HTTP 200 頁面；另將 sr6.html 列為 source-only、live HTTP 404 的部署差異。

## Pre-registration

1. **覆蓋率**：逐頁覆蓋 25/25 live pages；print-all 視為複本但仍檢查可達性與繼承錯誤。
2. **判定類別**：
   - 支持／可保留
   - 合理但需補限定
   - 需修正
   - 被現行數據反駁
   - 缺證據／不可判
3. **證據優先序**：
   - L1：機器 artifact、程式碼、hash、原始 count
   - L2：跨資料集或獨立重算
   - L3：validated report／已審查整合
   - L4：方法假說／spec／toy model
   - L5：教學直覺／未附來源敘述
4. **Hard claim ceiling**：
   - read 是分子觀測，不是 cell label；
   - local exact-PS×HP state/tree 不等於 cellular clone lineage；
   - methylation 目前是 association-only；
   - CN/LOH/purity/multiplicity/CCF 尚未整合進 exact topology；
   - tumour DNA fraction 不等於 cellular purity。
5. **反證條件**：若網站精確語句直接違反 authority manifest 的 forbidden claims、站內自相矛盾，或其數字無法在指定分母重算，即列需修正或被反駁。

## Step → Verify

1. 固定網站 commit 與 live route 清單  
   → **驗證**：Git commit=46b6f5…e605；25 routes HTTP 200、sr6.html HTTP 404。
2. 逐頁盤點 claim 與精確數字  
   → **驗證**：25/25 live pages 都有 page-level verdict；重點 claim 有 source line/link。
3. 對照 authority／KB／ledger  
   → **驗證**：authority artifacts hash 重算通過；每個 P0 claim 至少一個 L1/L2 source。
4. 重算關鍵比例與巢狀分母  
   → **驗證**：算術誤差 ≤ 0.01 percentage point；M1 screen 與 formal assay 不混分母。
5. 產出 Markdown＋standalone HTML  
   → **驗證**：HTML 無外部依賴、結構完整、number provenance／taboo audit／link audit 通過。

## 交付

- InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/pre-decision-audit.md
- InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/route_status.tsv
- InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/key_number_recheck.tsv
- InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/page_verdict_matrix.tsv
- InterSubMod/docs/reports/validated/2026/08/20260811_CCU_lab_tutorial網站科學敘述完整稽核_01.md
- InterSubMod/docs/reports/validated/2026/08/20260811_CCU_lab_tutorial網站科學敘述完整稽核_01.standalone.html

## 最終判定

- **CONDITIONAL PASS — MAJOR REVISION REQUIRED**。
- **P0**：index、M9、M12、M13、SR1、SR2b。
- **P1 工具事實**：M6/M7 的 HP tag vocabulary/type 與 platform support。
- **可保留核心**：M0、M1、M4、M8、M10、capstone；精確數字仍需 provenance/version/scope。
- **主要數據補充**：M9 screen 數字可重現，但 M1 與 formal assay 不可混用；88.26% 是 mathematical shape，不是 clone prevalence；paired +0.0112 是 read-classifier PoC，不是 variant-call F1。
