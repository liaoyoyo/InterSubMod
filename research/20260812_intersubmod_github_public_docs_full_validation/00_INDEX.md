<!--
建立時間: 2026-08-12
目標: 索引 InterSubMod GitHub README、Wiki、Pages 與 CCU lab-tutorial、演算法、CLI、數據的全面驗證證據鏈
處理範圍: GitHub main、遠端 feature/develop、獨立 Wiki、17 個 live Pages、CCU lab-tutorial 25 live routes+1 draft、C++/Python source、全樣本 authority artifacts
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/pre-decision-audit.md
  - InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md
-->

# InterSubMod GitHub 公開說明與教學完整驗證

## 任務契約

- Task type：B — Comprehensive validation。
- 服務目標：G2、G3、G4、G5。
- Scope：不得抽樣；涵蓋 GitHub 預設分支 README/QUICKSTART/PROJECT_SUMMARY、獨立 Wiki 全頁、GitHub Pages 全 17 頁、公開遠端 feature/develop 版本，以及被公開文字引用的演算法、CLI、runtime、數據與外部方法主張。
- 判定類別：`CONFIRMED`、`CONFIRMED_WITH_LIMITS`、`NEEDS_CORRECTION`、`CONTRADICTED`、`UNVERIFIABLE`。
- 科學證據層級：L1–L5；任何 cellular clone、lineage、causal、accuracy 或 prevalence claim 不得由低層級工程證據升格。

## Pre-registration（confirmatory）

| 預測 H | 否證條件 | Decision threshold |
|---|---|---|
| 公開說明存在版本漂移，且至少一個讀者首先接觸的 public plane 與目前 authority claim ceiling 衝突 | `main`、Wiki、Pages 內容與 commit 完全一致，且逐 claim 均未超出 authority | 任一 live/default plane 出現版本分裂或高影響矛盾即判 `NEEDS_CORRECTION` |
| 至少一個精確百分比存在 denominator／grain 錯置 | 所有比例重算均與文字宣稱的母體一致 | 任一比例的文字母體與計算分母不一致即判 `CONTRADICTED` |
| 公開 quickstart 不是外部 fresh clone 可直接重現 | fresh clone 具備所有示例輸入或自動下載步驟，且公開命令 exit 0 | 缺任一必要 input、依賴內部絕對路徑，或未揭露前置資料即判 `NEEDS_CORRECTION` |
| 工程測試 PASS 不足以確認 cellular subclone／lineage | 存在獨立 cellular truth、CN/LOH-aware 校正及跨生物樣本驗證 | 若 authority 明列 confirmed cellular subclone = 0，任何「已重建／已確認」用語均須降級 |

## Step → Verify

1. 鎖定四個公開版本面 → 驗證：記錄 commit、HTTP、SHA-256、頁面母體與部署漂移。
2. 拆解 atomic claims → 驗證：所有精確數字、絕對敘述、演算法與操作指令均有 inventory row。
3. Fresh build／test／CLI／source audit → 驗證：完整命令、exit code、stdout/stderr 片段與 source location。
4. 全樣本數值與 authority 重算 → 驗證：分母、效應量、CI、hash 與資料 grain 明列。
5. 外部原始來源交叉查證 → 驗證：論文或官方文件直接支持 claim，並區分推論。
6. 形成逐條確認／質疑／修正清單 → 驗證：每項均含 claim、證據、推論、反證條件與最小修正文案。
7. 產生 standalone HTML＋inline SVG → 驗證：來源連結、桌機／手機／no-JS／print 與禁忌文字 QA。

## 子任務契約

| 子任務 | 輸入 | 預期輸出 | 驗證標準 |
|---|---|---|---|
| 公開 claim inventory | README、Wiki、Pages、live routes | 全頁母體與 atomic claims 分類 | 無漏頁；精確數字與絕對語句全列 |
| 量化與 runtime 重驗 | authority、build、tests、公開數字 | fresh 重算表 | input／command／exit／actual output／denominator 齊全 |
| 演算法與 CLI 對照 | source、tests、教學文字 | claim→symbol→runtime evidence matrix | 每個結論有 file:line 或明示未找到 |

## 版本基準

- GitHub default `main`：`635437a65a33f8ba698acf85b22ebb069455c6cc`。
- 遠端 `develop`／`feat/lineage-tag-methylation-axes`：`ddd8909a838318d8a77969313e9561c8ff9d01c2`。
- Wiki git root：`6cfc990f1fbe0b0aa76b8adbe56bbe645ef7ee7b`。
- 本地稽核 HEAD：`73afaeac8e61c767241fa59c1ca6043a1c95290c`；公開文件與遠端 `ddd8909a` byte-identical。
- Live Pages：17/17 HTTP 200，逐檔 SHA-256 與遠端 feature 公開文件一致。
- CCU lab-tutorial baseline/current：`46b6f5b3016c187ad742fecbfa813f835b09e605` → `9eb1618d359e602d9c528675952b20d051fb2346`；25/25 正式 routes HTTP 200 且 byte-match，`sr6.html` draft 仍 404。

## 最終統計

- InterSubMod：158 個 deduplicated claim families；69 `CONFIRMED`、31 `CONFIRMED_WITH_LIMITS`、28 `NEEDS_CORRECTION`、26 `CONTRADICTED`、4 `UNVERIFIABLE`。
- 問題 family：58/158（36.71%）；直接 actionable correction：54/158（34.18%）；P0=34、P1=20、P2=4。
- Artifact worst-verdict：31 現存 artifacts 中 23 `CONTRADICTED`、6 `NEEDS_CORRECTION`、2 `CONFIRMED_WITH_LIMITS`；另 1 missing endpoint。
- CCU current-delta：32 findings；17 `OPEN`、6 `PARTIAL`、6 `REGRESSED`、3 `RESOLVED`，29/32 尚未完全解決。此為 problem-focused denominator，不是全站句子錯誤率。

## 交付物

- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/algorithm_cli_matrix.tsv`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/verification_results.tsv`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/chr2_18M_independent_audit_fresh.json`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/html_qa/qa_receipt.json`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/artifact_checksums.sha256`
- `InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md`
- `InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.standalone.html`
