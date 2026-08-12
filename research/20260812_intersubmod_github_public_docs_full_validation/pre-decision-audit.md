<!--
建立時間: 2026-08-12
目標: 在全面稽核公開 GitHub 教學前，驗證投入價值、完整性、衝突與 fail-closed 邊界
處理範圍: InterSubMod GitHub main、Wiki、Pages、remote feature/develop、repo authority 與 negative evidence
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/00_INDEX.md
-->

# Pre-decision audit

## 1. 決策與 Cynefin

- 決策：是否投入 full-scope、逐 claim、可重現的公開文件驗證，而非只做抽樣文字校對。
- Cynefin：**Complicated**。版本面、source、數據與外部文獻可由專家分析與 deterministic tests 解開；不需要先承諾新的生物假說。
- 任務類型：B — Comprehensive validation；subset 不可當 final evidence。

## 2. 完整性與版本 pilot

已先完成最低成本 pilot：

1. `git ls-remote --symref origin HEAD` → default branch 為 `main`，commit `635437a...`。
2. clone 獨立 Wiki → commit `6cfc990...`，共 8 個 Markdown 檔（7 content pages＋sidebar）。
3. live Pages 17 個 HTML 全部 HTTP 200，逐檔 hash 與遠端 feature/develop `ddd8909...` 的 `docs/explain/` 一致。
4. GitHub `main` README 仍宣稱「利用 read-level 甲基化模式區分不同細胞群體」，與 `confirmed cellular subclone=0`、methylation association-only 的 authority ceiling 衝突。
5. Wiki／新 README 的「全部突變中有三分之二是孤立單點」先驗算即發現 denominator 疑義：`170,131/469,849=36.21%`，66.52% 的分母其實是 strict components `255,752`。

結論：版本母體已足以界定，且 pilot 已找到可 falsify 的 P0 衝突；全面稽核具必要性。

## 3. Assumption map

| 假設 | 重要性 | 不確定性 | 快速驗證／處理 |
|---|---:|---:|---|
| GitHub `main`、Wiki、Pages 可視為同一版本 | 高 | 已否定 | 分成三個 public planes 各自判定 |
| live Pages 等同遠端 feature 公開文件 | 高 | 低 | 17/17 SHA-256 MATCH |
| authority manifest 足以支撐全部生物 claim | 高 | 高 | 僅支撐凍結範圍；外部/因果/細胞 claim 另查原始證據 |
| fresh runtime 可代表外部使用者體驗 | 中 | 高 | 分開「本機有依賴」與「fresh clone 自足性」 |
| 精確數字自洽即可支持敘述 | 高 | 已否定 | 另驗 denominator、grain、effect、CI 與 generalization |

## 4. 五維評分

| 維度 | 分數（1–5） | 理由 |
|---|---:|---|
| 科學潛力 | 5 | 公開 claim 直接影響研究可信度與 G5 |
| 機制清晰度 | 4 | claim→source/data 可追，但多個版本面分裂 |
| 可重現性 | 5 | source、tests、authority artifacts 與 live hashes 可機械驗證 |
| 新穎性 | 3 | 稽核方法不新，但可首次建立完整 public truth map |
| 成本效益 | 4 | 主要是 read-only 重算與 build；不需重跑大型 BAM pipeline |
| **合計** | **21/25** | 超過 GO 門檻 |

## 5. Conflict scan

- 既有 hard negatives：TO 安全約束下 FP removal=0%；TO 跨樣本 transfer 不成立；methyl robust association 3/811 且 association-only；confirmed cellular subclone=0；CN/LOH 未整合。
- 既有 retractions：KDE stale-binary claims、Thread B cross-sample filter、self-phasing／V5 attribution 均有正式修訂。
- 對公開主張的衝突：以 methylation 或 read cluster 宣稱 cell-group/subclone resolution、將 read-AF topology 寫成真 lineage、將同分子共現寫成無條件 direct cellular observation，均需限縮或反駁。

## 6. Red-team failure modes

1. **版本混用**：以本地 feature README 修正內容替 GitHub default `main` 辯護，漏掉外部使用者實際首先看到的舊 README。
2. **工程 PASS 升格科學 PASS**：hash/test/exit 0 只能確認可執行與資料完整，不能確認細胞亞群、演化方向或 accuracy。
3. **absence 誤判**：3/811 或 0 confirmed 不能自動解讀為生物訊號不存在；需區分 true negative、underpowered、not tested 與 fail-closed abstention。
4. **fresh clone 假象**：在內部伺服器成功不代表 GitHub clone 有範例 BAM/FASTA、可依 README 重現。

可推翻本次 GO 的條件：若全部 public planes 同一 commit、所有高影響 claim 均有匹配 grain/denominator/原始證據，且 fresh clone 指令全自足通過，則可停止擴大；pilot 已否定此條件。

## 7. Verdict

**GO_WITH_FAIL_CLOSED_SCOPE**。

- 允許：完整 claim audit、fresh build/test、小型 fixture runtime、authority 重算、外部原始來源核對、HTML/SVG 報告。
- 不允許：把本輪文件稽核當成新 biological validation；不重跑全量 BAM pipeline；不以未測項作 negative；不直接改公開內容或 push。
- 最終輸出必須分清：GitHub main、Wiki、Pages/feature、repo current science 四個層次。

