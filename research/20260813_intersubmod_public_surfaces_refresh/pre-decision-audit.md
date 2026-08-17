<!--
建立時間: 2026-08-13 17:36 +08:00
目標: 在修正 InterSubMod GitHub 與 Pages、產出 CCU 唯讀改進清單前，先驗證範圍、證據權威與否證條件
處理範圍: InterSubMod 公開 README/QUICKSTART/PROJECT_SUMMARY、Wiki、docs/explain Pages；CCU 僅讀取既有 audit，不修改任何 CCU 檔案或 patch
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/research/20260813_public_docs_p0_correction/00_INDEX.md
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
-->

# Pre-decision audit — InterSubMod 公開介面資訊更新

## 0. 判定

- **Verdict：GO（80/100）**。
- Task type：**B — Comprehensive validation**；不可用 subset 冒充全部公開介面。
- 服務目標：G2、G3、G4、G5。
- Cynefin：complicated；公開文字可由 frozen authority、claim inventory 與實際 source 對照，但 HTML generator／衍生頁同步需要專業查核。
- 高影響理由：錯誤公開 claim 會把 read-level molecular evidence 升格成 cellular clone／lineage，或讓外部讀者誤判可重現性。

## 1. 啟動研究任務 5 問

1. Thread D 相關：**是**；read-level epigenetic 與 LOH/CN 邊界會出現在公開方法說明。
2. Thread B 撤回範圍：**不重啟**；只把「永不再研究」改成「已測 formulation 為 negative，若以實質新假說重啟須先 audit」。
3. KDE-corrected：**不適用於本輪文件修正**；不得用未鎖定的 KDE 結果替換 frozen authority。
4. VCF caller AF：**不重算**；若提及 AF，明確區分 caller-VAF、read-AF 與 noncanonical TVAF transform。
5. 長計算／C++／搬移／NO-GO gate：**不涉及核心 C++ 或全 BAM 重跑**；涉及多檔案公開說明與 NO-GO 邊界，故需 claim-level gate。

## 2. 決策假說與否證條件

| 假說 | 支持線索 | 否證條件 | 處置 |
|---|---|---|---|
| H1：P0 修正後剩餘 24 個問題 claim 可在本機公開來源逐一 disposition | 2026-08-12 inventory 將問題分成 P0=34、P1=20、P2=4 | 任一 P1/P2 問題沒有 occurrence、owner、最小修正文案或明確 external action | 24/24 才能通過 |
| H2：InterSubMod GitHub／Pages 可統一到同一 scientific claim ceiling | exact-PS authority 已列允許與禁止 claim | 任一來源仍把 molecule co-occurrence、甲基 association 或 noncanonical CCF 稱為 cellular truth | residual guard 必須為 0 |
| H3：CCU 可與 InterSubMod 修改完全隔離 | 使用者明示「只整理重點改進清單，不改動 CCU」 | 任一 CCU source、live site、既有 patch 或 remote state 被修改 | 立即 FAIL；只允許在本 repo 新增唯讀清單 |
| H4：Pages 的 26 個缺 title/desc SVG 可無語意漂移地補齊 | 2026-08-13 P0 QA 已定位 02/03/04/05/06/08/09 | SVG 結構、圖形幾何、可見文字或 HTML 互動被意外改壞 | 17 頁 QA、SVG geometry 與 browser gate 全過 |

## 3. 關鍵假設與不確定性

- 高信心：`authority_manifest.json` 是本輪 cellular／lineage／methylation／CN-LOH claim ceiling。
- 高信心：24 個待修問題的母體定義是 verdict=`NEEDS_CORRECTION|CONTRADICTED|UNVERIFIABLE` 且 priority=`P1|P2`，不是所有 P1/P2 row。
- 中信心：`docs/explain` 是 Pages editorial upstream；Wiki 是人工同步 derivative。必須在文件中明說，不能宣稱自動 single-source。
- ⚠ 待確認：live GitHub default branch、Wiki、About 與 Pages 仍需另外發布才能閉環；本輪只驗證本地 source，不把 remote 狀態標成已修。
- ⚠ 待確認：tagged-BAM 7 檔 exact byte receipt 尚不存在；在 receipt 出現前只能保留 sidecar exact bytes，tagged-BAM total 與 ratio 標為未驗證。

## 4. Conflict scan／Red team

- Thread B：2026-04-26 撤回的是跨樣本 whitelist filter 用途；不能寫成所有 F1／TO／filter 研究永久禁止。
- 甲基：只允許 genetic-pattern-conditioned regional association；不得移動 topology edge、證明 clone 或宣稱因果功能。
- read/HP/PS/MM/ML：是 called/recorded molecular evidence，受 sequencing、alignment 與 tag error 影響；cell identity 未被觀測。
- read-AF／TVAF：目前未整合 CN/LOH；只能作 model-conditional ordering 或 historical approximate transform。
- HCC1395 chr2:18M：biological n=1、單 sample/locus；technical cross-basecaller replication 不是 biological replication。
- 工程 PASS：測試、schema、required-metric guard 只能支持工程契約，不保證科學正確性。
- 版本漂移：直接改 standalone HTML 可能與 generator 不一致；有 generator 的頁面必須同步 source 並重建／比對。

## 5. 最小安全路徑

1. 凍結 24 claim inventory 與目前 dirty-tree 邊界。
2. GitHub 公開文件與 Pages 分工修正，CCU 僅由既有 audit 產生清單。
3. 建立 P1/P2 claim registry 與 residual guard；逐 claim 記錄 local-fixed／external-action。
4. 跑 17 頁 HTML/SVG/browser QA、連結檢查與 fresh-reader。
5. 產出 validated MD + standalone HTML；寫 evidence ledger，但不 push／merge／部署／改 About。

