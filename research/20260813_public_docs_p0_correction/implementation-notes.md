<!--
建立時間: 2026-08-13
目標: 公開文件 claim remediation 的 living record
處理範圍: 158 claims；本輪 bounded write-set 為公開入口、Wiki/Pages source、generator 與 claim registry 工程
status: IN_PROGRESS_WITH_FAIL_CLOSED_GATES
-->

# Implementation Notes：公開文件 claim remediation

## 🔵 設計決定

### [2026-08-13] frozen inventory 是 ID 母體

- **Status**：Accepted。
- **決定**：158 個 claim ID 必須與 2026-08-12 frozen inventory 完全相同；修正不覆寫原稽核。
- **理由**：保留 before/after 與 locked public snapshot。
- **驗證**：registry validator檢查 duplicate、missing、extra與總數 158。

### [2026-08-13] source、live、release 三層狀態分離

- **Status**：Accepted。
- **決定**：source_status 只表示 checked-in source；live_status 只接受重新抓回的 public bytes；release 由獨立 Gate 判斷。
- **理由**：本機改字不能證明 GitHub main／Wiki／Pages／About 已更新。
- **驗證**：除有獨立API re-fetch receipt的C108外，所有claim的live state保持UNVERIFIED_AFTER_20260812_LOCKED_SNAPSHOT；模擬local source-live collapse必須失敗。

### [2026-08-13] P0 guard可 PASS，但整體 Source Gate維持 BLOCKED

- **Status**：Accepted。
- **決定**：33個local P0標SOURCE_READY；C108以獨立live receipt標EXTERNAL_LIVE_CONFIRMED_WITH_LIMITS；初始24個P1/P2中，C047/C048/C089/C114已有可追溯bounded evidence closure，其餘20個標SOURCE_EDITED_REVIEW_REQUIRED。
- **理由**：避免用 34 個 P0 的完成度吞掉完整 58 個問題 claim 母體。
- **驗證**：P0 source guard PASS；完整registry為69 confirmed、69 confirmed-with-limits、20 unverified。

### [2026-08-13] `p0_claim_registry.json` 是跨文件 guard 的唯一權威輸入

- **Status**：Accepted。
- **決定**：current 199-column schema、ISM／exact-PS solver 分界、denominator grain、dynamic test count、Python ≥3.10、tracked DEMO fixture 與 immutable/non-develop publication path 都由 `scripts/p0_claim_registry.json` 宣告；builder 不另藏第二份 policy。
- **理由**：只在 generated registry 或 validator 裡硬編規則會造成三份真相，且重新 build 可能默默遺失 gate。
- **驗證**：builder 先跑 P0 semantic validator，並將 guard registry path、SHA-256 與 PASS summary寫入 generated registry／receipt；完整 validator模擬 hash drift 必須 fail closed。

### [2026-08-13] Fixture「存在」升級為 Git blob 身分驗證

- **Status**：Accepted。
- **決定**：8 個 tiny fixture／runner path 必須同時是 non-symlink regular file、唯一 stage-0 Git index entry，且 working bytes 的 Git blob OID 等於 index OID；3 個 entrypoint另須 mode 100755。
- **理由**：本機有檔案不代表 fresh clone 收得到；只有路徑存在檢查會把 untracked fixture 誤判為已公開。
- **驗證**：回歸測試建立「磁碟存在但未 git add」fixture，validator必須拒絕。

### [2026-08-13] GitHub About C108獨立live receipt

- **Status**：Accepted／executed。
- **決定**：About description限縮為local mutation-state candidate analysis與methylation association；移除`subclone`／`phylogeny` topics。API action後立即re-fetch，保存repository identity、before/after、topics與時間。
- **理由**：remote surface不能由local source status代替，但單一remote PASS也不能吞掉default branch、Wiki、Pages或release gates。
- **驗證**：`github_about_c108_receipt.json` hash由builder與validator固定；hash drift regression必須fail closed。

### [2026-08-13] Browser QA變成repository-owned runner

- **Status**：Accepted。
- **決定**：17份 public explain HTML 加 4 份 handoff HTML，共 21 頁，固定跑 desktop 1440、mobile 390、no-JS 390與print 1440；禁止external runtime request，檢查document overflow、page/console error、HTML/XML/SVG parse，print只在記憶體產生。
- **觀察**：第一輪真實抓到page07 desktop `1506/1440`與mobile `1196/390` overflow；修正generator的container/grid/table/long-token CSS。將所有現存handoff HTML納入後，最終需以 84-case clean-source receipt 作驗收。
- **觀察**：擴充為完整21頁後的首輪為80/84 PASS；current 20260813 handoff與historical 20260806 handoff各在mobile/no-JS有table overflow。其餘80 cases、desktop/print、page/console error與external request均通過；此輪刻意寫成`FAIL_CLOSED`，修正並完整重跑前不得升格。
- **驗證**：runner與receipt皆入Git；receipt明示不驗science、accessibility、非Chromium browser或publication state。

### [2026-08-13] 第二輪 adversarial semantic audit 擴充跨文件 guard

- **Status**：Accepted／executed。
- **決定**：新增 cohort grain、PERMANOVA、cis heuristic、LOH mechanism、methylation-topology bridge與LongLineage blocker-set六類 corpus guard；不只比對單一頁面的固定句子。
- **理由**：第一輪P0 anchors可防已知原句回歸，卻不足以攔截同義改寫、否定詞遺失或把7 technical datasets寫成7 biological samples。
- **驗證**：42個guards在34個public source files執行500次document checks、285 required anchors、1,217 forbidden anchors，0 errors；回歸測試數由本次runner receipt動態記錄，不手抄固定數量。

### [2026-08-13] 生成圖必須綁定最新HTML source bytes

- **Status**：Accepted／executed。
- **決定**：8組SVG/PNG由`tools/extract_svg_for_github.py`重建，不接受手動編輯衍生圖；每份SVG記錄source document SHA-256與inline SVG SHA-256。
- **理由**：HTML文字改正後若沿用舊PNG/SVG，README仍可能公開已撤回的語意。
- **驗證**：8/8 generate-and-verify及後續8/8 verify-only均PASS；PNG尺寸與SVG viewBox比例一致，HTTP(S)載入與active SVG content受拒絕。

## 🟠 偏離之處

### [2026-08-13] 只讀移植 dirty source 的受控 28-file patch

- **Status**：Accepted after byte comparison。
- **觀察**：原工作樹含多項不相干 dirty work，不能整包搬入 release。
- **決定**：只逐檔移植README EN/ZH、Quickstart、Summary、8 Wiki sources、17 Pages HTML及指定generator；不搬stale receipts、cache、screenshots或其他dirty內容。
- **驗證**：28 個檔案與已稽核 dirty target bytes一致；寫入使用 apply_patch。

## 🟡 折衷

### [2026-08-13] P1/P2 已修 source，但保守不升格

- **Status**：Accepted。
- **方案**：已依mapping修正storage、identifiability、statistics、test count、code inventory、normal wording、performance與novelty等source文字；缺少閉合evidence或完整source guard者仍一律 UNVERIFIED/SOURCE_EDITED_REVIEW_REQUIRED，registry列出 exact targets與 minimum rewrite。
- **理由**：歷史 1.67 TiB／287×/293.34×／99%、hard-coded test count 與 code inventory 都曾有未閉合說法。現已確認 14 個 PRE-FIX BAM 的 paths/bytes/mtimes，並以 2026-07-11 sampled-identity receipt 重播確認其中 7 個 `paired_full`（1,840,983,466,353 bytes；7/7 MATCH）；但未計算 full-file SHA-256、也未建立與現行 sidecar 的 producer revision／內容等價，所以仍保留 `PARTIAL/NON_FINAL`。exact stat bytes quotient 294.2669×只可稱跨世代 footprint comparison；舊 287× 無效，不假稱 production 或因果壓縮率已驗證。
- **Revisit if**：每個 claim新增 source anchors、evidence receipt與 live publication receipt。

## 🔴 未決

- 其餘 20 個 P1/P2 source-edited claims需後續逐claim guard、evidence closure與review。
- main／Wiki／Pages發布與 re-fetch未執行。
- branch protection、full-history secret scan與完整 release gate不屬本檔可單獨閉合。

## 📚 Lore

- 170,131 / 255,752 = 66.5219% 的 grain 是 strict components；469,849 是 dataset-records，兩者不同 grain，**不可**用 170,131 / 469,849 推導 isolated dataset-record percentage。凍結 audit 中的 36.21% 是無效跨 grain 推導，只留在 `inventory_evidence` 作 historical/invalidated provenance，不適用於 current verdict。
- 63,506 / 71,955 = 88.2579% 是 frozen recurrence-allowed model下的 graph-shape統計。
- LongLineage tagged-BAM能力必須帶revision與privacy/status：private baseline/main snapshot `5daf50f`無writer；private public-preview candidate `b9aaa12`有`longlineage-tag-bam`，但仍`NOT_READY`／non-production且P3/P4/P5/P7/P8 BLOCKED。
- JSON regex中的 literal parenthesis要寫成合法JSON escape \\(，不能用非法 \(。
