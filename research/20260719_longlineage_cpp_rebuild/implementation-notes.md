<!--
建立時間: 2026-07-19 05:38 +08:00
目標: LongLineage C++/HTSlib 全流程重建實作過程 living document
處理範圍: P0-P8
cycle_id: cycle_20260719-longlineage-cpp-rebuild
spec_id: longlineage_cpp_rebuild
status: in_progress
advisory: on
關聯檔案:
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/00_INDEX.md
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/pre-decision-audit.md
  - /big7_disk/liaoyoyo2001/LongLineage/AGENTS.md
-->

# Implementation Notes: LongLineage C++／HTSlib 全流程重建

> **Purpose**: prospective紀錄設計決定、偏離、折衷、未決與gotchas；≤400行。

## 🔵 設計決定

### [2026-07-19 05:38] 專案、執行邊界與主張上限

- **Status**: Accepted（user-specified）。
- **決定**: We will建立新 Git history的 `/big7_disk/liaoyoyo2001/LongLineage/`，namespace/CLI
  使用 `longlineage`，production science全部為 C++；Python只讀 `VALIDATED_FROZEN` chart-ready
  artifacts建立圖片與雙語HTML。
- **影響**: production、evaluation、presentation分成不可逆向依賴的三個trust domains。
- **Revisit if**: 使用者明確 supersede，或公開前license audit要求變更。
- **Evidence tier**: L1（直接使用者決策）。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: production不得接受truth BED/VCF/label；latest HP/PS只可來自7/7 raw-all sidecar exact join；
不得建立第二份tagged BAM；claim ceiling為long-read sSNV co-occurrence與lineage-compatible
mutation-state families。
**DO NOT change**: 不得回退至raw BAM舊HP、其他sidecar、推測tag，亦不得把candidate稱真實clone/祖先/時間順序。
**Rationale**: 2026-07-19使用者提供的鎖定計畫。
<!-- END USER-SPECIFIED -->

### [2026-07-19 05:40] 資料契約與查詢是一級產品介面

- **Status**: Accepted（user-specified addition）。
- **決定**: We will為每種JSON/JSONL/TSV/BGZF artifact建立versioned schema、field dictionary、
  primary key、ordering、nullability、units、coordinate convention、enum與conservation rules；
  同時提供validated-only read-only query CLI與cookbook。
- **後果**: `run_id → artifact → schema → record key → input identity → producer/validator receipt`
  必須可雙向追溯；schema drift與未知欄位預設fail closed。
- **Revisit if**: 新artifact需要schema major version；只能用migration ADR升版。
- **Evidence tier**: L1。

### [2026-07-19 05:42] Phase state不得用scaffold冒充完成

- **Status**: Accepted。
- **決定**: We will以 machine-readable phase ledger記錄 `NOT_STARTED / IN_PROGRESS / BLOCKED /
  VERIFIED`；production kernel只要不是 VERIFIED即回傳明確非零退出碼。
- **後果**: 初期repo可build與測試，但不會輸出偽 `VALIDATED` receipt。
- **Revisit if**: 無。
- **Evidence tier**: L2（由既有 v2 NO-GO incident支持）。

## 🟠 偏離之處

### [2026-07-19 05:44] 本機開發工具路徑與production pin

- **Status**: Closed（2026-07-19 07:00 environment audit）。
- **規範要求**: production OCI pin HTSlib 1.18。
- **實測結果**: 本機 `pkg-config --modversion htslib` 為exact `1.18`；
  `/usr/bin/cmake`為3.22.1。PATH先找到的另一個CMake為3.14.4，因此所有canonical
  command固定寫 `/usr/bin/cmake`。
- **風險**: OCI完整build與cross-compiler matrix仍是release gate，不能由local build取代。
- **Revisit if**: dependency lock、container base digest或compiler matrix改版。
- **Evidence tier**: L1（本session命令輸出）。

## 🟡 折衷考量

### [2026-07-19 05:46] Packed native artifacts優先於legacy小檔

- **Status**: Accepted。
- **方案 A**: We will以BGZF chunked/indexed records作native SoT，legacy exporter按需轉換。
- **方案 B**: 延續約1,409,547個per-site小檔；拒絕，inode/I/O與原子發布風險過高。
- **方案 C**: Parquet；deferred，會增加跨語言runtime與schema依賴，v1先用HTSlib-native BGZF。
- **Revisit if**: native query效能或schema evolution證明BGZF不敷需求。
- **Evidence tier**: L2。

### [2026-07-19 05:48] 新repo與InterSubMod provenance分離

- **Status**: Accepted。
- **方案 A**: We will新建Git歷史，以 `ORIGIN.md`＋source-to-target manifest綁來源hash。
- **方案 B**: 直接subtree搬移歷史；拒絕，會混入大量無關研究與真實資料路徑。
- **方案 C**: 只貼程式不記來源；拒絕，無法做license/scientific audit。
- **Revisit if**: 公開前法務要求保留更完整commit history。
- **Evidence tier**: L1。

### [2026-07-19 07:10] Schema bundle、release attestation與run receipt不可依賴mutable state

- **Status**: Accepted。
- **決定**: Production manifest除scientific contracts外，必須綁定
  `schema/id_registry.json`與immutable `state/release_attestation.json` SHA-256。
  Registry逐schema保存physical SHA；catalog逐artifact保存record-schema SHA。
- **後果**: 同一SemVer schema原地變更會使preflight失敗。`longlineage run`不得直接讀日後會改變的
  `state/phase_ledger.json`判定P3/P4/P5；NOT_READY attestation固定回
  `KERNEL_BLOCKED`，不會產生production結果。
- **Revisit if**: schema major migration；必須新ADR與negative compatibility test。
- **Evidence tier**: L2（獨立contract red-team finding）。

### [2026-07-19 07:20] Record/query/receipt membership採closed-world DAG

- **Status**: Accepted。
- **決定**: Native record契約除欄位型別外，明定artifact/index key、semantic digest、
  per-unit lineage、receipt chain與run-file membership。Query只讀`VALIDATED_FROZEN`，
  filter是AND-only closed AST；incomplete scan的`matched_rows`必為null。
- **後果**:
  - topology unit保存FULL/PARTIAL input pattern、multiplicity與兩個per-unit evidence digest；
  - query response綁run/producer/validation/checksum與artifact/index digests；
  - tied topology不得輸出winner；
  - artifact catalog、lineage、semantic digest、producer/checksum/validation/final receipt
    的self-exclusion與creation order由machine registry驗證，避免hash cycle。
- **Revisit if**: P6 validator發現membership無法獨立重建；此時只能升schema major。
- **Evidence tier**: L2（雙人獨立契約稽核＋synthetic positive/negative schema validation）。

## 🔴 未決問題

### [2026-07-19 05:50] GPL來源公開前最終授權相容性

- **Status**: open。
- **Question**: InterSubMod來源片段、HTSlib、HiGHS與新repo的GPL-3.0-only組合是否通過逐檔audit？
- **Default if no answer**: 保持private，禁止public visibility與v1.0.0。
- **Revisit if**: 準備公開或引入第三方source。
- **Priority**: major。
- **Evidence tier**: L5。

### [2026-07-19 05:52] Direct HiGHS pinned source取得與build

- **Status**: open。
- **Question**: 是否能以可重現方式取得/固定 `scipy/highs@4a122958`並通過C++ API測試？
- **Default if no answer**: exact B&B/DP可繼續開發；需要HiGHS的case標 `ABSTAIN_DEPENDENCY_UNAVAILABLE`，
  不得release P5/P7。
- **Revisit if**: dependency audit完成。
- **Priority**: critical。
- **Evidence tier**: L5。

### [2026-07-19 07:25] HP family與Endpoint-B precedence缺frozen vectors

- **Status**: open／phase blocker。
- **Question**: normalized HP `0/1/2/...`到HP family的版本化映射，以及Endpoint-B O/X
  callability precedence，尚無可綁定的完整golden vectors。
- **Default if no answer**: P3/P4保持`BLOCKED`；不得臨時建立分組、合併O/X或從名稱推測family。
- **Revisit if**: authority與logical SHA完成freeze。
- **Priority**: critical。
- **Evidence tier**: L1（source/contract inventory缺口）。

### [2026-07-19 08:20] AI task與audit evidence採machine-enforced closed registry

- **Status**: Accepted／implemented at foundation。
- **決定**: 每個agent task記錄owner、parent/dependency DAG、typed write-set、
  bounded lease、heartbeat、scope與Step→Verify；active claim不得互相包含或重疊。
  Audit不接受Markdown宣稱，而以immutable envelope綁task、scope、Git commit、
  canonical tree SHA、完整argv與stdout/stderr SHA。
- **後果**: unknown field、`.`/`./` alias、glob、symlink ancestor、missing
  dependency、DAG cycle、expired/future/stale lease、偽audit與supersession cycle
  均由C++ governance負向測試拒絕。
- **限制**: envelope中的command digest不會自行執行；replay必須由明示授權命令
  重新執行並比對，不能把digest存在當作測過。
- **Evidence tier**: L1（Debug/Release/ASan-UBSan CTest 25/25）。

### [2026-07-19 08:25] 資料紀錄與查詢採六層可逆索引

- **Status**: Accepted／contract-implemented。
- **決定**:
  1. schema ID registry綁path＋physical SHA；
  2. artifact catalog綁record schema、primary/index key與canonical order；
  3. closed TSV registries管理type/status/lifecycle/transform/query operator；
  4. source-to-target manifest將presence、implementation、verification拆開；
  5. immutable receipt/audit chain綁run/task與semantic digest；
  6. production query只接受`VALIDATED_FROZEN`並使用closed typed AST。
- **後果**: canonical JSON欄位另有nested schema；Interval0 pair有machine-readable
  semantic group；未知欄位、missing/extra/duplicate binding、schema SHA drift與
  lifecycle倒退全部fail closed。
- **Evidence tier**: L1（22 offline schema IDs、16 artifacts、87 status/reason rows、
  positive/negative fixtures PASS）。

### [2026-07-19 08:40] Foundation測試通過不升級scientific phase

- **Status**: Accepted／verified for stated scope。
- **結果**: Debug、Release、ASan/UBSan各25/25 synthetic tests PASS；
  `check_all.sh`在Debug與Release均`FOUNDATION_PASS`，no-network test PASS。
- **阻擋**: strict gate仍以12個fixture-only／missing executable negative
  evidence項目退出1；P0/P1/P2/P6保持`IN_PROGRESS`，P3/P4/P5/P7/P8保持
  `BLOCKED`。
- **Evidence**:
  `InterSubMod/research/20260719_longlineage_cpp_rebuild/20260719_LongLineage_foundation實作與驗證紀錄_01.md`。
- **Evidence tier**: L1。

### [2026-07-19 09:20] 正式本機repo與雙層audit closeout

- **Status**: Completed for local foundation；remote仍blocked。
- **結果**: `/big7_disk/liaoyoyo2001/LongLineage/`已建立clean `main`，HEAD
  `3a789d3c8b384606dfad01ae0227834df01661ff`，199 tracked files且無build
  產物。`002` audit以相同task/scope supersede `001`，兩個Git tree digest均可
  從正式repo重播。
- **Fresh target verification**: 以正式路徑作source、`/tmp/LongLineage-target-smoke`
  作Release build；configure/build exit 0、CTest 25/25、foundation/no-network
  gate PASS。
- **仍未完成**: GitHub private remote；credential／create-repository authority
  不存在，保留`PRIVATE_REMOTE_NOT_YET_VERIFIED`。
- **Evidence**:
  `InterSubMod/research/20260719_longlineage_cpp_rebuild/20260719_LongLineage_foundation實作與驗證紀錄_01.md`。
- **Evidence tier**: L1。

## 📚 Lore

### [2026-07-19] Sidecar identity並非只有QNAME

- **Constraint**: canonical header為
  `#CHROM START0 END0 QNAME FLAG MAPQ CIGAR_B2 HP PS`；paired QNAME需依FLAG加 `/1`/`/2`，
  strand來自FLAG，CIGAR_B2為8-byte BLAKE2b digest。
- **Why it matters**: projection collision不可用QNAME-only或position-only解決。
- **Evidence**: `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/latest_tag_join.py`。

### [2026-07-19] MM/ML方向與unknown skip語意

- **Constraint**: MM以as-sequenced 5'方向編碼；`?`與`.`/省略的skip語意不同；ML為
  `[N/256,(N+1)/256)` bin；MN存在時必須等於SEQ length。
- **Why it matters**: reverse alignment或skip語意誤讀會改M1矩陣與所有下游決策。
- **Evidence**: KB `03_file_formats/bam-format.md`，last_verified 2026-07-12。

### [2026-07-19] 現有InterSubMod worktree高度dirty

- **Constraint**: 2026-07-19 `git status --short`顯示大量既有modified/deleted/untracked檔。
- **Why it matters**: 本cycle只能新增獨立topic與新repo，不可清理、還原或混入既有變更。
- **Evidence**: 本session實際命令輸出。

### [2026-07-19] Bootstrap subagent不可倒填事前lease

- **Constraint**: 三個subagent在machine lease/write-set機制完成前已開始寫入。
- **處理**: 明文記錄process deviation，不事後捏造task/lease；機制生效後由root
  單一lease完成收斂，未來delegation必須先登記child task與不重疊claim。
- **Why it matters**: audit trail要區分「當時存在的授權」與「後來建立的規則」。

### [2026-07-19] Sandbox中的LeakSanitizer失敗不是程式memory finding

- **Constraint**: 未設定`detect_leaks=0`時，LSan回報其不能在目前tracing boundary下
  運作；該次run被中止且不列正向證據。
- **處理**: 保留AddressSanitizer與UndefinedBehaviorSanitizer，僅關閉LSan後重跑
  25/25 PASS；leak-only檢查留給支援LSan的CI/container。
- **Why it matters**: 不可把環境不支援誤報成程式PASS或memory bug。

## Provenance

- InterSubMod baseline commit: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`。
- Pre-decision: `InterSubMod/research/20260719_longlineage_cpp_rebuild/pre-decision-audit.md`。
- Status: `in_progress`；未通過P0-P8前不得finalize。
