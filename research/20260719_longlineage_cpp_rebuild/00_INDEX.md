<!--
建立時間: 2026-07-19 05:38 +08:00
目標: 以獨立 LongLineage C++/HTSlib 專案重建 read-linked methylation、sSNV 共現與 lineage-compatible mutation-state family 流程
處理範圍: P0-P8；7 datasets / 6 biological samples / chr1-22；production truth-isolated
關聯檔案:
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/pre-decision-audit.md
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/implementation-notes.md
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/20260719_LongLineage_foundation實作與驗證紀錄_01.md
  - InterSubMod/state/cycles/cycle_20260719-longlineage-cpp-rebuild/audit.json
  - /big7_disk/liaoyoyo2001/LongLineage/AGENTS.md
-->

# LongLineage C++／HTSlib 全流程重建

> **Task type**: B — Comprehensive validation；subset 只可由 `probe` 產生並強制標 `PARTIAL`，不可作 release evidence。
>
> **服務目標**: G2、G3、G4、G5。核心價值是把 latest HP/PS、read-level methylation、sSNV
> co-occurrence、topology 與 provenance 做成可獨立驗證的 C++ authority。

## 研究問題

能否在不執行 Python 科學程式、production 不接觸 truth、也不複製第二份 tagged BAM 的前提下，
以 C++17＋HTSlib 重建現有 frozen 方法，並以獨立 validator 證明資料、決策與輸出可重現？

## 假設

1. 既有 Python source、frozen outputs、synthetic goldens 與 receipts 足以定義 parity oracle。
2. raw BAM 與 7/7 latest HP/PS sidecar 可用 exact identity fail-closed join。
3. M1/M2/topology 可被拆成 typed I/O、deterministic kernels、packed artifacts 與獨立驗證四層。
4. C++ port 只有在 decision/status/digest parity 通過後才可取代歷史 Python authority。

## 成功條件

- P0-P6 的 phase gates、fault injection 與 deterministic semantic digest 全數通過。
- P7 以 24/40 workers 跑 7 datasets，兩次 semantic SHA 相同，input before/after SHA 不變。
- production manifest、CLI、receipt 與 report builder 均無 truth 欄位或 truth-aware code path。
- P8 只讀 `VALIDATED_FROZEN` C++ artifacts；所有數字可回指 receipt。

## 失敗條件

- 任一 parity decision/status/read-set/family digest 漂移。
- 任一 missing/extra/duplicate、HP/PS multimatch、partial BGZF 或 forged PASS 未被 validator 阻擋。
- cap/deadline case仍產生 winner。
- Python presentation重算任何科學統計或 production 讀到 truth。

## 固定評估指標

- schema validation、record conservation、semantic SHA、input/output SHA、fault-injection catch rate。
- M1 decision parity、M2 gate parity、topology objective/family/parent/tree-count parity。
- wall/user/system、RSS/cgroup peak、I/O、faults、threads、queue/reorder wait、latency quantiles。
- 全量 census 使用使用者鎖定的 469,849 site keys 與各階段硬門檻；未跑前一律標 `PENDING`。

## 主要資料與限制

- production source：raw BAM 的 sequence/CIGAR/BQ/MM/ML ＋ 7/7 raw-all sidecar 的 latest HP/PS。
- topology Python authority SHA-256：
  `1def0de1952d127d8d013820ac7b0eabe33d6f1a66fd8c6d47a1985b14a32f77`。
- 既有 formal full co-occurrence result不存在；v2 M2 real pilot曾 timeout並正式 NO-GO。
- 無 single-cell／colony／spatial truth；claim ceiling 固定為 lineage-compatible mutation-state families。

## 預計輸出

1. `/big7_disk/liaoyoyo2001/LongLineage/` 私有新 Git repo。
2. C++ executables、schemas、synthetic fixtures、independent validator、query/export/evaluate tools。
3. AI governance、ADR、CURRENT_FOCUS、living notes、CI/release/secret/large-file gates。
4. P0-P8 receipts、驗證報告與 sanitized bilingual standalone HTML。

## 目前可驗證進度（2026-07-19）

- 新repo foundation已在隔離staging建置：AI治理、資料schema/lifecycle/provenance、
  typed I/O、deterministic runtime、small-q oracle與fail-closed CLI。
- Debug、Release、ASan/UBSan各25/25 synthetic tests PASS；Debug/Release完整
  foundation gate與no-network test PASS。
- Strict release evidence仍以12個缺少executable negative evidence的gate阻擋；
  P0-P2/P6不因scaffold測試通過而標為`VERIFIED`。
- 正式本機repo已落在`/big7_disk/liaoyoyo2001/LongLineage/`，clean `main`；
  以該路徑fresh Release configure/build/CTest與foundation gate均PASS。GitHub
  private remote仍因credential/authority缺口保留blocked。
- 詳細命令、輸入／輸出、資料查詢路徑與未完成邊界見
  `InterSubMod/research/20260719_longlineage_cpp_rebuild/20260719_LongLineage_foundation實作與驗證紀錄_01.md`。
