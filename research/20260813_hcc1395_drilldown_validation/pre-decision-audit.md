<!--
建立時間: 2026-08-13 01:26
目標: 在修改 drilldown generator 或推廣到多樣本前，先驗證目前 HCC1395_v1 的證據、風險與可否證條件
處理範圍: HCC1395_v1 全產物稽核、generator fail-closed 與 UI 誠實呈現；不產生跨樣本科學結論
cycle_id: cycle_20260813-0126-hcc1395-drilldown-validation
topic: 20260813_hcc1395_drilldown_validation
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/AGENTS.md
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/scripts/build_drilldown_dashboard.py
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/implementation-notes.md
-->

# Pre-Decision Audit：HCC1395 drilldown 全面驗證與多樣本 fail-closed

> **Verdict：GO（70/100）只適用於唯讀稽核、低成本 generator hardening 與 fixture 驗收；跨樣本科學結論維持 NO-CLAIM。**

## Frontmatter

- **Topic**: `20260813_hcc1395_drilldown_validation`
- **Triggered by**: new-spec + cross-sample gate
- **AI session**: 2026-08-13
- **Last updated**: 2026-08-13 01:26
- **Cycle ref**: `InterSubMod/state/cycles/cycle_20260813-0126-hcc1395-drilldown-validation/`

## §0 Cynefin Domain Gate

- **Domain**: Complex
- **Test**: 同一份 generator 在 HCC1395 可生成，不代表換 sample root 後會可預測地拒絕錯樣本；資料 provenance、科學口徑與瀏覽器互動會交互影響。
- **Probe-first**: 先用 `COLO829 --probe-only` 與小型 fixture 嘗試否證 sample isolation，再決定修補範圍。

## §1 Observation Completeness Checklist

| Observation | 狀態 | Tier | 來源 |
|---|---:|---:|---|
| HCC1395_v1 共有 11,590 regions、29,572 methylation rows、16,302 loci | ✓ | L1 | `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/receipt.json` + shard recount |
| 32,604 PNG 與 5,410 IGV payload 全數可解析 | ✓ | L1 | 本輪全量 asset scan；0 missing/orphan/parse error |
| Selfcheck 實際為 10 PASS / 1 FAIL / 1 SKIP | ✓ | L1 | `InterSubMod` 外部產物 `HCC1395_v1/SELFCHECK.md` |
| v1 capability 的 MLHP 解釋已過時且與 C12 矛盾 | ✗ | L1 | v1 `receipt.json` 與 current `InterSubMod/scripts/drilldown/sources/mlhp.py` 對照 |
| 換成 COLO829 時 extension root 會自動阻擋 HCC1395 資料 | ✗ | L1 | `--sample COLO829 --probe-only` 仍回報 HCC1395 LCA/ISM 數值 |
| 有 ≥3 個獨立樣本 drilldown bundle 可作 cohort 統計 | □ | — | `drilldown_out/` 目前只有 HCC1395_v1、HCC1395_v3 |
| SEQC2 truth/HC BED 已整合於本 bundle | ✗ | L1 | receipt 無 truth VCF、HC BED、som.py/hap.py 結果 |

## §2 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | fail-closed sample identity、provenance 與 UI state 是標準資料產品契約 |
| 觀察支撐 | 10 | HCC1395 v1/v3 與一個 COLO829 probe 有直接反例，但尚無完整多樣本 bundle |
| 機制清晰度 | 20 | hard-coded roots + linkage 永遠 PASS + receipt 不記 generator provenance 可直接解釋失敗 |
| 反例風險 | 10 | 可用 synthetic fixture 驗證，但大型真實 bundle 尚可能暴露額外 schema 差異 |
| 所需資源 | 10 | fixture/code/UI 修補 1–6 小時；完整七資料集生成遠超本輪 |
| **TOTAL** | **70 / 100** | **GO（限定工程 hardening；非跨樣本 science GO）** |

**Falsifier observable**：若修補不充分，`--sample COLO829` 仍會把 HCC1395 topology/LCA/ISM 標成可用，或新產物仍能在存在 FAIL/SKIP 時顯示「整體沒有問題」。

**Reality-test 三反例觀察**：

1. 錯樣本 fixture 必須非零退出或把該 capability 標為不可用，不能只顯示低 linkage。
2. 1 FAIL / 1 SKIP 必須產生 blocking/incomplete status，不能與 `6/6 完整` 並列。
3. 沒有 `filterMap` 的矩陣不能顯示 pointer 或宣稱點擊可篩選。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| CLI `sample` 與每一核心/擴充來源屬同一樣本 | HIGH | UNKNOWN（現況反例） | (2) MUST validate first |
| v1 可當作 immutable 失敗證據，不應原地重建 | HIGH | KNOWN | (1) verify quickly |
| bundle 指標足以作 F1/precision/recall 科學驗證 | HIGH | KNOWN FALSE | (1)；claim ceiling=observation-only |
| v3 能完全取代 v1 | LOW | UNKNOWN | (4) defer；v3 亦為 dirty/provenance 不完整 |

## §4 Quick Pilot（已執行）

1. 讀取 `InterSubMod/scripts/build_drilldown_dashboard.py` 與 `InterSubMod/scripts/drilldown/capability.py`。
2. 以既有 site config 執行 `--sample COLO829 --probe-only`。
3. 門檻：任何 HCC1395-only extension 必須 unavailable；實際仍回報 LCA 4.969×、ISM 1/20,937，因此 probe **FAIL**，確認需要 fail-closed 修補。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---:|---:|---:|
| sample identity contract 與 mismatch regression test | HIGH | 1–2 h | P0 |
| receipt generator commit/source hash/claim ceiling/output inventory | HIGH | 1–2 h | P0 |
| Selfcheck、假互動、mobile sticky/overflow 誠實呈現 | HIGH | 1–3 h | P1 |
| methylation BH-FDR、effect/n 定義與 test-family contract | HIGH | 2–6 h | P1 |
| 七資料集同 schema drilldown bundles + cohort manifest | HIGH | >6 h | P2（另立 cycle） |
| SEQC2 truth + HC BED benchmark 整合 | HIGH | >6 h | P2（另立 validation cycle） |

## §6 Evidence Conflict Scan

- Repository root 無 `MEMORY.md`，此項明確記為 unavailable；未以外部記憶補寫。
- `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md` 與本輪工程 hardening 無直接衝突。
- `InterSubMod/docs/CURRENT_FOCUS.md` 與 `InterSubMod/docs/reports/in_progress/2026/08/20260809_整體流程與軟體輸入輸出格式盤點_01.md` 均支持「lineage 只是一個軸、v3 dirty build/provenance 有缺口」，故限制 claim ceiling。
- KB HCC1395/benchmark workflow 指定 truth VCF、HC BED、som.py 與 TP/FP/FN/F1 口徑；v1 receipt 缺這些輸入，故禁止從 dashboard 推導 caller F1。

## §7 Decision Path

- **Verdict**: GO for audit + fail-closed/UI/reproducibility hardening；NO-CLAIM for cross-sample science。
- **Red-team failure modes**:
  1. 只修文案但不阻擋錯樣本，會把可信度問題藏得更深。
  2. 用 HCC1395 v1/v3 兩個版本冒充兩樣本，會造成 pseudo-replication。
  3. 把資產完整性誤寫成科學有效性，會越過 truth benchmark gate。
- **Red-team result**: 通過；上述三項都有可觀察 regression gate。
- **Decision lock**: Y。原始 v1 不覆寫；本輪只做可逆 source/test/report 修改。完整多樣本生成需新 cycle、sample-specific manifest 與用戶確認計算範圍。

## Provenance Footer

- **Commit hash at audit start**: `73afaeac8e61c767241fa59c1ca6043a1c95290c`
- **Skill version**: `/pre-decision-audit` v0.1
- **Audit JSON**: `InterSubMod/state/cycles/cycle_20260813-0126-hcc1395-drilldown-validation/audit.json`
- **Next**: `InterSubMod/research/20260813_hcc1395_drilldown_validation/implementation-notes.md`

