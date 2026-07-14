<!--
建立時間: 2026-07-15 03:20
目標: 在 canonical layered workstation 補回 GRCh38 全基因分布與樣本全貌統計前完成 evidence audit
處理範圍: layered-workstation-sample-overview
cycle_id: cycle_20260715-0320-layered-workstation-sample-overview
topic: layered-workstation-sample-overview
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_topology_workstation.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
-->

# Pre-Decision Audit：Layered workstation 樣本全貌層

> **Verdict**: GO — 只重用舊頁的閱讀形式；數據、分母與 ontology 全部綁定 canonical v5。
> **Task type**: B｜Comprehensive validation；7 datasets、chr1–22，不允許 subset 代替終驗。
> **Service goals**: G3 / G4 / G5。

## §0 Cynefin Domain Gate

- **Domain**: Complicated。
- **Test**: 相同行動（由 compact region index 生成離線圖、再以 Chromium 驗收）已有可預測結果；難點在資料契約對帳，不是未知機制探索。
- **Decision**: 可使用 generator-first 與全量 contract checks；不需先啟動新的生物假說 cycle。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---:|---:|---|
| 舊頁有按 GRCh38 真實長度縮放的 22 chromosome SVG | ✓ | L1 | `build_topology_workstation.py:279-353` + current-run Chromium screenshot |
| 新頁 `#genome-overview` 只有 chromosome count grid，沒有 SVG | ✓ | L1 | `build_layered_workstation_v5.py:901-919`；Chromium `genomeSvgCount=0` |
| 舊 HCC1395 7,143 / 6,288 / 35,332 屬 retired backbone | ✓ | L1 | `InterSubMod/docs/CURRENT_FOCUS.md` 2026-07-09、2026-07-14 |
| current HCC1395 為 W_tree 8,222、W_primary 7,932、complete 7,151 | ✓ | L1 | `current_layered_topology_v3_raw_all_v1.json` canonical sample record |
| compact index 已含 chrom/start/end、n_sSNV、C/Topo、evidence、HP、CN | ✓ | L1 | `build_layered_workstation_v5.py:167-191` |
| 舊頁 mobile 實際水平溢位 | ✓ | L1 | 390px viewport：document width 506px |
| 舊 single/linear/branched/star 能否安全映射 current region | □ | L4 | 待 current candidate/primary-unit grain audit |

## §2 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | GRCh38 座標與 current compact index 都是明確資料契約 |
| 觀察支撐 | 20 | 舊/新實際 render 與 7/7 canonical summary 已取得 |
| 機制清晰度 | 20 | 缺圖原因為 renderer 未實作，不是 browser failure |
| 反例風險 | 10 | 舊 morphology ontology 可能無安全一對一映射 |
| 所需資源 | 10 | 全 7 頁生成與多 viewport 稽核約 1–6h |
| **TOTAL** | **80 / 100** | **GO** |

**Falsifier observable**：若 compact index 不足以重建真實座標圖，至少一個 dataset 會缺 chrom/start/end、座標超出 GRCh38 長度，或 W_tree marks 無法閉合。

**Reality-test**：

1. 7/7 dataset 的 ideogram mark count 必須等於各自 W_tree。
2. 每個 mark 必須滿足 `1 <= start <= end <= GRCh38 chromosome length`。
3. 舊分類若無 current 等價 grain，必須顯示為 retired glossary，而不是移植舊 count。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---:|---:|---:|
| current INDEX 可支撐 coordinate ideogram | HIGH | KNOWN | 1 |
| GRCh38 autosome lengths 固定且可內嵌 | HIGH | KNOWN | 1 |
| 舊 7,143 母體可與 current W_tree 比較 | HIGH | KNOWN FALSE | 1 |
| single/linear/branched/star 可直接重算 | HIGH | UNKNOWN | 2 ⚠ |
| region-level CN 可當連續 CN segment | HIGH | KNOWN FALSE | 1 |
| 以 SVG marks 呈現 3.6k–9.7k regions 可在手機閱讀 | MED | UNKNOWN | 2 ⚠ |

右上象限先驗：current morphology grain、mobile density/interaction。

## §4 Quick Pilot / Checkpoint

1. 以 HCC1395 INDEX 生成 22 chromosome coordinate tracks。
   → 驗證：8,222 marks；0 越界；染色體寬度依 GRCh38 長度。
2. 加入 determinacy / evidence / HP / n_sSNV mode。
   → 驗證：每模式分組加總回 W_tree 或明示自己的分母。
3. 在 1440×1000、768×1024、390×844、320×720 檢查。
   → 驗證：body overflow=0；legend、labels、toggle、deep link 均可用。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---:|---:|---:|
| current morphology 的安全 grain 與分類公式 | HIGH | 1h | P0 |
| 7/7 coordinate/count contract test | HIGH | 1h | P0 |
| mobile ideogram density與點擊策略 | MED | 1h | P1 |
| current CN 只有 region-sidecar、非完整 segment track | MED | 0h（明示限制） | P1 |

## §6 Evidence Conflict Scan

| Prior conclusion | Relation | Source |
|---|---|---|
| 舊 `is_somatic` 23,810 / 35,332 backbone 已淘汰 | conflict：禁止搬數字 | `InterSubMod/docs/CURRENT_FOCUS.md` 2026-07-09 |
| current claim 只到 regional mutation-state candidates | dependent：禁止 clone/ancestry claim | `InterSubMod/docs/CURRENT_FOCUS.md` 2026-07-14 |
| 舊 topology workstation 已標 Archived / Deprecated | support：只可重用版型 | `build_topology_workstation.py:7` 與 generated banner |
| validated NEGATIVE 文件僅見 methylation hypothesis，與本 UI spec 無直接衝突 | neutral | `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md` |

Repo root 無 `MEMORY.md`，已記錄為 conflict-scan 限制；以 CURRENT_FOCUS 與 validated reports 補位。

## §7 Decision Path

- **Verdict**: GO。
- **Decision lock**: Y。
- **Next action**: baseline → current morphology audit → generator-first implementation → 7/7 build → Claude Code red team → Chromium full validation。
- **Revisit if**: canonical region schema、GRCh38 scope、primary HP contract 或 summary SHA 改變。

### Independent red-team gate

1. 最大風險：用熟悉的舊名稱包裝不同母體，造成數字看似可比。Gate：每圖固定顯示 grain + denominator。
2. 最大風險：8k marks 在手機縮成雜訊。Gate：保留全貌，手機以可水平捲動的固定繪圖寬度或 density aggregation，且 body 不溢位。
3. 最大風險：把 region-level CN 當完整 CN genome segment。Gate：CN 僅稱 region sidecar state。

紅隊未觸發降級；GO 維持。

## Provenance Footer

- **Base commit**: `2bec873`
- **Build time**: 2026-07-15 03:20 +08:00
- **Skill version**: `/pre-decision-audit` v0.1
- **Canonical summary SHA-256**: `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`
