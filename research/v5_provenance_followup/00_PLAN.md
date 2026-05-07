<!--
build_date: 2026-05-06
agent: plan mode → init-research scaffold (post-approval)
status: planned (P0)
priority: P0 (blocks PI report errata or confirmation)
related_audit: InterSubMod/docs/reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md
related_pi_report: InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
plan_source: ~/.claude/plans/nifty-enchanting-turtle.md (approved 2026-05-06)
-->

# V5 Provenance Followup — 專案計劃

## 為什麼有這個專案

2026-05-05 V5 Data Provenance Audit 發現 PI 報告 (4-29) 引用的 V5 數據實為 **Pass 1 only**（ploidy bug 讓 Pass 2 從未觸發）。`d0bcd8c` (4-30) 修復後 4-30/5-01 已產出 Pass 2 觸發的新 BAM，但**新 BAM 對比分析尚未完成**。本專案目的：補齊新 BAM 對比，回答「PI 報告 V5 結論是否站得住」+ 「self-phasing artifact 在 read-level 是否能找出具體個案 + 共通現象」。

## 核心目標（用戶 2026-05-06 確認的四層遞進）

1. **(找例子)** baseline self-phasing 在 read/位點層面找出**具體個案**證明問題真實存在
2. **(找共通)** 多案例間是否有**共通現象**（區域聚集、AF 分布、somatic density）
3. **(找原因)** 從個案與共通現象**推論機制**並用獨立證據鏈**驗證原因**
4. **(驗證修正)** V5 修正是否**真的修對了**？結果**合理可信**嗎？

## 子計劃概要（細節見 plan file）

| 階段 | 名稱 | 優先 | 估時 | 狀態 |
|------|-----|------|-----|------|
| **T1.1** | Three-BAM ISM Benchmark + Sanity Comparison | P0 | ~1 day | planned |
| **T1.2** | Read-level Vote countMap Audit + Regional Clustering（Path α）| P0 並列 | 1.5-3 days | planned |
| **T1.3** | 4-cell Ablation | P2（依 T1.1+T1.2 結果決定）| ~3 days | deferred |

## 三情境決策樹

T1.1 完成後依 ΔF1/ΔAUC/Δsanity 判定：
- **情境 A**（B vs A ≥ +5pp，B vs C ≥ +2pp）→ Pass 2 + V5 雙重加成 → PI 結論強化
- **情境 B**（B vs A 改善但 B ≈ C ±1pp 內）→ Pass 2 對 ISM verdict 冗餘 → 驗證用戶直覺
- **情境 C**（B vs A 差距縮小或反轉）→ V5 真實效益小於 PI 宣稱 → 需 PI errata

T1.2 完成後依 4 路徑（個案 trace / 區域聚集 / AF×density / 修正後消失）pass/fail 判定：
- ≥3 路徑通過 → 機制因果可確立
- ≤2 路徑通過 → 機制詮釋降級（從「self-phasing 是主因」降為「self-phasing 是部分原因」）

T1.3 觸發條件：
- 情境 A → 選做（補完整性）
- 情境 B → 強烈建議（量化 Pass 2 vs V5 flag 各自獨立貢獻）
- 情境 C → 必做（PI errata 證據）

## 關鍵檔案

詳見 [manifest.yaml](manifest.yaml) 與 [完整 plan file](~/.claude/plans/nifty-enchanting-turtle.md)。

**Input BAMs**（HCC1395 5kHz）：
- A baseline_09: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/baseline_09/tumor_tagged.bam`
- B v5_flag: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_tagged.bam`
- C v5_pass1only_4_12: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam`

**ISM binary**: `InterSubMod/build/bin/inter_sub_mod` (KDE-fixed, mtime 2026-04-21)

**待新建 scripts**：
- `InterSubMod/scripts/run_v5_followup_3bam.sh`（T1.1 wrapper）
- `InterSubMod/scripts/analysis/v5_provenance_3bam_compare.py`（T1.1 對比表）
- `InterSubMod/scripts/analysis/v5_provenance_vote_audit.py`（T1.2 投票分析）

**待修改 script**：
- `InterSubMod/scripts/analysis/v5_sanity_paired_check.py`（line 55-58 硬編路徑改 CLI argparse）

## 終點 deliverables

- `T1_1_3bam_compare/comparison_summary.tsv` + `comparison_report.md`（含三情境判定）
- `T1_2_read_level_audit/read_level_vote_audit.md`（含 4 路徑機制證據）
- 主報告：`InterSubMod/docs/reports/validated/2026/05/20260512_V5_pass2_validation_report_01.md`（建議命名，依 T1 完成日調整）
- evidence_ledger 條目 + memory 更新

## 啟動順序（建議）

1. **chr19 connectivity pilot**（30 秒）— 對 BAM B 跑 `run_vcf_all_snv.sh --mode chr19-verification`，驗 ISM binary 與 4-30 BAM 相容
2. **T1.2 binary patch + 編譯**（並行 T1.1，~半天）— 加 debug log + 編 3 版 testing-only binary
3. **T1.1 三 BAM 全基因組 ISM**（~36 min 順序跑）+ sanity check 參數化
4. **T1.2 chr19 投票 dump 三版本**（每版 ~3-5 min）+ 4 路徑分析
5. T1.1+T1.2 結果交叉判定 → 寫 comparison_report.md
6. 依結果決定 T1.3 是否啟動
