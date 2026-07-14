<!--
建立時間: 2026-07-15 03:55
目標: 請 Claude Code 對 layered workstation 全基因分布與樣本全貌層做獨立 pre-implementation red-team review
處理範圍: layered-workstation-sample-overview
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/research/20260715_layered_workstation_sample_overview/pre-decision-audit.md
  - InterSubMod/research/20260715_layered_workstation_sample_overview/implementation-notes.md
-->

# Claude Code pre-implementation review prompt

你是本次變更的獨立 read-only reviewer。請先遵守 repo 的 `AGENTS.md`、`.claude/CLAUDE.md` 與 `docs/CURRENT_FOCUS.md`，將任務分類為 B（comprehensive validation），服務 G3/G4/G5。

## 目標

檢視在 `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py` 新增以下功能的資料契約、資訊架構、語義風險、手機可讀性與驗證方式：

1. 真正依 GRCh38 chr1–22 長度與 current `region_index[].start/end` 繪製的 coordinate-proportional SVG，全基因 scope 永遠保留。
2. 四種 region-level 色彩模式：C／Topo determinacy、read-state evidence、primary HP multiplicity、CN region sidecar。
3. 五個 sample-wide overview panels：Topo count、C count、determinacy、primary HP × H3、region n_sSNV。
4. 原始 JSON 連結只留在收合的 evidence drawer。
5. 保留現有 chromosome aggregate grid 作精確 count lookup。

## 已驗證的界線

- 舊版 HCC1395 的 7,143 / 6,288 / 35,332 與 single/linear/branched/star、A/B/C/E 都是不同且已退役的資料宇宙／ontology，不能搬入 canonical v5。
- current HCC1395 的主分母是 W_tree=8,222、W_primary=7,932；詳細真值只可由 current machine summary 與 canonical region-view 讀取。
- current stored display trees 對部分 analytical candidate set 只保留前 32 棵，因此不能由 HTML stored trees 完整回推 morphology family counts。
- current n_sSNV 最大為 8，因 producer 的 MAX_SNV=8；8 桶必須標成「8（含 cap 到 8）」而不是舊版 `>8`。
- CN 只能稱為 region sidecar state，不能畫成連續 whole-genome CN segment。
- Incomplete 是不可評估，不能歸為 ambiguous；Topo=1 不是 biological truth、confirmed clone 或唯一時間順序。

## 請檢查的檔案

- `research/20260715_layered_workstation_sample_overview/pre-decision-audit.md`
- `research/20260715_layered_workstation_sample_overview/implementation-notes.md`
- `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py`
- `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`
- `docs/methodology/_assets/layered_workstation/README.md`
- `research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`

避免直接讀取 67 MB 的 archived standalone HTML；舊 renderer 只需看 `build_topology_workstation.py` 的相關函式。

## 回覆契約

請只輸出繁體中文 review，不修改任何檔案、不執行 commit、不啟動長計算。內容必須包含：

1. `VERDICT: PASS / PASS_WITH_CHANGES / FAIL`
2. P0 / P1 / P2 問題清單，每項附檔案與行號或欄位證據。
3. 建議的 panel grain + denominator 表。
4. 建議的 SVG interaction/accessibility contract。
5. 7/7 build 與 Playwright 最小驗收矩陣。
6. 對「為何舊名詞不應直接恢復」的一段使用者可讀 copy。

請特別質疑：region midpoint 是否比 start tick 更誠實、8k SVG marks 的效能、手機是否需要水平捲動、色盲與非色彩辨識、active chromosome 是否會同步 ideogram 與 grid，以及圖表分母是否在視覺上足夠醒目。
