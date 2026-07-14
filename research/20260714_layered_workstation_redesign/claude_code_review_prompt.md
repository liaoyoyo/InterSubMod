<!--
建立時間: 2026-07-14 Asia/Taipei
目標: 提供 Claude Code 唯讀交叉審查 layered workstation 全站重構
處理範圍: index + 7 sample pages；canonical v5；desktop/mobile
關聯檔案:
  - InterSubMod/research/20260714_layered_workstation_redesign/pre-decision-audit.md
  - InterSubMod/research/20260714_layered_workstation_redesign/implementation-notes.md
-->

你是 InterSubMod layered workstation 的獨立資深科學資訊架構與前端紅隊。只做唯讀審查，禁止修改檔案、禁止執行 shell、禁止提出會改動上游研究數據的操作。

請完整閱讀：

- docs/methodology/_assets/layered_workstation/index.html
- docs/methodology/_assets/layered_workstation/README.md
- docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
- docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation.py
- docs/CURRENT_FOCUS.md 的 2026-07-14 最新 canonical 段
- research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md
- research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
- research/20260714_layered_workstation_redesign/pre-decision-audit.md
- research/20260714_layered_workstation_redesign/implementation-notes.md
- 現行視覺參考 screenshots：docs/methodology/_assets/workstation_visual_audit/11_after_layered_index_desktop.png、12_after_layered_hcc1395_overview_desktop.png、13_after_layered_hcc1395_region_desktop.png、15_after_layered_index_mobile.png、16_after_layered_hcc1395_mobile_top.png、17_after_layered_hcc1395_mobile_region.png

任務是對「全站重新大改」提出可實作、可驗收的獨立意見。使用者明示：

1. index + 7 sample pages 全部頁面格式都要改善。
2. chr1-22 全基因視圖必須保留並更有用。
3. network/topology 概念要更清楚、美觀、按重要性佈局。
4. raw .json 連結要隱藏於折疊 source/evidence drawer，不可打斷數字閱讀。
5. 桌機與手機都需 Chromium/Playwright 驗收。
6. canonical 主結果只能是 LongPhase-S recalibrated FILTER=PASS v5，ClairS PASS v6 只能 sensitivity。
7. 不可把 regional candidate tree 說成 confirmed clone/ancestry、不可把 read AF 說成 CCF/posterior、L3 正確狀態是 not_evaluated/bounded auxiliary、PS只作QC。

請輸出：

A. 現行 P0/P1/P2 問題（每項含 evidence selector/line、風險、具體修法）
B. index 的 section/order/component spec
C. sample page 的 section/order/component spec，特別列 whole-genome → chromosome → region → HP-family candidate network 的閱讀鏈
D. network 視覺 encoding contract（observed/inferred、forced/variable、candidate complete/incomplete、C vs Topo、primary vs auxiliary），避免只靠顏色
E. mobile/keyboard/accessibility acceptance criteria
F. 你最反對或質疑的 3 個設計假設
G. 明確 verdict：GO/PROBE/NO-GO，以及 first smoke 應先驗哪 5 件事

務必精確，不要泛泛 UI 建議；所有科學詞彙必受上述 claim ceiling 約束。
