<!--
建立時間: 2026-07-15 04:18
目標: 請 Claude Code 對 layered workstation GRCh38 與樣本全貌改版做唯讀終審
處理範圍: 7 datasets × chr1-22 × 4 Chromium viewports
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_workstation_ui.py
  - InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_sample_overview_after/metrics.json
-->

# Claude Code post-implementation read-only review

你是 InterSubMod layered workstation 的獨立終審 reviewer。請只讀，不要編輯、建立、刪除或格式化任何檔案，也不要執行長計算。任務服務 G3/G4/G5，task type B comprehensive validation。

## 審核輸入

1. `AGENTS.md`
2. `docs/CURRENT_FOCUS.md`
3. `research/20260715_layered_workstation_sample_overview/pre-decision-audit.md`
4. `research/20260715_layered_workstation_sample_overview/implementation-notes.md`
5. `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py`
6. `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`
7. `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_workstation_ui.py`
8. `docs/methodology/_assets/layered_workstation/README.md`
9. `docs/methodology/_assets/workstation_visual_audit/20260715_sample_overview_after/capture_sample_overview_after.py`
10. `docs/methodology/_assets/workstation_visual_audit/20260715_sample_overview_after/metrics.json`

不要讀取 30–40 MB 的 generated sample HTML；以 renderer、validator、metrics 與 README 做 contract review 即可。

## 必查問題

1. GRCh38 圖是否真用 chr1–22 長度與 `floor((start+end)/2)` 落點，而非等寬格或假座標？
2. 每個樣本是否能以 W_tree 一 region 一 mark 閉合；座標是否 fail-closed？
3. Topo、C、determinacy、Primary HP × H3、Region × retained sSNV 是否各自明示正確 grain 與 denominator，並在 build 時閉合？
4. `n_sSNV` 是否明確限制 2–8，且 8 說成 current cap boundary、沒有誤宣稱自然 8 與 cap-compressed 已可拆分？
5. 是否避免把舊 7,143 / 6,288 / 35,332、single/linear/branched/star、A/B/C/E 直接當 current canonical？
6. CN 是否只稱 region sidecar，而未包裝為連續 copy-number segment？
7. JSON 是否仍可追溯但預設藏在 closed evidence drawer，不干擾主要數字閱讀？
8. 五種全基因模式、legend isolation、keyboard chromosome sync、mobile local horizontal scroll 與 no-body-overflow 是否由獨立 Playwright assertion 實測？
9. freshness gate 是否同時綁 renderer SHA、UI contract 與 source mtime，避免舊 HTML 被誤判為新？
10. 有沒有 P0 科學誤導、數值契約、無障礙或互動問題；以及值得在 commit 前修的 P1？

## 輸出格式

請用繁體中文、精簡但具體，依序輸出：

- `VERDICT: PASS | PASS_WITH_CHANGES | FAIL`
- `P0`、`P1`、`P2`（沒有就寫「無」）
- `CONTRACT CHECKS`：逐項列出 1–10 的 PASS/FAIL 與一行證據（檔案/函數/metrics key）
- `REMAINING LIMITS`：只列真正仍未解決、且 UI 已否定誤讀的限制
- `COMMIT READY: YES | NO`

不得只說「看起來正確」；每個問題都要引用可觀察證據。
