<!--
建立時間: 2026-07-14
目標: 以 Chromium/Playwright 與逐張目視建立 layered workstation 改版前視覺、互動與可及性基線
處理範圍: 全部 8 pages、desktop 1440x1000、mobile 390x844、32 screenshots
關聯檔案: metrics.json, capture_layered_workstation_before.py, screenshots/
-->

# Layered workstation 改前完整視覺稽核

用 PREP：**8/8 頁的桌面展示穩定、資訊層級與互動完整；手機初始頁也能重排，但候選詳情的位點表與樹圖在 7/7 samples 都不可完整閱讀，因此 before baseline 判定為「可操作、不可直接視為 mobile-ready」。（影響：高，信心：高）**

## 1. Audit scope

- Task type：B comprehensive validation；`partial=false`，沒有 sample subset。
- Surface：cohort index + HCC1395、COLO829、H1437、H2009、HCC1395_DORADO、HCC1937、HCC1954。
- Viewports：desktop 1440×1000；mobile 390×844；device scale factor 1。
- Browser：本機 headless Chromium 147.0.7727.15，`file://` 靜態 HTML。
- Capture rule：每個 page×viewport 先拍 initial full-page，再做 DOM 解讀與互動；sample 另選真實 tree-bearing region。
- Evidence：32/32 PNG 已逐張開圖，沒有空白、loading、錯頁或裁切整張畫面的無效截圖。
- 研究目標：服務 G4 reproducibility 與 G5 external verifiability。

## 2. User goal and accessibility target

工作站應讓讀者從「cohort census → sample → region → family tree/read evidence」一路理解觀察、推論與限制；桌面與手機都應保留證據完整性、可操作性、鍵盤可達性與可辨讀的 network/tree。

## 3. Overall result

| 面向 | 結果 |
|---|---:|
| Pages | 8/8 |
| Page×viewport | 16/16 |
| Screenshots | 32/32 |
| Automated checks | 273 pass / 16 findings / 289 total |
| Runtime errors | console 0 / page 0 / request 0 / response ≥400 0 |
| Link targets | 156 tested / 0 missing |
| Main controls | 0 failed |
| Details toggle/restore | 76 tested / 0 failed |
| Heading hierarchy | 0 issue |
| Keyboard focus indicator | 0 missing |
| Source HTML hashes | 8/8 unchanged |

16 findings 是三個問題在不同頁/viewport 的重複證據：7 個 mobile table clipping、7 個 mobile tree readability、2 個 index whole-genome scope warning。

## 4. Confirmed strengths

1. **首頁先說界線，再說數據。** Hero、三張 interpretation cards、L0→L3 evidence stack、census 與 denominator protocol 的順序清楚；`PENDING`、`ambiguous`、`capped` 沒有被包裝成確定生物結論。
2. **7 sample pages 使用一致骨架。** Header、provenance、五個 summary metrics、進階統計、filters、region list、detail 的位置穩定，跨樣本切換成本低。
3. **主要互動全可用。** overview、copy-link、determinacy / HP / chromosome / sort filters、region search、row selection、sample switch、tree previous/next、thumbnail 與 details 全通過。
4. **鍵盤基礎良好。** Skip link 能顯示並到達 `#cohort-table` / `#regions`；目前可見 focus targets 都有 outline 或 shadow；沒有 positive tabindex。
5. **手機初始頁不是桌面縮小版。** Header、status chips、summary cards、filters 與 list 都能單欄/雙欄重排；全頁沒有 body-level horizontal scrollbar。

## 5. UX risks

### P0 — 7/7 手機位點證據表被裁切

在 selected detail 中，位點表實際寬 551.4–693.6px，起點約 x=43px；390px viewport 只看到 S1–S4 與下一欄的一小部分。`.detail` 的 overflow 行為讓 body overflow 仍回報 0，但資料右半部實際不可達，這是「沒有 body overflow」不能代表 responsive 的典型反例。

影響：讀者無法檢查全部 7–8 位點、ALT/VAF 與 recurrence，會直接破壞證據鏈。證據見 steps 7、12、16、20、24、28、32。

### P0 — 7/7 手機 network/tree 標籤縮到 3.9–6.9px

樹 SVG intrinsic width 為 380–684px；手機統一壓到 293px，scale ratio 約 0.43–0.77。最壞 H2009 的估計 label size 為 3.9px。`.tstage` 宣告可水平捲動，但 SVG 的 `max-width:100%` 先把圖縮小，所以實際 `stage_scrollable=false`；使用者沒有放大或橫向探索路徑。

影響：線條拓撲仍看得到，但 genotype、hidden ancestor 與 recurrence label 不可讀。專用證據見 step 8；其他 samples 由 steps 12/16/20/24/28/32 與 metrics 交叉確認。桌面保留 intrinsic size，但 9px label 仍偏小。

### P1 — 首頁缺少明確 whole-genome scope

首頁有 `Cross-sample inventory`、`eligible regions` 與 `Region-level denominator`，sample 也顯示 `染色體：全部`，但 index 沒有明確寫出 `全基因組 / chr1–22 / whole-genome`。對外讀者仍要自己推斷 census 的 genomic scope。

### P1 — HCC1395_DORADO 的 no-germline 文案互相衝突

同一 detail 同時呈現 `no-germline`、`HP×0 (single-HP)`、`1 germline-HP 家族` 與 `HP family · unphased(none)`。這不一定是資料錯誤，但目前 taxonomy 對讀者不一致，容易把「unphased/no germline family」誤讀成 single-HP。證據見 steps 22、24。

### P1 — 大型 standalone HTML 的本機載入已達可感知延遲

8 個 HTML 合計 143,834,193 bytes；H2009 單頁 38,274,656 bytes，本機 load + rendered UI 約 5.14–5.15 秒。其餘 sample 約 1.2–3.0 秒。這仍可離線使用，但不適合作為輕量快速巡覽介面。

### P2 — mobile index 的右側決策欄需更強提示

首頁 table 有水平捲動、固定 sample 欄與文字 cue，功能正確；但 390px 初始畫面只露出前 2–3 欄，Validation / L3 / Inspect 都在右側。現有 cue 沒有 edge fade、欄位摘要或「已到最右側」回饋。

## 6. Accessibility risks

- 已確認：visible heading hierarchy 無跳級；skip link、keyboard progression、focus indicator 全通過。
- 需改善：tree 沒有可見的 figure-level 摘要；手機視覺標籤過小，即使 DOM 中保留 node/edge `<title>`，觸控與鍵盤使用者也難以取得同等資訊。
- 觸控尺寸：跨 16 狀態觀察到 80 個 `<24px` 高度/寬度項目，主要是 inline links、breadcrumb 與 nested summary；其中 inline text 可有 WCAG 例外，不能直接宣稱整體違規，但 mobile summary 應優先擴大 hit area。
- 本次不能從截圖 alone 證明完整 WCAG compliance；未做實際 screen reader announcement、contrast ratio 計算、200%/400% zoom 或 forced-colors 測試。

## 7. Prioritized opportunities

1. **P0：替 detail 位點表建立真正的 mobile reading mode。** 最小修正是外包可聚焦的 local scroller、加入左右滑提示與 edge fade；較佳方案是 mobile 將每個 S1–S8 改成垂直 evidence cards，同時保留 desktop table。
2. **P0：樹圖不要 fit-to-width。** Mobile 保留至少 intrinsic/min-width、讓 `.tstage` 實際水平捲動，label 提升到 11–12px；另提供 zoom/reset、文字化 edge list 與 figure-level `aria-labelledby` / description。
3. **P1：在 index hero/provenance 增加 scope badge。** 建議直接寫 `Genome scope: chr1–22；primary denominator: eligible U3 regions；7 datasets / 6 biological samples`，避免 dataset、sample 與 genome scope 混淆。
4. **P1：統一 no-germline taxonomy。** `HP×0` 不應再標 `single-HP`，並釐清 `unphased(none)` 是 lineage bucket、family 還是 fallback tree。
5. **P1：降低首次解析量。** 可保留 offline artifact，同時把 region payload 拆成按 chromosome/region lazy-load 的 companion data，或先載 census/list、選取後再解析 detail。
6. **P2：改善 mobile cohort table discoverability。** 保留 sticky sample column，增加右側漸層/箭頭、欄位 selector，或在 table 上方提供每 sample 的 compact decision summary。
7. **P2：把 mobile nested summary / breadcrumb 的 hit area 拉到至少 24px。** 不必放大文字，只需增加 block padding。

## 8. Screenshot step ledger

| Step | Description | Health |
|---:|---|---|
| 01 | index · desktop · initial full page | 良好 |
| 02 | index · mobile · initial full page | 注意：table 需水平探索，whole-genome scope 未明示 |
| 03 | HCC1395 · desktop · initial full page | 良好 |
| 04 | HCC1395 · desktop · selected region browser | 注意：detail 密度高、tree label 偏小 |
| 05 | HCC1395 · desktop · expanded candidate tree | 注意：拓撲清楚，但 9px genotype label 偏小 |
| 06 | HCC1395 · mobile · initial full page | 良好 |
| 07 | HCC1395 · mobile · selected region detail | 需修正：630px 位點表右側不可達 |
| 08 | HCC1395 · mobile · expanded candidate tree | 需修正：608→293px，label 約 4.3px |
| 09 | COLO829 · desktop · initial full page | 良好 |
| 10 | COLO829 · desktop · selected region browser | 注意：detail 密度高、tree label 偏小 |
| 11 | COLO829 · mobile · initial full page | 良好 |
| 12 | COLO829 · mobile · selected region detail | 需修正：693.6px 位點表右側不可達 |
| 13 | H1437 · desktop · initial full page | 良好 |
| 14 | H1437 · desktop · selected region browser | 注意：detail 密度高、tree label 偏小 |
| 15 | H1437 · mobile · initial full page | 良好 |
| 16 | H1437 · mobile · selected region detail | 需修正：630px 位點表右側不可達 |
| 17 | H2009 · desktop · initial full page | 良好 |
| 18 | H2009 · desktop · selected region browser | 注意：detail 密度高、tree label 偏小 |
| 19 | H2009 · mobile · initial full page | 良好 |
| 20 | H2009 · mobile · selected region detail | 需修正：630px 表被裁切；tree 最小 label 約 3.9px |
| 21 | HCC1395_DORADO · desktop · initial full page | 良好 |
| 22 | HCC1395_DORADO · desktop · selected region browser | 注意：no-germline / single-HP 文案衝突 |
| 23 | HCC1395_DORADO · mobile · initial full page | 良好 |
| 24 | HCC1395_DORADO · mobile · selected region detail | 需修正：630px 表被裁切，且 taxonomy 衝突 |
| 25 | HCC1937 · desktop · initial full page | 良好 |
| 26 | HCC1937 · desktop · selected region browser | 注意：detail 密度高、tree label 偏小 |
| 27 | HCC1937 · mobile · initial full page | 良好 |
| 28 | HCC1937 · mobile · selected region detail | 需修正：551.4px 位點表右側不可達 |
| 29 | HCC1954 · desktop · initial full page | 良好 |
| 30 | HCC1954 · desktop · selected region browser | 注意：detail 密度高、tree label 偏小 |
| 31 | HCC1954 · mobile · initial full page | 良好 |
| 32 | HCC1954 · mobile · selected region detail | 需修正：630px 位點表右側不可達 |

## 9. Inline visual evidence

### Steps 01–02 · cohort index

![Step 01 desktop index full page](screenshots/01_desktop_index_full_page.png)

![Step 02 mobile index full page](screenshots/02_mobile_index_full_page.png)

### Steps 03–08 · HCC1395

![Step 03 HCC1395 desktop full page](screenshots/03_desktop_hcc1395_full_page.png)

![Step 04 HCC1395 desktop selected region](screenshots/04_desktop_hcc1395_region_view.png)

![Step 05 HCC1395 desktop candidate tree](screenshots/05_desktop_hcc1395_candidate_tree.png)

![Step 06 HCC1395 mobile full page](screenshots/06_mobile_hcc1395_full_page.png)

![Step 07 HCC1395 mobile selected region](screenshots/07_mobile_hcc1395_region_view.png)

![Step 08 HCC1395 mobile candidate tree](screenshots/08_mobile_hcc1395_candidate_tree.png)

### Steps 09–12 · COLO829

![Step 09 COLO829 desktop full page](screenshots/09_desktop_colo829_full_page.png)

![Step 10 COLO829 desktop selected region](screenshots/10_desktop_colo829_region_view.png)

![Step 11 COLO829 mobile full page](screenshots/11_mobile_colo829_full_page.png)

![Step 12 COLO829 mobile selected region](screenshots/12_mobile_colo829_region_view.png)

### Steps 13–16 · H1437

![Step 13 H1437 desktop full page](screenshots/13_desktop_h1437_full_page.png)

![Step 14 H1437 desktop selected region](screenshots/14_desktop_h1437_region_view.png)

![Step 15 H1437 mobile full page](screenshots/15_mobile_h1437_full_page.png)

![Step 16 H1437 mobile selected region](screenshots/16_mobile_h1437_region_view.png)

### Steps 17–20 · H2009

![Step 17 H2009 desktop full page](screenshots/17_desktop_h2009_full_page.png)

![Step 18 H2009 desktop selected region](screenshots/18_desktop_h2009_region_view.png)

![Step 19 H2009 mobile full page](screenshots/19_mobile_h2009_full_page.png)

![Step 20 H2009 mobile selected region](screenshots/20_mobile_h2009_region_view.png)

### Steps 21–24 · HCC1395_DORADO

![Step 21 HCC1395 DORADO desktop full page](screenshots/21_desktop_hcc1395_dorado_full_page.png)

![Step 22 HCC1395 DORADO desktop selected region](screenshots/22_desktop_hcc1395_dorado_region_view.png)

![Step 23 HCC1395 DORADO mobile full page](screenshots/23_mobile_hcc1395_dorado_full_page.png)

![Step 24 HCC1395 DORADO mobile selected region](screenshots/24_mobile_hcc1395_dorado_region_view.png)

### Steps 25–28 · HCC1937

![Step 25 HCC1937 desktop full page](screenshots/25_desktop_hcc1937_full_page.png)

![Step 26 HCC1937 desktop selected region](screenshots/26_desktop_hcc1937_region_view.png)

![Step 27 HCC1937 mobile full page](screenshots/27_mobile_hcc1937_full_page.png)

![Step 28 HCC1937 mobile selected region](screenshots/28_mobile_hcc1937_region_view.png)

### Steps 29–32 · HCC1954

![Step 29 HCC1954 desktop full page](screenshots/29_desktop_hcc1954_full_page.png)

![Step 30 HCC1954 desktop selected region](screenshots/30_desktop_hcc1954_region_view.png)

![Step 31 HCC1954 mobile full page](screenshots/31_mobile_hcc1954_full_page.png)

![Step 32 HCC1954 mobile selected region](screenshots/32_mobile_hcc1954_region_view.png)

## 10. Reproducibility and actual output

輸入：

```text
/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/layered_workstation/*.html
```

完整擷取命令：

```bash
python3 docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/capture_layered_workstation_before.py
```

重算命令：

```bash
python3 docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/capture_layered_workstation_before.py --reuse-screenshots
```

最終 actual excerpt：

```json
{
  "status": "fail",
  "exit_code": 1,
  "metrics": "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/metrics.json",
  "screenshots": 32,
  "pages": 8,
  "page_viewports": 16,
  "checks_total": 289,
  "checks_failed": 16,
  "source_hashes_unchanged": true
}
```

來源 SHA-256 前後一致；`git status --short` 只顯示新增 audit directory，沒有 layered workstation source modification。完整 per-page metrics、candidate region、control results、focus sequence、link targets 與 screenshot metadata 在 [metrics.json](metrics.json)。

## 11. Evidence limits

- 這是本機 `file://` artifact 稽核；load time 代表本機 parse/render，不代表 HTTP/CDN latency。
- 自動化實際操作了 control 與 keyboard path，但沒有宣稱完整 WCAG compliance。
- Tree readability 的 px 是 SVG scale 推估，再由逐張目視確認；不等同正式視力/使用者研究。
- 本輪依任務要求不修改 HTML、generator 或 data，因此所有問題保留為 before evidence。
