<!--
類型: Harness 自我失敗記錄 (self-harness failure archive)
建立: 2026-07-09
定位: ADAS「archive of discoveries」的失敗版(人工閘、restraint-compliant)。
      補 §8.3 Reflexion(研究方向失敗) + §9.2 Postmortem(NO-GO) 沒涵蓋的
      「harness/tooling/process 層失敗」— 讓下個 session 不重踩同一坑。
build_branch: research/subclonal-reconstruction-202606
-->

# Harness 失敗記錄（Self-Harness Failure Archive）

> **為什麼要夠豐富**：ADAS/AI-Scientist 的自我改進靠「archive of discoveries」；失敗記錄若只寫「壞了→修了」，下個 session 拿不到可行動資訊會重踩。**每條必含：症狀 / 觸發 / 根因(機制) / 最小重現 / 偵測訊號(含「什麼沒抓到」) / 修法 / 永久防線 / 連結。**
>
> **分工**：本 log = harness/tooling/process 失敗（可跨任務重現）。研究方向失敗 → `/scientific-rigor §8.3` REFLECTION.md；NO-GO → `docs/postmortems/`。
> **用法**：新 session 開工前掃本 log 對應類別；每次 harness/tool 踩坑後**當下**追加一條（趁 repro 資訊還在）。

## Schema（每條必填）
`ID | 日期 | 類別 | 症狀 | 觸發 | 根因(機制) | 最小重現 | 偵測訊號(含盲點) | 修法 | 永久防線 | 連結`

---

## H-001 · 2026-07-06 · [artifact-render] workstation 卡片「無圖」
- **症狀**：verify-workstation HTML 每張卡顯示「無圖」，圖不出現。
- **觸發**：spec 的 `figures` 欄我給成 list `[{path,caption}]`。
- **根因(機制)**：`build_workstation.py:123` 的 `figHTML` 期望 `figures={mode:'png', imgs:[{path,caption}]}`；拿到 list → `f.mode` undefined → 走 fallback `<div class="nofig">無圖</div>`。schema 不符不報錯、靜默降級。
- **最小重現**：spec item 加 `"figures":[{"path":"x.png"}]` → build → 開頁 = 無圖。
- **偵測訊號(含盲點)**：🔴 `curl -w %{http_code}` = 200 **抓不到**（HTTP OK ≠ 渲染 OK）；pageerror=0 也抓不到（非 JS 錯是靜默降級）。**只有 render→截圖→目視 / naturalWidth 檢查抓得到**。
- **修法**：`figures={"mode":"png","imgs":[...]}`（commit fe91c1e）。
- **永久防線**：**G1 閉環 QA**（`tools/render_html_shot.py` render→PNG→Claude Read 目視 + broken-image 偵測）；不再用 curl 200 當「圖 OK」證據。
- **連結**：commit ba32c0d/41a1401(加圖)→fe91c1e(格式修)；`tools/render_html_shot.py`。

## H-002 · 2026-07-06 · [artifact-render] 卡片無圖的前置 — 根本沒給 figures
- **症狀**：同 H-001 前一輪，卡片只有 metrics 無任何圖。
- **根因**：spec 完全沒 `figures` 欄（以為 metrics 就夠）。
- **偵測訊號**：目視「無圖」；`build_workstation.py:123` `if(!f)return '無圖'`。
- **修法**：對每候選 matplotlib 產 per-read 甲基分布 PNG → figs/ → 塞 figures。
- **永久防線**：資料密集觀察卡**預設要有視覺**（分布圖），不只數字；納入 G1 QA rubric「每卡有圖?」。
- **連結**：`gen_q5_figs.py`。

## H-003 · 2026-07-06 · [io-format] VCF tabix 用 .csi 非 .tbi → 靜默 0 結果
- **症狀**：`methyl_bimodality_cn_rate.py` 秒完成、`bimodal_rate_by_cn={}`（0 區）。
- **觸發**：`pysam.TabixFile(p)` 預設找 `p.tbi`。
- **根因(機制)**：本專案 filtered_snv VCF 只有 `.csi` 索引無 `.tbi` → TabixFile 拋例外被 except 吞 → ref_alt_map 回空 → 每區 `len(som)<2` skip → 0。**錯誤被 try/except 靜默**。
- **最小重現**：對只有 .csi 的 bgzip VCF `pysam.TabixFile(path).fetch(...)` → index not found。
- **偵測訊號(含盲點)**：🔴 exit 0 + 「跑很快」是唯一線索；不看內容會以為成功。**instant-complete + 空結果 = 紅旗**。
- **修法**：`idx = p+'.csi' if exists else p+'.tbi'`; `TabixFile(p, index=idx)`（methyl_auxiliary 的 `_tabix` 早有此處理，我重寫時漏了）。
- **永久防線**：新寫 tabix 存取**先驗一個已知有結果的 region fetch≠空**再全跑；「秒完成+空」必查。
- **連結**：`methyl_bimodality_cn_rate.py` tabix()。

## H-004 · 2026-07-06 · [output-schema] 加了記錄但忘了寫進輸出 JSON
- **症狀**：改 cnrate 加 `det_ct` 記錄，但重跑後輸出無 det_ct（merge 會空）。
- **根因**：`det_ct` 在迴圈有累加，但 `json.dump({...})` 沒列 `det_ct` key。
- **偵測訊號**：code review 當場抓到（重跑前）；否則 merge 才發現空、白跑一輪。
- **修法**：json.dump 補 `"det_ct": {...}`（改完 py_compile + 重跑）。
- **永久防線**：加新統計欄 → **同一次 edit 同時改「累加處」+「輸出處」**，改完 grep 輸出含新 key 再跑。
- **連結**：`methyl_bimodality_cn_rate.py`。

## H-005 · 2026-07-06 · [performance] 單進程 BAM 掃太慢
- **症狀**：cnrate 單進程跑 47min 未完（load 1）。
- **根因**：methyl_auxiliary 原本 15-way 分片平行；我單進程 = ~15× 慢。BAM 隨機讀 + GMM per-group 是瓶頸。
- **偵測訊號**：`ps etimes` 47min + 99% CPU 單進程 = 該分片。
- **修法**：加 `SM_CHUNK_TOTAL/IDX` strided 分片 → 8-way×3樣本=24 平行 → ~10min。
- **永久防線**：per-region BAM 掃**預設分片平行**（見 memory `reference_machine_resource_ceiling_parallelization`：批次要吃滿 44-46 核）；單進程跑久直接 kill+chunk。
- **連結**：`methyl_bimodality_cn_rate.py` main() chunk 段。

## H-006 · 2026-07-06 · [analysis-method] count 富集當因果（base-rate 假象）
- **症狀**：先說「77-86% residual 在 CN-gain = CN 假象」，後被自己的 rate 分析推翻。
- **根因(機制)**：只看 bimodal groups 的 CN 分布(分子)，無分母 → gain 區佔基因組大 → 大多數群本來就在 gain = base-rate 效應，非 CN-gain 造成更多雙峰。正確檢定 = 雙峰**率** by CN(帶分母)。
- **偵測訊號**：宣稱「X% 在 Y = Y 造成」時，問「Y 的 base rate 是多少?」缺分母 = 紅旗。
- **修法**：算 rate by CN(H2009 gain 0.170≈neutral 0.163=1.04x 平坦)→ 推翻。
- **永久防線**：`/scientific-rigor §4 DAG` + 「富集宣稱必附分母/base-rate」；納入 claim checklist。
- **連結**：memory `project_methyl_lineage_cn_gain_confound`；commit 1a23d74。

## H-007 · 2026-07-07 · [access] server 用 localhost 給遠端用戶
- **症狀**：用戶「HTML 沒有顯示」，但 server 8765 curl 200。
- **根因**：給的 URL 是 `localhost:8765`；用戶若在別機，localhost=用戶自己的機器非 server。
- **偵測訊號**：server 端 curl 200 但用戶看不到 = 網路/host 問題非檔案問題。
- **修法**：`--bind 0.0.0.0` + 提醒「遠端把 localhost 換 hostname/IP」。
- **永久防線**：給 server URL 附「遠端存取換 host」註；優先 SendUserFile(直接送檔) 免 host 問題。

## H-008 · 2026-07-09 · [orphan-figure] standalone dashboard 引用不存在的 figs 目錄（全破圖）
- **症狀**：`20260623_geometry_divergence_observation_dashboard.standalone.html` 每個 locus 卡片破圖。
- **觸發**：批次 render QA(48 交付物)flag 出來。
- **根因(機制)**：HTML 引用 `figs_geomdiv/geom_TP_*.png`(400 張外部相對 PNG)，但 **`figs_geomdiv/` 目錄根本不存在**(0 PNG) → 400/400 破。「standalone」名不副實(圖沒 base64 內嵌也沒隨檔 bundle)。
- **最小重現**：`render_html_shot.py` 該檔 → n_broken=400。
- **偵測訊號(含盲點)**：目視破圖 icon;`naturalWidth==0`。🔴 但**首輪批次只報 28/400**(見 H-010 lazy-load 盲點)。
- **修法**：未修(07-02 舊 pilot 非活躍 flagship;修需重生 400 張 geom 圖，待用戶決定)。
- **永久防線**：「standalone」HTML 交付前跑 render QA 確認 0 broken;圖要嘛 base64 內嵌要嘛同資料夾 bundle。`/pipeline-manifest` orphan-figure 稽核可涵蓋。
- **連結**：`docs/methodology/_assets/20260618_subcluster_pilot/20260623_geometry_divergence_observation_dashboard.standalone.html`。

## H-009 · 2026-07-09 · [tool-bug] batch_render_qa PNG 檔名碰撞(stem 相同覆蓋)
- **症狀**：批次 QA 中 `layered_workstation/H2009.html` 與 `topology_workstation/H2009.html` 產同名 `H2009.png` → 後者覆蓋前者。
- **根因**：`name = Path(h).stem` 只取檔名，跨目錄同名(各樣本)碰撞。
- **偵測訊號**：summary 有兩列 H2009 但 qa_batch/ 只一個 H2009.png。
- **修法**：`name = parent.name + '__' + stem`(namespace by 父目錄)。
- **永久防線**：批次輸出檔名含來源目錄。
- **連結**：`tools/batch_render_qa.py`。

## H-010 · 2026-07-09 · [tool-bug] broken-image 偵測漏 lazy-load(28/400 vs 真 400/400)
- **症狀**：figs_geomdiv 全缺(400 該破)但批次只報 28 broken。
- **根因(機制)**：`document.images.filter(i.complete && naturalWidth==0)` — fold 下 lazy-load(`loading="lazy"`)的 img `complete==false` 不算 → 只抓到首屏已嘗試載入的。networkidle 也不等 fold 下圖。
- **偵測訊號**：n_broken 遠小於視覺可見破圖比例 = 漏抓。
- **修法**：render 後 **JS 逐屏 scrollBy 到底觸發 lazy-load** → 再等 → 才數 broken(改後 400/400 正確)。
- **永久防線**：任何「全頁 img 完整性」檢查必先滾動觸發 lazy-load 再判。
- **連結**：`tools/batch_render_qa.py` scroll 段。

## H-011 · 2026-08-06 · [silent-corruption] SVG `<text>` 內用 HTML `<b>` → 其後元素全部掉出圖外
- **症狀**：手刻 SVG 圖只畫出前半，後半元素變成散落純文字堆在圖下方（數字黏在一起如 `30 8 6 25`、`HP1HP2`）。
- **根因(機制)**：HTML parser 在 SVG 命名空間內遇到 HTML-only 標籤（`<b>`/`<i>`/`<br>`…）會**中止 SVG 解析上下文**，該標籤之後的所有 `<rect>/<text>/<path>` 都被當成 HTML 文字節點處理。SVG 內要加粗必須用 `<tspan font-weight="700">`。
- **最小重現**：`<text ...>前段<b>後段</b></text>` 之後再放任意 `<rect>` → rect 不會被畫出。
- **偵測訊號(含盲點)**：🔴 **所有既有訊號全部失效** — curl 200、`no pageerror`、`no broken images`（因為根本沒有 `<img>`）、頁面高度可能不變（掉出的文字若在摺疊 `<details>` 內完全看不出）。**唯一可靠訊號 = 逐個 svg 區塊做 XML 解析，或實際截圖用眼睛看**。
- **影響量化(本輪實測)**：`docs/explain/02_ism-core` SVG#7「稀疏表 vs 密集表」**損失 16/32 元素(50%)**、SVG#2 損失 2/25(8%)；`04_subclone-reconstruction-chr2-18M` SVG#3 污染在最末元素故損失 0。**自 2026-06-12 建立起壞了近 2 個月無人發現。**
- **修法**：`<b>x</b>` → `<tspan font-weight="700">x</tspan>`；本輪並批次修 4 處。
- **永久防線**：`tools/explain_page_qa.py`（新建）—— 逐 svg 區塊 XML 解析 + 正則抓 `<text|tspan|desc|title>` 內的 HTML 標籤。交付前必跑。

## H-012 · 2026-08-06 · [font] emoji 在目標環境是豆腐方框（CJK 字型缺 U+1F534 等）
- **症狀**：頁面中 🔴 ⭐ ✅ ❌ 🟡 🟢 📌 全部顯示為空心方框，交付給教授的頁面出現大量豆腐字。
- **根因**：`Noto Sans CJK TC` 等字型不含 Emoji block（U+1F300–1FAFF）與部分 Misc Symbols（U+2B50 ⭐、U+2705 ✅、U+274C ❌）。**但 U+2714 ✔ / U+2718 ✘ / U+25CF ● / U+2605 ★ / U+26A0 ⚠ / ①②③ 皆正常**。
- **偵測訊號(含盲點)**：🔴 **grep 看不出來**（字元存在且合法）；render 無錯誤。只有**實際截圖用眼睛看**才發現。
- **修法**：HTML 內 🔴→`<span style="color:var(--c-dead)">●</span>`、⭐→★、✅→✔、❌→✘；**SVG 內不能塞 `<span>`**，改為直接移除（顏色語意已由 `fill` 屬性承載）或換等價符號。
- **影響範圍(本輪)**：解釋中心 17 頁中 **10 頁受影響**，最多的 `02_ism-core` 有 32 個。
- **永久防線**：`tools/explain_page_qa.py` 內建 `TOFU_RISK` 字典（含替代建議）。新增符號前先在目標環境截圖確認。

## H-013 · 2026-08-06 · [naming-lie] 檔名叫 `.standalone.html` 但實際 `<link>` 外部 CSS
- **症狀**：`docs/explain/index.standalone.html` 等頁移動位置或單獨寄給他人時**完全失去樣式**。
- **根因**：檔名宣稱 standalone，實際靠 `<link rel="stylesheet" href="_assets/explain.css">` 相對路徑。在原目錄開沒問題 → 長期未被發現。
- **偵測訊號**：只有把檔案複製到別處或單獨寄出才會暴露。
- **修法**：內嵌 CSS 全文，保留註解指明編輯 SoT 仍是 `_assets/explain.css`。本輪修 3 頁（index / 01 / 03）。
- **永久防線**：`tools/explain_page_qa.py` 檢查外部依賴（並先剝除 `<style>` 與 HTML 註解，避免內嵌 CSS 註解中提到 `<link>` 造成假陽性）。
