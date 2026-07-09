---
name: verify-workstation
description: >
  觀察驗證介面 — 把「一批需要人工逐項確認的結果」（loci / variants / regions / features / 任意 items）
  + 每項的 metrics + 觀察輸出圖（ISM read×CpG 甲基熱圖、read×read 距離矩陣等，或 inline SVG）
  生成「互動判讀工作站」standalone HTML：每項可記人工判讀（同意/存疑/否定 + 原因）存 localStorage、
  JSON/CSV 匯出/匯入再校正、⌖ provenance 來源 badge、修正過程紀錄(changelog) section、modal 細看大圖、
  §13-A 由構造防捏造（每個 metric 由 spec 注入、缺必填值即 refuse exit 3、分類須資料算出不可手填）。
  USE WHEN：「判讀工作站」「位點判斷與驗證」「逐位點/逐項人工確認」「review station」「per-locus/per-item judgment」
  「把結果做成可互動肉眼確認介面」「let me confirm/verify each locus/item」「想記錄我對每個結果的判讀並匯出」
  「補上修正過程紀錄的結果檢視頁」、需人類逐項 confirm/verify 一批結果並記錄判讀+輸出時。
  SKIP WHEN：單張方法解釋圖(用 /methods-example)、純 .md→HTML 排版無 per-item 判讀(用 /html-report-build)、
  AI 自動算的特徵 verdict 分析無人工 UI(用 /feature-layered-observation)、決策性散文報告(用 /results-report)、
  純檔案組織稽核(用 /data-audit)、單一結論看板(用 /research-dashboard)、程式碼層 build/test 驗證(用 /verification-loop)、
  假說 L1-L4 驗證層級規劃(用 /validation-protocol)。
metadata:
  category: 視覺化
  node_type: skill
paths:
  - "**/display*/**"
  - "**/*workstation*"
  - "**/*judgment*"
  - "research/**/scripts/*.py"
  - "docs/**/*workstation*.html"
---

# /verify-workstation — 觀察驗證介面（互動逐項判讀工作站）

> **一句話**：給一批要人工肉眼確認的結果，生成一個可逐項記錄判讀+匯出、嵌觀察圖、附修正過程紀錄、且數字由構造防捏造的互動 HTML 工作站。
>
> **定位**：`/methods-example` 產「解釋方法的一張圖」；`/html-report-build` 把「一份 .md 排版成 HTML」；本 skill 產「一批結果的互動人工裁決介面」。三者共用 §13-A 注入引擎但 deliverable 不同。

## 何時用（與鄰近 skill 的硬邊界）

| 你要的 | 用 | 不是 |
|---|---|---|
| 逐項人工判讀 + 記錄 + 匯出 + 嵌輸出圖 | **本 skill** | — |
| 解釋「方法怎麼運作」的單張圖 | `/methods-example` | 本 skill |
| 一份 .md 變漂亮 HTML（無 per-item 判讀） | `/html-report-build` | 本 skill |
| AI 自動算特徵 verdict（無人工 UI） | `/feature-layered-observation` | 本 skill |
| 程式碼正確性 / build / test | `/verification-loop` | 本 skill |
| 假說 L1-L4 驗證層級 | `/validation-protocol` | 本 skill |

## 8 步工作流 W0–W7

> **鐵則承襲 CLAUDE.md §13.0 / §13.7**：先有驗證過的數字才寫文件；產數字的分析與寫工作站的 generator **絕不同批平行**。

- **W0 LOCK-AND-GATHER**：分析先跑完 → 數字全落 `*.json/*.csv` → Read 讀回確認非 error/未完成。每個 metric 記下 `src`（檔案:行）。**物理隔離**：分析 Bash 與 generator 在不同 tool-call batch。
- **W1 CARD SCHEMA**：定義 per-item 卡片欄位 = `id` / `title` / `subtitle` / `badges[]`（PASS/Tier/✓驗/FP/CN…）/ `metrics[]`（每筆 `{k,v,src}`）/ `figures` / `meta_metrics`（卡頭一行顯示哪幾個）/ `modal_stats`（細看時的完整欄）。
- **W2 FIGURE BINDING**：每項綁觀察圖。**圖策略門檻**：≤ ~數百項 → `figures.mode="svg"`（inline，單檔可攜）；> ~數百項或真實資料點陣圖（如 60,700 張 ISM PNG）→ `figures.mode="png"` 外部相對路徑 + gitignore + 重生指令（見 `references/section_modules.md`）。
- **W3 COMPUTED CLASSIFICATION**：分類（CLEAN/CONFOUNDED、PASS/可能漏掉/篩掉、TierA/B…）**在 W0 資料層算出**寫進 spec，**generator 不重算、人不手填**。判準（如 normal FDR<0.05=confounded、PERMANOVA gate）要可由 src grep 到。
- **W4 VERDICT RECORDING**：`item_config.verdict_states` 定判讀態（預設 同意/存疑/否定）+ `reason_options`。generator 產 localStorage 判讀按鈕（key=`meta.lskey`）。
- **W5 PROVENANCE + CHANGELOG**：每 metric 帶 `src` → ⌖ badge。`changelog`（修正過程紀錄）= data-driven：`{data_status, binary_status, phase, corrections:[{id,what,status,effect,src}], audit_conclusion}`，每條 `src` 須親驗（commit 用 `git log` 確認 ancestor-of-HEAD、doc 確認存在）。狀態用 `in-HEAD/planned-phase2/not-done/authoritative/superseded`。`meta.build_commit` 用 `git rev-parse --short HEAD`。
- **W6 §13-A REFUSE-ON-MISSING RENDER**：跑 `tools/build_workstation.py <spec.json>`。`item_config.required_metrics` 列必填欄；缺 → **exit 3 refuse**（不 render dash，杜絕 page-04 Fig4 式 v1/v2 漂移）。後盾：`scripts/number_provenance.py audit`、`scripts/fill_report.py`。
- **W7 EXPORT + VISUAL QC（閉環，2026-07-09 落地工具）**：**跑 `python3 tools/render_html_shot.py <html> -o <png>` → 用 Read 工具看 PNG**（Claude=VLM，AI-Scientist-v2 VLM-feedback 模式、人工閘）→ 檢 tofu(CJK)/溢出/錯位/**卡片「無圖」**/無來源數字 → 修 → 重 render 直到 pass。工具自動報 pageerror + broken-image(naturalWidth==0)。🔴 **curl 200 ≠ 渲染 OK**（HTTP OK 抓不到「無圖」靜默降級，見 `docs/harness/HARNESS_FAILURE_LOG.md` H-001）。再確認 ⌖ 來源、判讀按鈕、匯出 JSON/CSV 可動；必要時 `number_provenance.py` 稽核報告級數字。

## 資產

- `tools/build_workstation.py` — 通用 generator（共用 case：≤ 數百項、inline SVG 或外部 PNG）。讀單一 `spec.json` → 出 standalone HTML。`req()` 缺必填值 exit 3。
- `templates/workstation_spec.schema.json` — spec 結構（meta / banner / changelog / sections / item_config / items / provenance）。
- `references/antifab_contract.md` — §13-A 注入/refuse/computed-class 規則。
- `references/section_modules.md` — 可選模組（funnel / 即時門檻試算 / 圖庫 filter / charts / SEQC2-confound）+ **genome-scale 密集陣列模式**（>~1000 項用 compact column array，見 HCC1395 ASM script 102）。
- `references/provenance_changelog.md` — changelog 資料契約 + src 驗證規則。
- `examples/` — 兩個 worked 實例（見下）+ `selftest_spec.json`。

## 兩個 worked 實例（兩種尺度）

1. **chr2:18M subclone（小尺度・inline SVG・富敘述）** — `InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/build_workstation_html.py` → `docs/explain/07_subclone-judgment-workstation-chr2-18M.standalone.html`。~10 項、inline SVG 圖、computed methyl class、conflict 校正表。本 skill 的 generator 即從此 pattern 抽出。
2. **HCC1395 ASM 篩選（genome-scale・外部 PNG・30,350 項）** — `InterSubMod/research/tsg_promoter_asm_reviewer/scripts/102_build_workstation.py` → `…/display_v2/20260610_judgment_workstation_01.html`。22-col compact array、60,700 外部 PNG、Tier A/B + 5 層驗證 + 即時門檻 + redesign meta-loop + §⓪b 修正過程紀錄。genome-scale 參考此 script（非本 generator）。

## 依賴（DEPEND ON，不吸收）

- `/html-report-build` 設計 token（暗色盤、sticky、card）。
- `/methods-example` 的 SVG primitive（當圖是示意而非真實資料點陣）。
- `scripts/{fill_report.py, number_provenance.py, provenance_stamp.sh}` — §13 機械後盾。

## 診斷

- generator exit 3 = §13-A refuse（缺必填 metric）→ 補 W0 資料或調 `required_metrics`。
- 圖空白 = `figures.mode="png"` 路徑相對於 HTML 所在目錄；外部圖需與 HTML 同層 `figs/`（gitignore→重生）。
- 判讀沒存 = 不同 `meta.lskey` 互不共用；換檔/換機請用匯出 JSON 再匯入。
- 建新 skill 改本檔後須同步 CLAUDE.md §3 計數（creation_guard / skill_registry_sync hook 會提醒）。
