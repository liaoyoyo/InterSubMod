<!--
建立時間: 2026-07-14
目標: 對 canonical v5 layered_workstation 的 1 個 index + 7 個 dataset HTML 做 post-implementation semantic / visual red-team
處理範圍: Task B comprehensive validation；8 pages × desktop 1440×1000 / mobile 390×844；chr1–22 全量靜態重算
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/11_layered_workstation_ia_redteam.md；InterSubMod/docs/methodology/_assets/layered_workstation/；InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
-->

# Layered workstation post-implementation semantic red-team

用 SCQA（Situation → Complication → Question → Answer）：canonical v5 的資料 custody 與四個舊 P0 已關閉，但 2026-07-14 21:01–21:02 生成的 HTML 仍暴露數個 renderer／可及性矛盾；最新版 source 已修本輪找到的 renderer P1，唯有 region-level PS 查詢受上游 payload 限制，因此必須重生後再跑 final freeze QA。

**TL;DR：P0=0；資料與 C／Topo 主語意可接受，但本次 frozen HTML 不可直接發布；最新版 source 已修本輪 9 類 renderer／可及性問題，renderer P1=0，嚴格 §10 尚有 1 個 upstream acceptance gap（region-level mixed-PS membership），且尚未重生驗證（影響：中高，信心：高）。**

## 1. 任務分類、目標與 verdict

- **Task type**：B — Comprehensive validation；不是 subset、demo 或 exploratory pilot。
- **完整範圍**：`index.html` + 7 個 dataset pages；7 datasets / 6 biological samples；每頁 chr1–22；桌機與手機各一次，共 16 browser contexts。
- **服務目標**：G3（read-level epigenetic／evidence role 不越界）、G4（多 dataset 一致性與可重現）、G5（可被外部 reviewer 以 hash、selectors、receipts 重驗）。
- **稽核對象凍結點**：HTML 生成時間為 2026-07-14 21:01–21:02；後續 source 修補**不算**本批 HTML 的驗收證據。
- **發布判定**：**CONDITIONAL NO-SHIP UNTIL REGEN + FINAL QA**。P0 已全關，但 audited HTML 有 P1；source 與 artifact 尚未同版。

| Gate | 結果 | 證據摘要 |
|---|---|---|
| Canonical v5 custody | PASS | current summary SHA、7/7 region-view SHA、154/154 chunks、51,815/51,815 index/detail identities 全通過 |
| C／Topo／completeness | PASS | complete regions 全部 `C≥Topo≥1`；incomplete / no-primary 的 C、Topo 均為 unavailable |
| Primary／auxiliary | PASS with renderer wording fix | 72,994 primary units 守恆；H3/H4/none/reference 不進 W_primary、C、Topo |
| Partial-only | PASS | 44,672 primary units 保留；沒有被當作 no-data/no-clone |
| Forced-edge scope | Data PASS；audited wording FAIL | unit-level scope 重算正確，但 region-level incomplete 文案會遮蔽完整 sibling 的有效判讀 |
| PS／CN／L3 | PARTIAL | aggregate PS=5,623 且 0 PS-derived edge；CN/L3 不排序；region-level PS membership 不在 payload |
| Browser / responsive | PASS with P1/P2 | 16/16 擷取、0 page error、0 request failure、0 body overflow；1 favicon 404；舊互動焦點與假 tab model 失敗 |
| Raw JSON drawer | PASS | 8/8 頁皆只有 1 個預設收合 drawer；正文中 raw JSON links 不可見 |

## 2. 稽核 artifact 身分

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `InterSubMod/docs/methodology/_assets/layered_workstation/index.html` | 38,976 | `fb1442c70d65fd2ca6f5ababfa899a776b486a4e0644d01ca1bd9ef1102ced98` |
| `COLO829.html` | 37,263,161 | `376f86498d3c50c01269c5605307bbee22d7a301c8d169d3fb09932d2c21239b` |
| `H1437.html` | 42,687,770 | `ef58b898b53c18e3f74d39e2c51b843d8bc4ba1221a6da5732635cfdd8eef9c3` |
| `H2009.html` | 52,071,009 | `493ccac3363c277dd446f4ee0dfb79db54ddcf281fd037d60c9eb2cc5b13b32b` |
| `HCC1395.html` | 38,737,069 | `fad300c6b171900b572fd8e81d4276e04c74a71f256515ef0587ac649924915d` |
| `HCC1395_DORADO.html` | 40,356,899 | `e85c6095b6db6115a32fcb1136d96e19f00e30334816511b65ce881c57624fe4` |
| `HCC1937.html` | 13,743,826 | `80a482049956848fe58c44a1782ee5376a959e90a487aa45d8198a58197037be` |
| `HCC1954.html` | 21,206,310 | `d56d537260b925f1067c3796ac08a661bfa4de8b07f6377c2e01bbb51a652a49` |

Machine authority：`InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`，SHA-256 `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`。

## 3. Scientific / semantic acceptance

### 3.1 Canonical totals 與 scope

全量靜態重算與 current machine summary 一致：

| Metric | Recomputed |
|---|---:|
| W_tree regions | 51,815 |
| W_primary regions | 50,215 |
| No-primary regions | 1,600 |
| Primary units | 72,994 |
| Primary full + partial units | 28,322 |
| Primary partial-only units | 44,672 |
| Primary full-only units | 0 |
| Complete / incomplete regions | 42,240 / 7,975 |
| C=1, Topo=1 | 11,582 |
| C>1, Topo=1 | 10,737 |
| C>1, Topo>1 | 19,921 |
| Impossible C=1, Topo>1 | 0 |
| Mixed-PS regions | 5,623 |
| Recurrence-facet regions | 158 |

精確檢查：

- `#chunk-chr1` … `#chunk-chr22`：7 pages × 22 = **154/154** `data-sha256` 與實際 payload digest 相同。
- 每個 `.region_index[*].chunk_index` 對回該 chromosome chunk：**51,815/51,815** region identity 相同。
- 所有 complete primary units 由 stored full exact set 重算 unlabeled shape groups，與 `n_distinct_shapes_exact` 全部相同。
- 所有 complete region 均滿足 `C≥Topo≥1`；不可能狀態為 0。
- 所有 incomplete / no-primary region 的 joint C／Topo 都是 null／畫面 `—`，沒有把 unavailable 當 0。

### 3.2 Primary／auxiliary 與 partial-only

- Primary 僅 `is_primary_lineage===true` 的 mutation-bearing HP1／HP2；auxiliary H3/H4/none/reference 不增加 HP multiplicity。
- 44,672 個 partial-only primary units 仍有 candidate set；`n_full_pops=0` 不被推成「沒有 clone」或「沒有資料」。
- 可見 body 對 `softmax`、`posterior`、獨立 `VAF`、`clone-c` 的命中均為 0；`read ALT fraction` 只出現在「descriptive, not CCF」邊界說明。

### 3.3 Forced-edge scope 的反例重算

Region 是 incomplete 不代表每個 sibling unit 都 incomplete。全量重算找到：

| Dataset | Incomplete regions with complete primary sibling | 其中有可畫 forced edge的 regions | Forced edges |
|---|---:|---:|---:|
| COLO829 | 498 | 126 | 127 |
| H1437 | 947 | 272 | 322 |
| H2009 | 1,833 | 581 | 665 |
| HCC1395 | 345 | 145 | 158 |
| HCC1395_DORADO | 289 | 118 | 121 |
| HCC1937 | 21 | 16 | 17 |
| HCC1954 | 135 | 65 | 67 |
| **總計** | **4,068** | **1,323** | **1,477** |

Audited HTML 的 `#detail .verdict p` 把 incomplete region 一律說成「不評估 forced edge」，會否定完整 sibling unit；最新版 source 已改為「未完整 unit 不評估，完整 sibling 仍可在 unit scope 評估」。資料本身的 unit-level forced scope 沒有錯。

## 4. 兩個強制特殊案例

### 4.1 H2009 recurrence + capped overlap

Filters `#ftopo=incomplete` + `#fsignal=recurrence` 精確回傳 **4** regions；再填 `#fq` 精確座標回傳 1。四區全部：candidate incomplete、C/Topo null、recurrence facet=true，而且不進 complete denominator。

| Region | Recurrence primary family | Capped primary family | Stored / reported prefix |
|---|---|---|---|
| `chr8:79992384-80149786` | 2 | 1 | 1/1；16/16 |
| `chr9:275701-337149` | 2 | 1 | 1/1；2/2 |
| `chr13:93837736-93888639` | 1 | 2 | 3/3；1/1 |
| `chr15:31733893-31800487` | 2 | 1 | 2/2；7/7 |

Audited `chr8` detail 有 4 SVG、14 條 ROOT edges；舊 `edgeAcquisition()` 對 ROOT 回空字串，而且 `tree.recurrence` 未進 exact SVG，因此可見圖與 accessibility tree 都無法說清 repeated acquisition。最新版 source 已增加 ROOT→state 的 `+S_i`、directed arrow、edge adjacency/status description、candidate-level recurrence note 與 exact-edge repeated-site 標記。

### 4.2 HCC1395_DORADO no-primary

- `#ftopo=no_primary_lineage` 精確回傳 **720** regions。
- 目標 `chr1:190064024-190196077`：Primary lane=0、Auxiliary lane=2、candidate status=not applicable、C/Topo=`—/—`。
- 舊 `#detail .site-table caption` 固定顯示「primary HP lineages 加總」，但 `siteEvidence()` 其實 fallback 到 auxiliary／control，形成直接誤標。
- 最新版 source 已讓 caption 明示「no-primary region：auxiliary／control units 加總；僅描述觀察，不進 C／Topo」。

## 5. Audited HTML 的 P1 renderer findings 與 source disposition

| ID | Audited artifact 證據（exact selector / count） | 風險 | 最新 source disposition |
|---|---|---|---|
| P1-01 no-primary caption | `#detail .site-table caption`；DORADO **720** no-primary regions | auxiliary 證據被誤稱 primary | 已修，待重生 |
| P1-02 fake tabs | `.shape-tab,.tree-tab[role=tab]`；`[role=tab][aria-controls]`=0、`[role=tabpanel]`=0、explicit tabindex=0；ArrowRight 不動 | 無效 ARIA tab grammar | 已改為 ordinary button group + `aria-pressed`，待重生 |
| P1-03 incomplete region wording | `#detail .verdict p`；4,068 mixed-completeness regions、1,477 valid sibling-unit forced edges | 把 region-level incomplete 誤套到所有 units | 已修 scope wording，待重生 |
| P1-04 network accessible grammar | `#detail .network-scroll svg` 的 aria snapshot 只有圖名／節點／部分 `+S`，無 directed adjacency/status；H2009 target 有 14 ROOT edges | 看不出方向、邊穩定度與重複 acquisition | 已加 arrow、ROOT `+S_i`、sr-only edge list、recurrence exact annotation，待重生 |
| P1-05 focus loss | `.result-row` Enter 後 7/7 desktop `document.activeElement=BODY`；多 tree 案例點 `.tree-tab` 後也為 BODY | keyboard user 失去位置 | 已修 detail/result/tree/Prev/Next focus route，待重生 |
| P1-06 capped prefix labels | **8,014/8,021** capped units 顯示 `Stored / total x/x`；只有 7 個 stored<reported | x/x 被誤讀成 exact universe 完整 | 已拆 Exact／Stored exact／Prefix 以及 exact-total／enumerated-prefix 文案，待重生 |
| P1-07 recurrence visibility | Renderer 未消費 `tree.recurrence`；H2009 四區只能看到 region facet | repeated acquisition 只存在 JSON，圖上不可定位 | 已修 candidate note + exact SVG repeated-site，待重生 |

普通 button group 不再假裝 tab，因此 Enter／Space 可用且 ARIA 合法；原 §10 的 Arrow-key carousel checklist 改為 **N/A by component redesign**，不是宣稱仍為 tablist。

## 6. 最新 source diff-level review

快速凍結檢查的 renderer SHA-256 為 `09e4f9551a78bae87c30a5b5b1b955f3363e7a0e35b595fa041f660c61c5370a`；driver SHA-256 為 `e71885d903d6494496971d4570b5f33f49a7999bc9cb04f3f7eec264cebee727`。`python3 -m py_compile` 與抽出 inline JS 後的 `node --check -` 都以 exit code 0 完成。

最後一輪 diff-level red-team 另找到兩個 P1，主線已在上述 hash 修正：

- **Shape subset wording**：analysis-complete 但 display-subset 時，不再一律宣稱「unshown exact shapes remain」；現在以 `storedShapeCount===n_distinct_shapes_exact` 分流「所有 exact shapes 已涵蓋、只缺部分 exact trees」與「仍有未展示 shapes」。
- **Per-lineage accessible names**：`CANDIDATE_CACHE` 保存 `famLabel`；shape/tree groups、buttons、network scroll/SVG title、raw read group/caption 都加入 family + shape + candidate index，避免 HP1/HP2 元件同名。
- **URL history**：filters／region 使用 push/replace 分流；`restoreFromHash()` + `popstate` 回復 controls、active chromosome、selected region 與 detail。Refresh、copy URL、Back、Forward 的 source contract 已具備。

### 唯一仍開的 upstream acceptance gap — region-level mixed-PS membership

- Aggregate `mixed_PS_regions=5,623` 已在 index / dataset sidecar 顯示，且 PS 明確「QC only / not topology edge」。
- 但 `compact_index()` 無 PS membership，filters 也無 PS facet；source 目前誠實顯示：canonical region payload 沒有逐 region PS membership，因此刻意不提供 PS filter。
- 這避免了假資料與假 filter，卻仍未達原 acceptance「5,623 mixed-PS regions 可查」。要真正關閉，必須由 upstream region-view 增加 hash-bound `mixed_PS`／PS-count 欄位；renderer 不應自行猜測。
- 此缺口不造成現行 C／Topo、edge 或 filter 的錯誤，因此不是 P0，也不是可由 renderer 合理補猜的問題。

## 7. Visual / responsive audit

![七個 dataset 桌機 detail contact sheet](12_layered_workstation_postimplementation_redteam_assets/contact_sheet_desktop.png)

![七個 dataset 手機 detail contact sheet](12_layered_workstation_postimplementation_redteam_assets/contact_sheet_mobile.png)

![Index desktop](12_layered_workstation_postimplementation_redteam_assets/screenshots/index_desktop.png)

![Index mobile](12_layered_workstation_postimplementation_redteam_assets/screenshots/index_mobile.png)

Chromium 16-context 實際結果：

- 16 screenshots；HTTP status 200；page errors=0；request failures=0。
- `max(document.body.scrollWidth, document.documentElement.scrollWidth) <= innerWidth`：16/16，沒有 body-level horizontal overflow。
- Wide site/read tables 與 networks 都在有 `role=region`、accessible name、`tabindex=0` 的 local scroller；手機仍可到達右側欄位與圖。
- 8/8 raw evidence drawers 初始收合；可見 raw JSON links=0。
- 唯一 console error 是 index desktop 的 `/favicon.ico` 404；最新版 source 已加入 data favicon，待重生確認歸零。
- Index 的層級、scope funnel、C/Topo 比例、dataset launchers 與 evidence boundary 在 1440/390 都保持清楚；手機 cohort table 使用 local horizontal scroll，未裁掉整頁。

### P2 / final QA residuals

- 320px smoke 尚未在本輪 frozen artifact 執行；final freeze runner 必須補 8 pages × 320。
- 尚未以 NVDA／VoiceOver／TalkBack 實機朗讀；本報告只聲明 DOM／keyboard／accessibility-tree contract。
- HTML 單檔為 13.7–52.1 MB；雖已用 lightweight index + chr chunks 延後 JSON decode，仍需 final warm-server 量測 H2009 interactive median≤2s、filter p95≤100ms。
- `.detail-toolbar button` 40px、`.shape-tab/.tree-tab` 38px，符合 WCAG 2.5.8 的 24px AA 最低值，但未達原設計目標 44px；列為 touch polish，不影響科學語意。
- Index 小字 `#66777b` 在 canvas `#f2efe7` 約 4.07:1；位於白色 panel 約 4.60:1。Canvas 上的小型 muted text 應於 final visual polish 提升對比。

## 8. Raw JSON drawer 與 provenance

- Index：`details.evidence-drawer`=1、open=0、drawer 內 JSON links=8。
- 每個 dataset：`details.evidence-drawer`=1、open=0；audited artifact 有 3 個 source metadata links。
- 最新 source 另加入 hash-bound backbone sensitivity comparison JSON，但仍只出現在收合 drawer；不會污染主敘事或改 C/Topo denominator。
- Receipt hashes：
  - `static_semantic_receipt.json`：`ad55f939a9c1048a8409e32eda528e43bde821972e74c9f2b8c4988556b7b3ab`
  - `playwright_semantic_receipt.json`：`bda246f6768c8da04b7d085cbca1c19483757dcc5e6a93f4bc0e6db14adde86b`
  - Chromium capture script：`1ab8f312d71e1fd3266d429c4c33d1d9bab3bdf5d8ec96829adac5307986172f`

## 9. Step → Verify 與實跑命令

### 9.1 靜態全量語意重算

**輸入**：

- `InterSubMod/docs/methodology/_assets/layered_workstation/{index.html,7 sample HTML}`
- `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`

**命令**：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 /tmp/postimpl_static_audit.py
```

**輸出**：`InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/static_semantic_receipt.json`

**實際輸出摘要**：errors=`[]`；7 datasets / 6 biological samples；154/154 chunk hashes PASS；51,815/51,815 identities PASS；all C/Topo invariants PASS；H2009 overlap exactly 4；DORADO no-primary exactly 720。

### 9.2 Chromium / Playwright 16-context 擷取

**輸入路徑**：`InterSubMod/docs/methodology/_assets/layered_workstation/`

**執行命令**：

```bash
python3 /bip7_disk/liaoyoyo2001/.codex/skills/webapp-testing/scripts/with_server.py \
  --server "python3 -m http.server 8881 --bind 127.0.0.1 --directory /big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets" \
  --port 8881 -- \
  python3 research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/capture_postimplementation_semantic_audit.py \
  --base-url http://127.0.0.1:8881/layered_workstation
```

**輸出路徑**：

- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/playwright_semantic_receipt.json`
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/screenshots/`

**實際輸出片段**：

```text
PASS capture index desktop: index_desktop.png
PASS capture H2009 desktop: H2009_desktop.png
PASS capture HCC1395_DORADO mobile: HCC1395_DORADO_mobile.png
page_viewport_runs=16 screenshots=16 page_errors=0 request_failures=0 horizontal_page_overflows=0
```

### 9.3 Contact sheets

**執行命令**：

```bash
montage -thumbnail 480x333 -tile 2x4 -geometry 480x333+10+10 \
  research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/screenshots/{HCC1395,COLO829,H1437,H2009,HCC1395_DORADO,HCC1937,HCC1954}_desktop.png \
  research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/contact_sheet_desktop.png

montage -thumbnail 234x506 -tile 2x4 -geometry 234x506+10+10 \
  research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/screenshots/{HCC1395,COLO829,H1437,H2009,HCC1395_DORADO,HCC1937,HCC1954}_mobile.png \
  research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam_assets/contact_sheet_mobile.png
```

**實際輸出**：desktop `1000×1412`；mobile `508×2104`；已逐張視覺檢查。

## 10. Evidence limits 與最終交接

1. 本報告稽核的 HTML hash 已凍結於 §2；之後 source 修補不能回填為這批 artifact 的 PASS。
2. 沒有修改 HTML、generator 或 canonical data；本子任務只新增 audit report、receipts、capture script 與 screenshots。
3. `P0=0` 只代表 scientific/data contract 沒有發布阻斷；不等於尚未重生的 source 已完成 browser acceptance。
4. Final release 必須在一次性 regeneration 後固定 8 HTML SHA，再跑 1440/390/320、keyboard、URL history、no-network、chunk hash 與特殊案例 regression；若 SHA 再變，QA 證據失效。

**結論**：canonical v5 的核心展示方式已從「舊 clone/ranking UI」升級成可審核的 region → primary HP unit → exact candidate → topology shape → independent sidecars。最新版 source 的 renderer P0/P1 已歸零；嚴格 §10 唯一未關的是 upstream 沒有 region-level PS membership。所有 source 修正仍需由重生後 HTML 與固定 SHA 的 final QA 證明，不能以本報告的舊 artifact 截圖代替。
