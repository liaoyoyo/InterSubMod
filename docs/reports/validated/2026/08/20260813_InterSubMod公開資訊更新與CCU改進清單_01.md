<!--
建立時間：2026-08-13
目標：說明 InterSubMod GitHub／Pages 公開資訊的本地修正、驗證與發布邊界，並整合 CCU 教學站唯讀改進清單
處理範圍：158 個公開 claim families、58 個 problem claims、17 個 Pages、37 個 inline SVG；CCU 32 個 frozen findings；不含 commit/push/deploy 或 CCU 修改
關聯檔案：
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/report_artifact.json
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/remote_state_receipt.md
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/20260813_CCU教學站重點改進清單_01.md
-->

# InterSubMod 公開資訊校正與 CCU 改進稽核

> **敘述框架：Claim–Evidence–Disposition ＋ SCQA。**
>
> **TL;DR：InterSubMod GitHub／Pages 的本地來源已完成 58/58 問題 claim 處置，17/17 Pages、37/37 SVG、68/68 browser checks 通過；但 live 僅 GitHub About 已更新，default main、Wiki、Pages 仍待另行授權發布。CCU 僅交付 16 項未解問題的改善清單，沒有改動 CCU。（影響：高，信心：高）**

本文件是 HTML 報告的 supporting technical source；主要交付為同名
`20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html`。HTML 的三張量化圖由 canonical artifact 原生產生 SVG，沒有另寫平行 chart runtime。

---

## 1. 完整判定

### 1.1 158 個公開 claim families

| Verdict | 數量 | 本輪解讀 |
|---|---:|---|
| `CONFIRMED` | 69 | 證據足以支持目前文字 |
| `CONFIRMED_WITH_LIMITS` | 31 | 保留明示限制後可用；不算待修問題 |
| `NEEDS_CORRECTION` | 28 | 必須限縮或補足分母／版本／假設 |
| `CONTRADICTED` | 26 | 現有證據直接反駁原敘述 |
| `UNVERIFIABLE` | 4 | 缺可重算 receipt，不得保留精確主張 |
| **合計** | **158** | **100 可保留；58 屬問題集合** |

58 個問題的 priority closure 為 `P0 34 + P1 20 + P2 4 = 58`。P0 registry 與 P1/P2 registry 的 claim-ID 聯集和 problem set 完全相等；兩個正向 guard 都 PASS，移除一個 claim 的負控制會 FAIL。因此這是 fail-closed correction contract，不是只挑容易修的 subset。

### 1.2 修正與 live 狀態不可混用

| Surface | 本地來源 | Live 狀態 | 還需什麼 |
|---|---|---|---|
| GitHub About | remote metadata | **`RESOLVED_LIVE`** | 保留 regression check |
| Default `main` README | 已修正 | 舊版；待發布 | 明確授權後 merge/push，再抓 live bytes |
| GitHub Wiki | 已修正 | remote head 仍舊 | 發布 Wiki 後核對 remote head/content |
| GitHub Pages | 17/17 已修正且 QA PASS | 17/17 routes 仍為舊 bytes | 部署後重驗 HTTP、hash、browser QA |

本輪沒有 commit、push、merge、PR、Wiki publish 或 Pages deploy；不能把 working-tree PASS 寫成使用者已在 live site 看到。

---

## 2. InterSubMod 的科學敘述上限

| 證據層 | 目前可以確認 | 目前不能確認 | 若要升級需要補 |
|---|---|---|---|
| Read／alignment | 同一 physical molecule 上的 allele、HP/PS、MM/ML called evidence；仍可能有 basecalling、mapping、phase/tag error | 來源細胞、同一 cellular clone、親源或實體 chromosome copy | cell barcode/single-cell 或 orthogonal truth、版本化 error/QC model |
| Genetic structure | 局部、允許 recurrence、模型條件下的 mutation-state candidate arborescence／topology signature | confirmed clone、lineage、完整演化史、canonical CN/LOH-corrected CCF | allele-specific CN/LOH、purity、multi-region/timepoint、held-out phase、truth calibration |
| Methylation | genetic-pattern-conditioned regional association；p/effect 依指定 null 與資料設計解讀 | 獨立確認 subclone、因果甲基化、functional mechanism、universal cell clock | 獨立標籤、held-out design、multiplicity/dispersion control、功能驗證 |
| 效能 | 綁定 commit、輸入、硬體、命令、重複次數與 exit 的版本限定測量 | 無 receipt 的線性 speedup、固定 latency 或跨版本泛化 | benchmark artifact、分布與失敗率、controlled scaling、binary/input identity |

公開入口現在統一使用：

> **physical-molecule evidence → local recurrence-allowed model-conditional mutation-state candidate → pattern-conditioned methylation association**

不再把它升格為 confirmed cellular clone／lineage／causal methylation。

---

## 3. InterSubMod 修正與確認清單

| 主題 | 發現的問題 | 已採取處置 | 終局判定 |
|---|---|---|---|
| 研究定位 | molecule co-occurrence 被寫成細胞／clone／lineage 直接觀測 | EN/ZH README、Wiki、Pages 改為 local mutation-state candidate；細胞層明列 model-dependent | `LOCAL_CORRECTED` |
| 分母 | `170,131` 對 component 與 record 分母混用；88.26% 被解讀為演化史 | 明列 170,131/255,752=66.52%；170,131/469,849=36.21%；63,506/71,955=88.26% 只限 rankable units | `LOCAL_CORRECTED` |
| 儲存量 | 1.67 TiB、287×、SEQ/QUAL >99% 缺 exact receipt/census | 舊精確主張撤回；只保留 7 sidecars `6,256,168,164 bytes = 5.826510641724 GiB` | `CORRECTED_WITH_GAP` |
| HCC1395 BAM | 舊 `292 GB` 不是 exact byte accounting | `stat -L` 為 `283,071,595,503 bytes = 263.63096712436527 GiB` | `LOCAL_CORRECTED` |
| 統計 | p、Cramér's V、PERMANOVA 被解讀成 truth、因果或 latent structure | p 限於 specified-null incompatibility；V 是 association magnitude；PERMANOVA 只測 label-associated centroid separation，需 PERMDISP | `LOCAL_CORRECTED` |
| Methylation | 用 mutation-defined labels 後再稱甲基化獨立確認 subclone | 只稱 pattern-conditioned association，且不能移動 genetic candidate topology | `LOCAL_CORRECTED` |
| LongLineage | HCC1395 0 topology units 被泛化；main/feature capability 混用 | 0 units 限 frozen HCC1395 dataset-gate；main `5daf50f` 與 feature `b9aaa12` 分開 | `LOCAL_CORRECTED` |
| Schema／版本數字 | 歷史 column counts 與 tracked 199-column runtime 混為單一全檔版本；測試/code counts 過時 | 歷史 schema 明示 scope；tracked source 199 columns；270 tests/39 suites、2,147/291 `.py` 綁日期、commit、corpus | `LOCAL_CORRECTED` |
| 效能 | 2.9 秒、<300 ms、32-core linear speedup 缺 benchmark artifact | 數字撤回；未來須有 hardware/commit/input/repetitions/distribution/failure receipt | `UNSUPPORTED_REMOVED` |
| Fail-closed generator | 「從根本防止捏造」超過程式實際能力 | 限縮為缺 declared required metric 時 exit 3；不驗 truth、denominator、method 或 science | `LOCAL_CORRECTED` |
| Pages HTML/SVG | 26/37 SVG 缺 title/desc；mobile overflow；page11–16 wrapper 不完整 | 37/37 SVG contract、17/17 structural PASS、68/68 profiles PASS；8 組 GitHub SVG/PNG 重建 | `LOCAL_QA_PASS` |
| Publication | local corrected 容易被寫成 live fixed | About 單獨 `RESOLVED_LIVE`；default main、Wiki、Pages 保持 pending | `1_LIVE / 3_PENDING` |

### 3.1 Pages 的可觀察驗證

- HTML：17/17 頁，title/h1/doctype/lang/ID/local-link contract 通過。
- Inline SVG：37/37 均含非空 `<title>`、`<desc>`、`role="img"`、有效 `aria-labelledby` 與 `viewBox`。
- Browser：17 pages × desktop/mobile/no-JS/print = 68 checks；HTTP 200、無水平 overflow、無 console/page/local request errors。
- 公開圖例：8 組 SVG/PNG；8/8 XML 可解析、8/8 title/desc、8/8 PNG 可解碼。
- 內容 guard：P0 34/34 與 P1/P2 24/24 均 PASS；P1/P2 negative control exit 1。

這些結果確認文件與呈現 contract，不會自動驗證生物 truth 或算法外部效度。

---

## 4. CCU 教學站：只列問題，不修改

### 4.1 32 項 closure

| 分區 | 數量 | 精確狀態 |
|---|---:|---|
| 既有 patch targets | 13 | `PATCH_VALIDATED_ON_PINNED_CLONE / NOT_APPLIED / NOT_DEPLOYED` |
| 仍待處理 | 16 | `OPEN_CHECKLIST / NO REPAIR CLAIM` |
| 先前已解 | 3 | `PRIOR_RESOLVED / REGRESSION_GUARD` |
| **合計** | **32** | **13 + 16 + 3 = 32；0 unclassified，0 duplicate** |

### 4.2 16 項下一輪優先清單

| 序 | 優先 | Finding | 核心問題 | 要補的資料／算法／統計 | 最低驗收 |
|---:|---|---|---|---|---|
| 1 | P0 | `OLD-P1-SR2C` | 64.89% 無 numerator/denominator，殘留 cell wording | artifact、unit/grain、version、SHA；找不到就移除 | 不再無 receipt 顯示 64.89%；不升格 cell |
| 2 | P0 | `DELTA-NEW-011` | n=60、δ=1% 直接推 0.07 detection bound | detector、H0/H1、α、power、sidedness、multiplicity、dispersion | 有版本化 power receipt，否則移除 0.07 |
| 3 | P0 | `DELTA-NEW-010` | H2009 max k=12 被寫成 19×成本原因 | matched profiler、k/coverage/unit/search/I/O、硬體/commit、confound analysis | 僅 controlled evidence 可稱 cause |
| 4 | P0 | `DELTA-NEW-009` | LOH 被寫成必然 quietly wrong | allele-specific CN/LOH truth、purity strata、loss simulations、abstention | 無絕對語氣；loss unknown 時 flag/abstain |
| 5 | P0 | `DELTA-NEW-007` | four-gamete 被說成唯一 tree／免 search | theorem scope、unique/non-unique/recurrence/loss/error fixtures、candidate search | 非嚴格假設不得稱唯一；五類 tests |
| 6 | P0 | `DELTA-NEW-003` | 未實作輸出被寫成可估 clone/CCF/cell proportion | current/proposed schema、CN/LOH/purity model、truth/calibration | persistent `SPEC / NOT IMPLEMENTED / NOT VALIDATED` |
| 7 | P0 | `OLD-P1-M06` | HP tag 型別與平台支援事實錯誤 | LongPhase-S/TO pinned release、runtime help、tag fixture | 每項由 fixture/command 重現；錯誤字串歸零 |
| 8 | P0 | `OLD-P1-M07` | 混用 LongPhase-S/TO tag，忽略本地 negative result | versioned benchmark、caller AF、TP/FP/FN、TP-loss guard | 工具語意分開；rescue 限定資料與版本 |
| 9 | P0 | `OLD-P1-SR6` | source 存在但 live 404／狀態不明 | route/navigation manifest、source/generated/live hash、deploy receipt | 正式 deploy 200+parity，或明示 draft 並排除導航 |
| 10 | P1 | `OLD-P1-M02` | 0.909/0.9093 缺 technical-pair receipt | dataset IDs、region、version、n、Spearman、CI | 不稱 accuracy／biological replication |
| 11 | P1 | `OLD-P1-M03` | 精確 F1 表缺 row provenance 與 TP/FP/FN | truth/BED/FILTER、confusion counts、command/hash | 每列可重算；缺件不顯示精確值／排名 |
| 12 | P1 | `OLD-P1-M05` | haplotype concentration 被當獨立確認 | marker-held-out、orthogonal anchor、CN/LOH/purity strata、controls | 未通 independence/systematics gate 時 abstain |
| 13 | P1 | `OLD-P1-M11` | CLI/threshold 被寫成跨版本事實 | release/commit、binary SHA、build flags、runtime help、benchmark | pinned binary 可重現；不稱 universal/optimal |
| 14 | P1 | `OLD-P1-SR2` | research guide 未持續標未實作／未驗證 | estimand、truth、pre-registration、implementation state | 首屏 persistent proposal banner；無 result tense |
| 15 | P1 | `OLD-P1-SR4` | algebra/spec 被讀成已校準 estimator | truth simulation、held-out phase、CN/LOH/purity、calibration | 首屏 `ALGEBRA ONLY`；不稱 validated/robust |
| 16 | P1 | `OLD-P1-GLOSSARY` | glossary 缺 recurrence/CN/LOH/purity/self-phasing 邊界 | terminology matrix、tag fixtures、authority mapping | 跨頁 fixture 一致；physical-copy/cell promotion 歸零 |

完整逐項的原問題、推論、必要資料、演算法／模型、統計、最低修文與 acceptance condition 見：
`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/20260813_CCU教學站重點改進清單_01.md`。

> **CCU 終局狀態：`NO CCU CHANGE / CHECKLIST ONLY / NOT PROOF OF REPAIR`。**

---

## 5. 驗證方法與限制

1. Frozen inventory 定義母體，避免修正後移動 denominator。
2. Registries 做 exact-set closure；required/forbidden anchors 加負控制，避免永真 guard。
3. 可重算精確數值保留版本、分母、命令；缺 receipt 不補猜測值。
4. README/Wiki/Pages 分流修改；Pages 另跑結構、可近用性、browser 與圖例 QA。
5. GitHub/Pages live state 以 API、remote refs、HTTP 與 byte hashes 獨立查核。
6. CCU 只讀 frozen TSV、receipt 與既有 patch；四個輸入在前後 bytes/mtime/SHA 不變。

限制：本輪沒有新增 SEQC2 truth-set performance、CN/LOH-corrected CCF、single-cell lineage、causal methylation、跨樣本 runtime 或 controlled speedup evidence。文件校正不等於方法性能已提升；About 之外的 live surfaces 也尚未更新。

---

## 6. 建議執行順序

1. 在乾淨 publication branch 重跑 P0/P1/P2 guards、Pages structural/browser QA 與圖例生成。
2. 取得明確授權後，分別發布 default README、Wiki、Pages；發布後重新抓 remote refs 與 17 頁 hashes。
3. 可選擇把 GitHub repository homepage 設為 Pages URL；這是 remote metadata action，本輪未做。
4. CCU 在沒有額外授權前維持 checklist-only；若授權，先修 9 項 P0，再補 7 項 P1。
5. 下一個研究補件應優先處理 CN/LOH/purity、independent/held-out phasing、truth-set benchmark 與版本化效能 receipt。

## 7. 最終判定

- **InterSubMod local source**：`PASS / CORRECTED`。
- **InterSubMod Pages local QA**：`17/17 pages, 37/37 SVG, 68/68 browser checks PASS`。
- **InterSubMod live publication**：About `RESOLVED_LIVE`；default main/Wiki/Pages `PENDING`。
- **CCU**：`CHECKLIST ONLY / NO CHANGE`。
- **科學結論上限**：molecule evidence、local model-conditional candidate、pattern-conditioned association；不升格 cellular clone/lineage/causality。
