<!--
建立時間: 2026-07-15 16:33 +0800
目標: 記錄七資料集拓撲比較與 HCC1395 技術配對 dashboard 的設計決定、偏離、折衷與 QA。
處理範圍: current canonical v5；GRCh38 chr1-22；7 datasets／6 biological samples。
cycle_id: cycle_20260715-1633-sample-topology-comparison
spec_id: sample_topology_comparison
status: validated
advisory: on
關聯檔案:
  - InterSubMod/research/20260715_sample_topology_comparison/pre-decision-audit.md
  - InterSubMod/docs/methodology/_assets/layered_workstation/index.html
-->

# Implementation Notes：七資料集拓撲比較與 HCC1395 技術配對

> **Purpose**：即時記錄比較 grain、統計判讀、資訊架構與驗證結果。
> **Task type**：B Comprehensive validation；全 chr1–22 × 7 datasets。

## 🔵 設計決定（Design Decisions）

### 2026-07-15 16:33 — 同時保留 7-dataset 與 6-biological-sample 口徑

- **Status**：Accepted
- **決定**：We will report dataset-level micro totals and a biological-sample macro profile that first averages HCC1395／HCC1395_DORADO, then gives each of six biological IDs equal weight.
- **理由**：避免同一 HCC1395 biological sample 因兩個 dataset 被重複加權，並保留使用者逐 dataset 比較需求。
- **Revisit if**：加入新的 technical replicate 或 biological replicate。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-15 16:33 — 相似性必須是多層證據

- **Status**：Accepted
- **決定**：We will separately show marginal composition TVD, exact-region overlap, matched label agreement／κ, and prior current-v5 site／allele-set／shape evidence; no single scalar will be called biological similarity.
- **理由**：比例接近不代表同一 genomic regions 得到相同結構；strict shape 也不能被 marginal agreement取代。
- **Revisit if**：未來有 cell-resolved ground truth 可定義 calibrated accuracy。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-15 16:33 — 事先固定 effect-size 門檻與 HCC rank 規則

- **Status**：Accepted
- **決定**：We will use the operational bands in `pre-decision-audit.md §4` and require HCC profile rank in the top 3 of 21 pairs before calling it relatively close.
- **理由**：防止看到數字後才移動判定門檻。
- **Revisit if**：累積足夠同樣本 technical-pair calibration cohort。
- **Evidence tier**：L2 ⭐⭐⭐⭐

<!-- BEGIN USER-SPECIFIED -->
**Decision**：展示必須能清楚整理各樣本實際拓撲狀況、種類比例、跨樣本比較，並驗證 HCC1395 兩份技術資料是否相似。
**DO NOT change**：不得只提供單樣本圖或只以視覺印象宣稱相似。
**Rationale**：2026-07-15 user request。
<!-- END USER-SPECIFIED -->

## 🟠 偏離之處（Deviations）

### 2026-07-15 16:33 — 「樣本相似」改成分層的「技術再現性」

- **Status**：Accepted
- **規範要求**：比較兩個 HCC1395 是否相似來驗證。
- **實作偏離**：介面會回答 composition／region／topology 各層是否接近，並稱 technical/cross-platform reproducibility；不稱 biological replication、accuracy 或 clone truth validation。
- **理由**：兩列來自同一 biological sample，且沒有 single-cell clone ground truth。
- **回退路徑**：若未來加入獨立 biological replicate 與 cell-resolved truth，可另建 accuracy layer。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

## 🟡 折衷考量（Trade-offs）

### 2026-07-15 16:33 — 主距離用 TVD，JSD/Hellinger只保留機器資料

- **Status**：Accepted
- **方案 A**：We will use TVD as the primary displayed composition distance because it has a direct percentage-mass interpretation; optional JSD/Hellinger may remain in the machine artifact for sensitivity.
- **方案 B**：只顯示 composite cluster；rejected，會遮蔽各 dimension 的差異。
- **理由**：7 個 dataset 的精確 pairwise lookup 適合 heatmap＋table，使用者可切換 structural／read-AF／morphology。
- **Revisit if**：樣本數增至可穩定做 cohort clustering。
- **Evidence tier**：L2 ⭐⭐⭐⭐

### 2026-07-15 16:33 — CI 以 chromosome block 而非 region IID

- **Status**：Accepted
- **決定**：We will bootstrap 22 autosomes with replacement for HCC pair matched-region metrics; region-level binomial CI 不作 headline。
- **理由**：同染色體區域不是獨立 biological replicates；block CI 較不會製造虛假精確度。
- **Revisit if**：可建立更合理的 spatial block／PS-aware resampling unit。
- **Evidence tier**：L2 ⭐⭐⭐⭐

## 🔴 未決問題（Open Questions）

### 2026-07-15 17:43 — sidecar exact-coordinate join 與各維度 HCC rank

- **Status**：resolved
- **Question**：兩側 primary-both exact-coordinate 母體有多大、三維 agreement／κ 為何、HCC在21 pair的距離名次是否落 top 3？
- **Resolution**：primary-both intersection=6,402、Jaccard=0.696；structural／read-AF／morphology agreement=69.79%／65.56%／75.18%，κ=0.520／0.446／0.557；三維 full-profile TVD=0.114、rank=2/21，conditional-evaluable TVD=0.085、rank=2/21。相對排名通過 top-3 規則，但 full profile 跨過 0.10 明顯差異門檻，因此只判定 partial／moderate technical reproducibility。
- **Priority**：P0
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

## ✅ 實作與驗證紀錄

### 2026-07-15 17:00 — 可重現分析與 machine artifacts

- **Input**：7 份 current-v5 sidecar、20260714 unknown-K artifact、HCC exact-signature artifact。
- **Command**：`python3 InterSubMod/research/20260715_sample_topology_comparison/scripts/analyze_sample_topology_comparison.py --bootstrap-iterations 5000 --seed 20260715`
- **Output**：`InterSubMod/research/20260715_sample_topology_comparison/artifacts/sample_topology_comparison.json`＋3 TSV＋validation receipt。
- **Result**：receipt PASS；7 datasets／6 biological IDs／21 pairs；21 個 dimension partition 全守恆；pair matrices 對稱；HCC confusion denominators 全為 6,402。

### 2026-07-15 17:28 — HTML dashboard 與原頁回歸修復

- **Implementation**：首頁加入 7-profile bars、三維切換、精確 table、7×7 TVD matrix、pair inspector、HCC agreement／κ／confusion／evidence ladder；JSON links 保持在 closed details。
- **Regression found**：首次全回歸發現 `cohort_read_af_morphology()` 回傳區塊被新函式插入位置截斷，首頁舊兩張 distribution cards 顯示為 `None`。
- **Fix**：將 morphology aggregate 與 HTML return 精確移回原函式，移除 unreachable dead block；重建後兩卡各自回加 W_primary=50,215。
- **Output**：`InterSubMod/docs/methodology/_assets/layered_workstation/index.html`，SHA256 `fce3b88e71c1bc3063da732556637408be301eb208120fed38a0d3899555fd49`。

### 2026-07-15 17:43 — Chromium／Playwright 最終收據

- **Focused comparison**：desktop＋mobile，30/30 checks，4 screenshots，exit=0。
- **Full workstation**：首頁＋7 sample pages × desktop/mobile，共 16 runs、219/219 checks、18 screenshots，exit=0。
- **Unit tests**：14/14 PASS。
- **Visual correction**：QA 截圖原先受 `scroll-behavior:smooth` 中途動畫影響；稽核內改用 deterministic scroll，並加入 7/7 profile row geometry 與 screenshot alignment assertions。全回歸另加入 Chromium captureScreenshot 有限重試，排除環境型 transient failure。
- **Durable report**：`InterSubMod/research/20260715_sample_topology_comparison/20260715_七資料集拓撲結構比較與HCC1395跨技術重現性驗證_01.md`。

## 📚 Lore

### 2026-07-15 — Sidecar 三維度同分母但不同問題

- **Constraint**：structural determinacy、read-AF selection 與 mutation-state morphology 都回加 `W_primary`，但不可相互代稱。
- **Why it matters**：Topo=1、read-AF unique-first 與 direct/sister morphology 分別回答候選完整性、描述性順位與 graph shape。
- **Evidence**：current-v5 sidecar schema 1.0.0。

### 2026-07-15 — Retained position 與 allele-site key 差一筆

- **Constraint**：HCC pair 以 `(chrom, position)` 交集為 25,121；以 `(chrom, position, REF, ALT)` 交集為 25,120。
- **Root cause**：`chr1:70439946` 在 HCC1395 為 C>A、DORADO 為 C>T；不是計算漂移。
- **Why it matters**：UI headline 採較嚴格 allele-site 25,120／84.87%，exact-signature position universe 只作 shape pairing；兩種 count 不混用。
- **Evidence**：兩份 SHA-bound `ssnv_site_ledger_*.tsv.gz` 的 retained rows 逐 key 比對。

## Provenance Footer

- **Commit hash**：pending；working-tree base `b95f19e7b0ffbbe9322ad91df7858d8c696a4036`
- **Build time**：2026-07-15 17:43 +0800
- **Skill**：`/implementation-notes` v0.1
- **Pre-decision**：`InterSubMod/research/20260715_sample_topology_comparison/pre-decision-audit.md`
- **Cycle**：`InterSubMod/state/cycles/cycle_20260715-1633-sample-topology-comparison/`
