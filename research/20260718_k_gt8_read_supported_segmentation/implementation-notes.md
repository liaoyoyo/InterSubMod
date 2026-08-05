<!--
建立時間: 2026-07-18 06:10 +08:00
目標: 追蹤 HCC1395 k>8 read-support segmentation 的實作決策、偏離、折衷與未決事項
處理範圍: research-only Python implementation、probe、HCC1395 chr1-22 full validation
關聯檔案:
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/pre-decision-audit.md
-->

# Implementation Notes：k>8 read-supported segmentation

## 設計決定

- 主算法固定為 ordered hypergraph block-reward DP，不使用 pairwise clique-expanded min-cut。
- 每條 unique canonical molecule 在單一 `HP family × exact known PS × legacy component` 內形成一個 constraint；只保留 fixed `R/A` calls，至少兩點才提供 linkage。
- 分塊為連續、非重疊、每塊 `k≤8`；同一 constraint 跨多 boundaries 仍只損失一次。
- tie-break：retained molecule weight → retained exact-pattern count → block 數少 → cut gap 總和大 → cut coordinates lexicographic。
- downstream VAF 只允許比較候選樹的相對可能性，不參與 segmentation。
- `log1p(molecule count)` 權重保留原始有限 Decimal，不以固定小數位量化；DP 依本次所有非負權重自動建立足以精確加總的 local Decimal context，避免加總順序改變末位與極近 tie 排序。
- full runtime 定義包含 chr1–22 的 sparse read extraction 與 k>8 segmentation；chr21 即使沒有 k>8 target 仍執行 extraction，避免 runtime 分母偷減。
- `k≤8` 與實體 span 是兩個正交限制：primary 報告保留只限制 `k≤8` 的 exact optimum；另以
  `50/100/200 kb` inclusive hard span cap 做 cached sensitivity，不把任一固定 bp 門檻宣稱為唯一生物真值。
- 區塊是否可進候選樹推論另設 evidence gate：只有
  `retained_exact_pattern_count>0 AND raw_retained_molecule_weight>0 AND primary_active_site_count>=2`
  才標 `TREE_READY_LOCAL`；其餘一律 `ABSTAIN_ZERO_LOCAL_EVIDENCE`。此 gate 只代表具備局部推論輸入，
  不代表真實或唯一演化樹。

## 偏離

- 2026-07-18 06:10：使用者要求 full HCC1395；由於即時 load/CPU=1.66 超過 gate 1.25，先完成程式與 probe，full launch 延後到資源通過。這是排程偏離，final scope 不縮小。
- 2026-07-18 06:24：首個真實 chr22 extraction 揭露 BAM 內有相同 canonical coordinate identity、不同 RG 字串的重複 alignment row。Frozen manifest 明定 `collapse_redundant_rows_with_identical_HP_PS`，但舊 extractor v1.2 未實作；新增 v1.3 extractor，只在 target evidence、BQ、MAPQ、HP、PS 全相同時 collapse，任一衝突即 fail-closed。
- 2026-07-18 06:27：首個 chr22 partition 在 log1p sensitivity 觸發 objective/disposition assertion。根因是 28 位 Decimal 不同加總順序造成 `1E-26` 差，不是切點或 constraint fate 錯。以 exact local context 修正並新增回歸測試。
- 2026-07-18 06:36：load1 降至 25.2，通過 load≤60 gate；啟動完整 HCC1395 chr1–22 sequential run。
- 2026-07-18 07:31：執行中 spot check 的 load1=67.42，`/big7_disk` 對應 sdc 在 1 秒
  `iostat -xz` 視窗為 100% util；full wall 因此是 shared-host、低優先序的實際 operational
  time，不是隔離式微基準。不得把與舊流程的比值稱為精確 speedup。
- 2026-07-18 07:48：獨立 red-team 找到 legacy `log1p` sensitivity 的 28-digit
  transcendental rounding 反例：兩側 exact product 都是 72，舊 Decimal log sum
  卻因約 `5E-28` 差選到不同 cut。Primary raw-molecule partition 不受影響；舊
  `weight_stable` 標籤降級，另以 exact integer product `∏(count+1)` 重算後才允許報告。
- 2026-07-18 09:22：portable reader 正式 verifier 在 classic-scrollbar Chromium
  iframe 揭露 `.analytics-top-bar width:100vw` 邊界；`clientWidth=1425`、`scrollWidth=1433`。
  新增 repo-local runtime wrapper，透過 plugin 官方 `runtimeHtml`／custom builder 擴充點只覆寫
  top-bar width/margins，不修改 plugin cache、不手改 artifact 或最終 HTML。

## 折衷

- 非重疊 blocks 可直接餵給現行 `k≤8` solver，但單一 read constraint 若涵蓋超過 8 個 site-span，任何此類分塊都必然有 loss。此類輸出保留作 `UNAVOIDABLE_READ_GT8/PARTIAL_CUT`，不偽裝 lossless。
- overlapping bags 可作 sensitivity，但只有每個 hyperedge 都完整落入至少一個 bag、running-intersection 與 separator candidate 相容時，才可稱 exact overlap；不作本輪 primary。
- read count 可能被 coverage 支配；正式結果比較 raw molecule、`log1p(count)`、equal-pattern 三種權重，切點不一致時標 `WEIGHT_SENSITIVE_UNRESOLVED`。
- 前述 `log1p(count)` robustness 必須採 exact product 等價比較；不得沿用 legacy
  Decimal-28 的 stable/sensitive 結果。
- DP 先最大化完整保留的 read-support，再以 block 少作 tie-break；因此實際 block 數可高於 `Σceil(k/8)` 理論下限。chr22 理論下限 15 塊、實際最佳解 16 塊，這不是 conservation failure。
- 加入實體 hard span cap 可阻止由 `≤50 kb` 相鄰邊串聯而成的超長 component/block，但會切斷一部分
  具直接 molecule 支持的遠距 linkage；因此應報 sensitivity curve，而不是默默固定一個 cap。
- 100% site assignment 只代表每個原 k>8 位點都有工程 fate；零 retained local evidence 的 block
  必須 abstain，不能以位點已被切入小區塊為由產生樹。

## 未決

- bridge threshold 1/2/3/5 只保留在 extractor diagnostics；本輪 segmentation constraint 直接使用 exact known-PS primary read patterns，不用 threshold 決定切點。
- 本輪 full 為了 runtime 分母與 frozen authority 完整，使用全 chr PASS sSNV，不改成 derived k>8 VCF。derived target-only extraction 只可作未來優化版，須另立 input identity。
- legacy tree solver adapter 另列第二階段；primary schema 保持 HP×known PS，禁止因 adapter 重新 pooling phase sets。

## 已完成驗證

### 單元、反例與 evidence gate

- ordered hypergraph DP、I/O、duplicate collapse、span-cap oracle、full summarizer、evidence gate、
  runtime／unavoidable constraint 分解、report artifact 與 cached sensitivity runner最終合計
  `95 passed`（Python）；scrollbar-safe portable wrapper 另有 `3/3 PASS`（Node）。
- 500 組小型 brute-force oracle：0 mismatch。
- 250 組含實體 span cap 的小型 brute-force oracle：0 mismatch；`cap=None` 200 組與 primary DP 完全一致。
- 3,574 sites＋50,000 synthetic constraints：2.71 秒，peak RSS 76,304 KiB。

### HCC1395 chr22 真實 probe

- extraction：543 sSNV；63,424 canonical eligible rows before collapse；9,660 identical duplicate identities collapsed；53,764 unique molecules；70,163 sparse site calls；所有 receipt checks PASS。
- extraction wall：2:56.31；peak RSS 141,676 KiB。
- k>8：6 components、98 sites；舊 densest-8 保留 48 sites、排除 50；新法保留 98/98，產生 16 blocks，全部 `k≤8`。
- raw read-pattern support：新法 `1,037/1,151 = 90.10%`；舊 densest-8 `787/1,151 = 68.38%`。
- weight stability：4/6 components 的 raw/equal/log1p cuts 相同；2/6 sensitivity。
- deterministic rerun semantic SHA-256 兩次相同：`77d6128cefb42902d4b3598fd08a11a46ce5eddd02bdf63fe1458ddb2f617297`。
- partition wall：1.50 秒；peak RSS 24,652 KiB。

### HCC1395 chr1／chr22 span-cap probe

- chr1 no-cap：27 blocks；raw read-pattern support `98.44%`；24/27 blocks 的實體 span >50 kb，
  最大 224,877 bp。
- chr1 50 kb cap：57 blocks；raw read-pattern support `91.92%`；相較 no-cap 降低 6.52 個百分點，
  但消除 >50 kb block。
- chr22 no-cap：16 blocks；raw read-pattern support `90.10%`。
- chr22 50 kb cap：29 blocks；raw read-pattern support `80.71%`；舊 densest-8 為 `68.38%`。
- chr22 evidence gate：12/16 blocks、73/98 block-sites 為 `TREE_READY_LOCAL`；其餘
  4 blocks、25 block-sites abstain。
- chr22 不可完整保存 constraint 分解：read-pattern 本身 `n_fixed_ra>8` 為 64 patterns／68 molecules；
  `n_fixed_ra≤8` 但 endpoints 橫跨超過 8 個 ordered sites 為 23 patterns／25 molecules。

### HCC1395 chr6 extreme-chain 中途驗證

- 83 個 k>8 components、25,657 sites 被守恆分配到 3,252 個 `k≤8` blocks。
- Frozen inventory 中 chr6＋chr16 合計 146/408 components（35.78%），但承載
  43,819/47,570 k>8 sites（92.11%）；全基因組平均必須搭配這兩個極端集中區分層解讀。
- raw read-pattern retention：新法 `65.88%`，舊 densest-8 `33.41%`。
- 只有 45 blocks／310 block-sites 通過 `TREE_READY_LOCAL`；3,207 blocks／25,347 block-sites
  必須 abstain。這個反例確認「切得下」與「有足夠局部 read 證據建樹」必須分開報告。
- chr6 不可完整保存 constraint 分解：read-pattern 本身 `n_fixed_ra>8` 為 827 patterns／845 molecules；
  `n_fixed_ra≤8` 但 endpoints 橫跨超過 8 個 ordered sites 為 150 patterns／154 molecules。
- chr6 partition child 的 pattern load／aggregate 為 10.871 秒；receipt 中歷史命名為
  `ordered_hypergraph_dp` 的 9.397 秒其實是 three-weight component loop，包含
  raw/equal/log 三次 DP、densest-8 counterfactual、diagnostic rows 與 aggregation，
  **不是 pure DP timer**。兩者都是 21.992 秒 partition stage 的內含子區段，不與
  stage wall 重複相加。

### Exact-log robustness remediation

- 新增 read-only exact-product audit：以任意精度整數 `∏(molecule_weight+1)` 精確比較
  legacy `Σlog1p(count)` sensitivity，primary raw-molecule cuts 只重現驗證、不改寫。
- chr22：legacy log 與 exact product 6/6 components 相同；corrected stable 4/6，
  stability 改判 0。
- chr6：legacy log 與 exact product 83/83 components 相同；corrected stable 82/83，
  stability 改判 0。
- Full 408-component remediation：8/8 gates PASS；legacy log 與 exact cut
  `408/408` 相同、差異 `0`；corrected stable `372/408`、sensitive `36/408`，
  stability 改判 `0`。wall `18.89 s`、RSS `78,948 KiB`。

### HP/PS observed-unit retention red-team

- chr6＋chr22 共 89 個 components；57 個 observed HP/PS units，其中 44 個符合
  `molecule weight≥20 AND exact patterns≥5`。
- 4,683 molecule×component incidences 中保留 3,364、cut 227、unavoidable 1,092；
  weighted retention 71.8343%。
- Eligible unit retention median 94.2868%，但 12/44 低於 50%；這證明 aggregate
  retention 可能被高支持 units 拉高，正式報告必須保留 unit-level distribution。
- 17 個 HP1/HP2 paired components 的 absolute retention delta median 1.25 pp；
  12/17 <5 pp，2/17 ≥25 pp，最大 38.26 pp。
- 沒有 observed constraint 的 70 components 不造假成 0% 或 100%；本 audit 不重建
  未被 partition outputs 識別的 HP/PS opportunity universe。

### Full run

- 輸入 manifest：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/release_contract/input_contract/canonical_manifest.json`
- 執行命令：
  `nice -n 10 ionice -c2 -n7 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/run_hcc1395_full_segmentation.py --manifest <manifest> --expected-manifest-sha256 16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45 --python /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python --output-root <full/HCC1395_chr1_22_v1>`
- 輸出：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1`
- 狀態：**COMPLETE**；outer exit `0`、`_SUCCESS` 與 runner receipt 已驗證。
- 正式輸出：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1`
- 完整 autosomal inventory：`79,687` PASS biallelic sSNV；`16,501` positional components
  （`8,279` singleton＋`8,222` multilocus）；其中 `408` 個 k>8 components、`47,570` sites。
- 新法輸出 `6,237` blocks，全部 `1≤k≤8`，`47,570/47,570` sites 指派完成；
  block-k 分布為 `41/59/63/98/110/162/324/5,380`（k=1…8）。
- 舊 densest-8 僅選 `3,264/47,570=6.8615%` sites，排除 `44,306`；這是 site
  assignment counterfactual，不代表舊流程只分析 3,264 個全樣本位點。
- 同一 BASEQ20 exact-pattern incidence denominator：
  - 新法 `52,281/58,812=88.8951%`；
  - densest-8 counterfactual `34,294/58,812=58.3112%`；
  - 改善 `17,987` incidences、`+30.5839 pp`、`1.5245×`。
- exact patterns：新法 `14,062/19,933=70.5463%`，舊法
  `8,578/19,933=43.0342%`。286 個有分母 components 中，255 改善、31 相同、0 較差；
  另 122 個沒有 observed denominator，不造假成 0%。
- `TREE_READY_LOCAL=682` blocks／`4,044` member-sites，其中直接 primary-active sites
  為 `3,972`；`ABSTAIN=5,555` blocks／`43,526` member-sites。ABSTAIN 仍含
  `1,062` active sites，不能解讀為「完全沒有任何觀測」。
- unavoidable 為 `4,088` patterns／`4,290` molecule weight：
  `n_fixed_ra>8` 為 `3,333/3,522`；`n_fixed_ra≤8` 但 ordered-site span>8 為
  `755/768`。另有可避免但被最佳分割切斷 `1,783` patterns／`2,241` weight。

### Full runtime 與比較界線

- 實測 outer wall：`8,198 s = 2:16:38`；peak RSS `625,892 KiB`；exit `0`。
- extraction：`8,058.586 s = 134.310 min`（outer 的 `98.299%`）。
- partition stage：`120.095 s = 2.002 min`；其中 pattern load/aggregate
  `87.586 s`，three-weight component loop `22.742 s`。後者含三組權重、densest-8、
  diagnostics 與 materialization，**不是 pure DP**。
- runner overhead：`19.319 s`。
- 舊 terminal-success v5 historical proxy：`5,086.484 s = 84.7747 min`；
  新／舊 operational wall index `1.6117×`，新 run 多 `51.8586 min`。
- 這不是演算法 speedup/slowdown：舊值是 filesystem timestamp proxy、BASEQ0、8,222
  multilocus regions＋tree/region/ledger；新值是 `/usr/bin/time`、BASEQ20、408 k>8
  components、沒有 candidate-tree/VAF，且 shared-host 負載不同。

### Full hard-span sensitivity

- no cap：`6,237` blocks、retention `88.8951%`、最大 span `224,877 bp`。
- 50 kb：`6,859` blocks、`84.4811%`（相對 no-cap `−4.4141 pp`），最大 `49,953 bp`。
- 100 kb：`6,336` blocks、`88.4666%`（`−0.4285 pp`），最大 `99,976 bp`。
- 200 kb：`6,238` blocks、`88.8934%`（`−0.0017 pp`），最大 `184,422 bp`。
- 工程建議：100 kb 是較均衡預設；50 kb 是較保守的短區塊 sensitivity；200 kb 幾乎等同
  no-cap。任何 cap 都不是唯一生物真值。

### Full HP/PS audit 與 portable HTML

- HP/PS audit：`237/237` gates PASS；723 observed units、528 eligible；
  overall `52,281/58,812=88.8951%`，eligible median `98.7977%`。
  258 個 HP1/HP2 paired components 中，208 的差異 <5 pp、233 <10 pp、8 ≥25 pp。
  HP3 不在本輪資料契約，不能推論 HP3/not-HP3。
- 最終 portable report：
  `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/report_delivery/20260718T014238Z_full_v7/report.html`
- canonical artifact：10 datasets／638 rows／30 manifest blocks／12 cards／8 charts／3 tables；
  Data Analytics validation `ok=true`。
- browser verifier：41 rendered blocks／12 metrics／8 charts／3 tables；
  validation/package/verification 全 PASS，source dialog PASS，1440/390 viewports PASS。
- HTML SHA-256：
  `f719b534ba5a02313b3d28451d2ae6f535bc7dceb607e1b39c79b0fd0c4f7efa`。
- Claim ceiling：以上證明 deterministic engineering conservation 與 read-supported local
  k≤8 segmentation 可執行；沒有執行候選樹 solver 或 VAF likelihood，不能宣稱唯一真實演化樹、
  clone/subclone 身分或 HCC1395 vs DORADO 跨來源拓撲重現。
