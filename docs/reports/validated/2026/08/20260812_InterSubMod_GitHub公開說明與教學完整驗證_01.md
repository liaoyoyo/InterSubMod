<!--
建立時間: 2026-08-12 Asia/Taipei
目標: 完整驗證 InterSubMod GitHub README/Wiki/Pages 與 CCU lab-tutorial 現行公開教學的科學敘述、數字、演算法、CLI 與重現性
處理範圍: InterSubMod GitHub main@635437a、develop/feature@ddd8909、Wiki@6cfc990、Pages@fbdf7c7；CCU lab-tutorial deployed main@9eb1618；7 technical datasets/6 biological IDs/chr1-22 authority；fresh build/test/CLI；chr2:18M 獨立重算；外部原始論文
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/00_INDEX.md
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/algorithm_cli_matrix.tsv
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/verification_results.tsv
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv
branch: feat/lineage-tag-methylation-axes
commit: 73afaeac8e61c767241fa59c1ca6043a1c95290c
worktree: dirty; 僅新增本次稽核產物，未覆寫或清除既有使用者變更
status: VALIDATED_WITH_MAJOR_CORRECTIONS_REQUIRED
report_class: B_COMPREHENSIVE_VALIDATION
framework: SCQA + Claim-Evidence-Verdict
-->

# InterSubMod GitHub 公開說明與 CCU 教學站完整驗證

用 SCQA + Claim–Evidence–Verdict：**完整驗證完成 — InterSubMod 158 組去重公開主張中，58 組有實質問題或缺乏可驗證證據；現行 CCU 教學站雖已修掉部分舊問題，仍保留可直接反駁的工具敘述與分子→細胞升格，因此整體只能「有條件通過、需重大修訂」（影響：高，信心：高）。**

> **整體裁決：CONDITIONAL PASS — MAJOR REVISION REQUIRED**  
> 工程核心可建置、270/270 tests 通過，許多數值也可精確重算；但工程 PASS 不能替代 cellular subclone／lineage 真值。公開入口目前混合不同 commit、錯誤分母、不可攜命令與超出證據上限的生物敘述。

## 1. 30 秒結論

### Situation

InterSubMod 的確能整合 long-read sSNV、HP/PS、同分子共現與 read-level methylation；現行 authority 也支持一條可重算的 exact-PS funnel，以及局部、model-conditional 的 mutation-state candidate topology。

### Complication

公開 GitHub `main`、feature/develop、Wiki 與 Pages 並非同一版本；更重要的是，多處把「同一 molecule 的 linkage」寫成「同一 cell lineage」，把局部候選 topology 寫成已重建的 clone tree，或把 association-only methylation 寫成獨立 clone confirmation。另有可機械反駁的分母、schema、distance default、測試數、檔案數、storage ratio 與 quickstart 問題。

### Question

哪些內容可信、多少部分有問題、需要補什麼資料才能把說明升格為可對外防守的科學主張？

### Answer

- **158 組去重 claim families**：69 confirmed、31 confirmed with limits、28 needs correction、26 contradicted、4 unverifiable。
- **58/158（36.71%）** 至少需修正或補證據；其中 **54/158（34.18%）** 已可直接判定矛盾或需改寫，另 4 組目前不可驗證。
- **34 組 P0**：應先修 claim ceiling、版本入口、分母、演算法/CLI、公開重現性與 chr2 因果敘述。
- 以「每個公開 artifact 的最嚴重 claim」計，31 個現存 InterSubMod artifact 中 **29 個至少有一個 actionable 問題**，2 個只有需要限縮的 claim；另有 `main/README.zh-TW.md` 公開端點 404。
- CCU 教學站 current-delta 共追蹤 **32 項 findings**：17 OPEN、6 PARTIAL、6 REGRESSED、3 RESOLVED；也就是 **29/32 尚未完全解決**。這是 problem-focused finding denominator，不是「網站 90.6% 的內容錯誤」。
- 這個比例**不是「36.71% 的句子都是錯的」**：同一命題在英文、中文、Wiki、Pages 重複時只算一次；頁面判定只表示該頁至少有一個最嚴重問題，不代表整頁內容皆錯。

## 2. 稽核契約、版本與統計口徑

### 2.1 Task type 與服務目標

- Task type：**B — Comprehensive validation**；未以 1–3 chromosome 或 1–2 samples 代替完整母體。
- 服務 G2：區分 paired、tumour-only 與 read-level evidence 能支持的幅度。
- 服務 G3：限定 methylation 是 association-only sidecar，指出升級所需資料。
- 服務 G4：重算 7 technical datasets／6 biological IDs／chr1–22 authority 與 technical pair reproducibility。
- 服務 G5：以公開版本、fresh commands、hash、numerator/denominator 與原始文獻形成可外部重播證據鏈。

### 2.2 鎖定的公開版本

| 公開面 | 鎖定版本 | 實際狀態 |
|---|---|---|
| InterSubMod GitHub default `main` | `635437a65a33f8ba698acf85b22ebb069455c6cc` | 舊 README；首先進站者仍會看到過強 claim |
| InterSubMod remote `develop`／feature | `ddd8909a838318d8a77969313e9561c8ff9d01c2` | README/Wiki/Pages 的較新內容來源 |
| GitHub Wiki | `6cfc990f1fbe0b0aa76b8adbe56bbe645ef7ee7b` | 手動同步衍生物，不是自動同源輸出 |
| InterSubMod live Pages | deploy `fbdf7c7d229b0b37930778cbcd928cf9e15ecb17` | 17/17 HTTP 200；live bytes 與部署 source 相同 |
| 本地算法／authority | `73afaeac8e61c767241fa59c1ca6043a1c95290c` | dirty；build-input 與 remote feature tracked source byte-equivalent |
| CCU lab-tutorial live | `9eb1618d359e602d9c528675952b20d051fb2346` | 25 個正式 route HTTP 200 且 25/25 hash match；`sr6.html` 是 source-only draft／404 |

完整頁面母體與 SHA-256：[`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv)。

### 2.3 Claim 計數規則

1. 把可真偽判定的數字、絕對語句、演算法行為、CLI 行為與科學結論拆成 atomic proposition。
2. 英中鏡像、Wiki、Pages 的同義 proposition 合併為同一 claim family，但保留所有 occurrence。
3. 判定為 `CONFIRMED`、`CONFIRMED_WITH_LIMITS`、`NEEDS_CORRECTION`、`CONTRADICTED`、`UNVERIFIABLE`。
4. 「有問題」定義為後三類；`CONFIRMED_WITH_LIMITS` 不計入 58，但頁面仍應顯示 caveat。
5. 完整 158-row evidence matrix：[`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv)。

| Verdict | Claim families | 比例 | 意義 |
|---|---:|---:|---|
| CONFIRMED | 69 | 43.67% | 目前版本與 scope 下有直接證據 |
| CONFIRMED_WITH_LIMITS | 31 | 19.62% | 核心可保留，但須明列版本、grain 或 evidence ceiling |
| NEEDS_CORRECTION | 28 | 17.72% | 方向可能有價值，現行措辭超出證據或定義不精確 |
| CONTRADICTED | 26 | 16.46% | 可被 source/runtime/arithmetic/authority 直接反駁 |
| UNVERIFIABLE | 4 | 2.53% | 缺 receipt、benchmark 或系統性文獻證據 |
| **合計有問題** | **58** | **36.71%** | `NEEDS_CORRECTION + CONTRADICTED + UNVERIFIABLE` |

### 2.4 科學證據層級（L1–L5）

| Tier | 本報告中的證據 | 可支持的語氣 |
|---|---|---|
| L1 | frozen machine artifact、source symbol、fresh runtime、hash、numerator/denominator | 確認該版本與 grain 的直接事實；不自動支持 biological truth |
| L2 | 跨 technical datasets、held-out/LOSO、獨立重算、技術重複 | 支持再現／有限泛化；HCC technical pair 不是 biological replicate |
| L3 | validated authority/report 與多來源整合 | 支持 bounded integrated verdict；不得凌駕 L1 反證 |
| L4 | spec、proposal、toy model、future method | 只能用「設計／假說／預期／候選」，不能用 production/validated result 語氣 |
| L5 | 教學直覺、沒有 receipt 的精確數字或 novelty claim | 不可用作確認；應補來源或移除 |

## 3. 證據上限：最大問題不是程式不能跑，而是把工程證據升格成生物真值

```text
BAM / VCF / MM / ML
       ↓ L1：直接觀測
per-read allele、HP/PS、同分子 co-occurrence、methylation pattern
       ↓ L2：有模型條件的推論
bounded local molecular/mutation-state candidate + recurrence-allowed graph
       ↓ L2：conditioned association
genetic-label-associated regional methylation
       ╫ 目前缺口
allele-specific CN/LOH + purity/multiplicity/CCF + cellular barcode/truth
       ╫
confirmed cellular clone、biological lineage、whole-sample phylogeny、functional silencing
```

目前 authority 明列：confirmed cellular subclone = 0、confirmed linear biological ancestry = 0；methylation role 是 `pattern-conditioned association-only sidecar`。因此：

- 可確認：同一 read／molecule 是否同時承載特定位點的 allele。
- 可條件推論：在 exact PS、HP、bounded block 與指定 loss 下的局部 candidate topology。
- 不可直接確認：這些 reads 一定來自同一個 cell clone、其祖先關係、whole-sample phylogeny 或某基因已功能性關閉。

## 4. P0 修正清單：先做這 34 組

以下用「問題 → 反證／推論 → 最小修正 → 完成標準」整理；完整逐項 occurrence 與 source location 仍以 158-row TSV 為準。

### 4.1 Public source of truth 與 claim ceiling

- [ ] **C107／C108：首頁與 About 的 “subclonal reconstruction/resolution” 超過現行輸出。**  
  證據：authority 只允許 local recurrence-allowed mutation-state candidate，沒有 cellular truth。  
  修正：標題改為「single-molecule sSNV co-occurrence 的局部 mutation-state candidate reconstruction」；About 改為「read-level integration、local candidate analysis、methylation association」。  
  驗收：GitHub About、default README、feature README、Wiki 首頁、Pages 首頁使用同一 claim ceiling。

- [ ] **C109／C112：same molecule 被寫成 same lineage，candidate topology 被簡稱為 tree/reconstruction。**  
  推論：molecule ID 不帶 cell barcode；tree 是指定觀測模型與 loss 下的 candidate。  
  修正：所有首次出現處使用「直接觀測同分子共現；cellular co-membership/lineage 仍是模型推論」及「local recurrence-allowed minimum mutation-state candidate arborescence」。  
  驗收：全文搜尋 `cell lineage`、`clone tree`、`subclone`，每處要麼有獨立 truth，要麼明示 candidate／local／model-conditional。

- [ ] **C124–C131：chr2 教學把 allele association、normal-screen-negative methylation、approximate TVAF 與 rescue 升格成 clone confirmation。**  
  證據：新鮮重算確認局部分子共現，但沒有 cell truth；methyl groups 由 genetic label 條件化，不是獨立證據。  
  修正：`hard evidence`→`strong locus-specific association`、`tumor-acquired`→`tumor-associated/normal-ASM-screen-negative`、`subclone state`→`regional molecular-state candidate`、`correct tag`→`unvalidated future hypothesis`。  
  驗收：頁面 04–08 不再把 association 稱為 orthogonal cellular confirmation。

- [ ] **C135／C136／C140：default main 與 Project Summary 仍聲稱已解析 cellular subclones。**  
  反證：authority 的 biological claim counts 為 0；read clustering 不是 cell clustering。  
  修正：先修 feature 文案，再把 bounded README 原子合併到 main；cluster heatmap 改稱 read-level methylation clusters。  
  驗收：新訪客只讀 default branch 也不會得到「已確認 cellular clones」的結論。

### 4.2 Denominator、schema 與 algorithm contract

- [ ] **C113／C142：170,131 的母體被錯寫成全部 mutations。**  
  反證：`170,131/255,752=66.5219%` 是 k=1 strict components；`170,131/469,849=36.2097%`；current positional singletons 是 `50,432/469,849=10.7337%`。  
  修正：「k=1 strict read-linkage components：170,131/255,752（66.52%）」；禁止稱為「約三分之二的所有 mutations」。  
  驗收：數字旁固定顯示 numerator、denominator、unit/grain 與 artifact ID。

- [ ] **C059：`significance_summary.csv` 不是 193 columns／無 schema。**  
  反證：current runtime/source header 為 199 columns；包含 `VerificationSchemaVersion=2` 與 `RegionStratificationSchemaVersion=1`。  
  修正：版本釘到 `73afaea` 並寫 199；若需要 whole-file schema，另加單一 layout version。  
  驗收：文件數字由 header test 自動產生，不手抄。

- [ ] **C143／C144：effective default 不是 BERNOULLI。**  
  反證：ArgParser 會清除／重建 distance list；help 與 fresh minimal run 均顯示 NHD。  
  修正：「no-argument effective default=NHD；BERNOULLI 僅在明示指定時使用」；另註 Config initializer 不等於 CLI effective default。  
  驗收：01–03、Wiki、README 與 `--help` 一致。

- [ ] **C151：`inter_sub_mod` 不是 exact-PS 核心 funnel 的產生器。**  
  反證：core funnel 來自獨立 research exact-PS topology/AF solver 與 runners；`inter_sub_mod` 產 per-region methylation/statistics。  
  修正：公開 I/O 圖分開兩條 provenance chain。  
  驗收：每個數字都可追到實際 executable/script、command、input、output。

### 4.3 Public reproducibility 與 anti-fabrication

- [ ] **C134／C150：『每個數字、命令、格式都已驗證』的 blanket assurance 不成立。**  
  反證：測試數、Python file counts、schema、distance default、quickstart 與 storage claim 已出現矛盾。  
  修正：改為 per-claim receipt table，欄位至少含 commit、command、environment、input、output、date、scope、known failures。  
  驗收：文件 CI 讀 receipt，不允許沒有 artifact ID 的精確數字。

- [ ] **C133：report generator 不能讓 fabrication『結構上不可能』。**  
  反證：它能阻擋宣告為 required 的缺值，但 `CV=999999` 或不存在的 source 仍可被渲染；它不驗證科學正確性、分母或被省略欄位。  
  修正：「防止 declared required metrics 被靜默省略；不保證 claim truth」。  
  驗收：加入 range、source existence、hash、denominator 與 schema validators。

- [ ] **C138／C139／C149：公開 quickstart 不是 fresh clone 10 分鐘可跑。**  
  反證：BAM、FASTA、VCF 與 `one_snv.vcf` 4/4 必要 inputs 不在 Git；script 依賴 `/big7`、`/big8` 私有路徑。  
  修正：明標 `INTERNAL_ONLY`，或發布可授權 tiny fixture、checksum/download step 與全參數 portable wrapper。  
  驗收：空白環境 clone→build→fixture run exit 0，CI 保存 stdout 與 outputs hash。

- [ ] **C141：default branch 的中文 README 404。**  
  修正：英文與繁中 corrected README 同一 PR 原子合併。  
  驗收：兩個 raw endpoints HTTP 200，互相連結，claim-registry CI 均通過。

- [ ] **C148：`verify_pipeline_numbers.py` 沒有重算所有公開數字。**  
  反證：現行 validator 只涵蓋歷史 35,332-site pipeline，沒有重算 exact-PS、LongLineage、storage、code count、test count。  
  修正：把 validator scope 寫清楚，依 claim family 增設 validator。  
  驗收：158-row inventory 中每一個精確 numeric claim 具有 validator ID 或 `UNVERIFIED` 標誌。

### 4.4 chr2:18M 的邏輯反證

- [ ] **C152：LOH 不能推出每個 read difference 必為 somatic、也不能排除 sequencing error。**  
  證據：matched normal 沒有六個預期 ALT，但仍有少量其他 DEL/T errors；pos4 有 G/T/DEL ambiguity。  
  修正：「LOH 排除簡單的 two-parent-haplotype explanation；error、CN、stochastic methylation 等替代解釋仍需檢查」。  
  驗收：logic page 顯示 competing-explanation table。

- [ ] **C153：α/β zero observed co-occurrence 不是 sister branches 的唯一解。**  
  推論：finite coverage、block 間無 genetic bridge、recurrence allowance 及 missingness 都可能產生相同觀測。  
  修正：「在 parsimony 下與 α/β branching candidate 相容；非唯一識別」。  
  驗收：頁面同時列 support、counterexample 與 non-identifiability caveat。

- [ ] **C154／C155：methylation 不能直接說明基因被關閉，也尚未證明能 rescue 正確 tag。**  
  反證：沒有 expression/functional assay，亦無 per-read tag truth 或 quantified rescue benchmark。  
  修正：前者限於 local methylation state；後者移到 unvalidated future hypothesis。  
  驗收：沒有 observed-result 圖把 predicted/rescued tag 畫成 ground truth。

### 4.5 版本邊界

- [ ] **C043：『InterSubMod 與 LongLineage 都不能寫 tagged BAM』已被 branch-specific source 反駁。**  
  證據：InterSubMod 不寫 BAM；LongLineage feature `b9aaa12` 有 `longlineage-tag-bam`，public main `5daf50f` 沒有。  
  修正：句子必須釘 repo、branch/commit 與 binary。  
  驗收：feature 與 main 的 capability matrix 分列。

## 5. P1/P2 修正與補證據清單：24 組

| Claim | 問題 | 推論／直接證據 | 最小修正與驗收 |
|---|---|---|---|
| C047/C089 | 1.67 TiB tagged-BAM total 無 registry | current raw BAM 1.76384 TiB 是不同物件，不能代替 | 發布 7 個 tagged BAM 的 path、exact bytes、compression、hash、sum |
| C048 | 1.67 TiB／5.83 GiB 不是 287× | 以顯示單位計為 `1.67×1024/5.83=293.34×` | 由 exact byte receipt 重算；若沒 receipt 不報倍率 |
| C090 | `>99% payload、0% downstream use` 不可驗證 | 無 field-level byte census／utility audit | 只說「九欄 sidecar 不保留 SEQ/QUAL」，或補逐欄 census |
| C110 | short-read problem 被寫成普遍、可證明無解 | PairClone 可利用 proximal same-read pairs；single-bulk reconstruction 可在額外假設下 benchmark | 限定為「per-locus marginal VAF alone、無 linkage／extra assumptions」並引精確 theorem |
| C111 | single bulk 無法分開四種 methyl cause 被寫成普遍定理 | 是目前 measurement set 的 identifiability ceiling，不是全領域定理 | 加「在目前資料、無 orthogonal data/extra assumptions」 |
| C114 | eligibility flags 的 artifact 指向含糊 | 欄位在 `cohort_receipt.json`、`summary/all7_summary.json`，不在 authority top-level | 寫出 artifact 名、path、schema、hash |
| C115 | Wiki 被稱為 docs/explain 自動產物 | Wiki 是人工同步 derivative，可 drift | 更正文案並加 sync CI |
| C116 | 29 SVG 無 counting definition | deploy raw inline `<svg>` count=37 | 改成 37，或發布 29 個 semantic diagram registry |
| C117/C118 | `DEAD/never reopen/bulletproof` 過度絕對 | AGENTS 允許新假說重新設計 F1/TO/filter 實驗 | 改為「evaluated formulation 的多項 convergent negative tests；重開需新假說/pre-audit」 |
| C119 | paired caller 可「精準分開」somatic/germline | 所有 caller 仍有 TP/FP/FN | 改為「通常改善分離，仍須 truth benchmark」 |
| C120 | heuristic 結果稱 `true cis` | `1.8×drift+0.02` 是 HCC-only heuristic，不是 causal truth test | 稱 `cis-compatible/copy-screen-negative`，附公式與單樣本限制 |
| C121 | significant PERMANOVA＋gated V=0 稱 latent truth | PERMANOVA 是 label-associated centroid separation，且受 dispersion 影響 | 稱 candidate structure；加 PERMDISP |
| C122 | p「控制 truth」 | p 是與指定 null 的不相容度，不是 truth probability | 改正統計定義 |
| C123 | PERMANOVA「確認 grouping 有效」 | test 本身不驗證 biological labels | 改為測 centroid separation，並說明 design/dispersion limits |
| C132/C158 | first/unique novelty 未驗證 | 無 dated systematic search 或 matched head-to-head benchmark | 移除優先權 claim；或補 search protocol、納排標準、比較表 |
| C145 | tests/suites 已過期 | fresh run=270 tests/39 suites，全過 | 改 270/39 並釘 commit/date；最好由 command 產生 |
| C146/C147 | Python file count 錯 | deploy `fbdf7c7`：全 repo 2,147、`scripts/` 291 | 修數字並附 counting command、include/exclude rule |
| C156 | normal reads 全部 REF | 只確認六個 expected ALT 都未見；另有少量 non-reference errors | 使用較窄句子，不把「未見目標 ALT」等同「每條 read 全 REF」 |
| C157 | 32-core linear speedup、<300 ms 無 benchmark | 缺 hardware、inputs、repetitions、distribution、scaling curve | 移除，或發布 locked benchmark artifact |

## 6. 已確認、但要保留 scope 的核心清單

### 6.1 可確認的數據

- [x] **資料母體**：7 technical datasets、6 biological IDs、chr1–22、469,849 PASS biallelic autosomal sSNV dataset-records。HCC1395 與 HCC1395_DORADO 是同一 biological cell line 的 technical pair。
- [x] **strict components**：255,752；其中 k=1 為 170,131（66.5219%）。這個比例只屬 component grain。
- [x] **bounded funnel**：98,955 final groups；85,941 mutation-bearing；75,224 family-complete；71,955 ranked-complete。
- [x] **mathematical shape**：39,648 unique + 23,858 tied/same shape = 63,506；`63,506/71,955=88.2579%`。只代表 ranked、model-conditional rooted-unlabeled shape，不是 biological prevalence／accuracy。
- [x] **methyl formal assay**：1,045 formal、811 evaluable、3 robust；`3/811=0.3699%`。只支持 genetic-pattern-conditioned regional association。
- [x] **paired-pure read classifier**：baseline 0.861066、external 0.872226、ΔF1=+0.011159；95% bootstrap CI `+0.004437` 到 `+0.018808`；73,230 sampled reads／537 regions，屬 read-level PoC。
- [x] **canonical HCC1395 variant-level**：paired ΔF1=+0.000981；tumour-only ΔF1=+0.000274；不可與上列 read-level ΔF1 混用。
- [x] **technical reproducibility**：HCC1395/HCC1395_DORADO VAF Spearman≈0.909；支持技術再現，不是 accuracy 或 cellular truth。

### 6.2 可確認的工程行為

- [x] current build-input-equivalent source 可 Release configure/build，exit 0；編譯時間 242.81 s，有 4 類 warning。
- [x] direct GoogleTest：270 tests／39 suites，全過，exit 0；CTest：270/270，全過，exit 0。
- [x] one-locus fresh run：exit 0、1.90 s、85 reads、11 CpGs、NHD valid=3,443、invalid=127；因 reads<100 觸發 KDE unavailable／75× fallback，不能把此 run 當 KDE-corrected performance benchmark。
- [x] effective distance default=NHD；`--threads` help 顯示 1，但實際 Config fallback 為 16，這本身也應文件化。
- [x] LongLineage feature `b9aaa12` clean build 通過，47/47 tests 通過；`check_all.sh` 是 FOUNDATION_PASS，同時明示 release blockers 尚未解除，不能稱 production-ready。
- [x] Pages 17/17 HTTP 200 且 byte-match；feature/source raw inline SVG=37。

## 7. Fresh algorithm／CLI 驗證結果

完整 35-row matrix：[`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/algorithm_cli_matrix.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/algorithm_cli_matrix.tsv)。

| 稽核項目 | 實際結果 | 公開文字影響 |
|---|---|---|
| Config schema | 199 columns，含兩個 component schema fields | 193／no schema 是直接錯誤 |
| Distance selection | CLI effective default NHD；第一個 metric 驅動 primary clustering | 不能寫 BERNOULLI production default；多 metric 不等於共同驅動 |
| Multi-lineage axis | parser 可接多個，但 current logic 在第一個有效 axis 後 `break` | 文件須說 current behavior，不可暗示完整 multi-axis integration |
| MM/ML grammar | 現行支援範圍有特定 grammar／type 限制 | 不能泛稱所有合法 MM/ML encoding 均等價支援 |
| `tree.nwk` | read-cluster dendrogram／輸出結構 | 不是 cellular clone phylogeny |
| `methylation.csv` | 第一欄是 integer row id，非 read name | 文件 I/O schema 要改 |
| `subclone_structure` | deprecated region stratification | 不能當目前 cellular-subclone output |
| LongLineage tagged BAM | feature 有、main 無 | capability 必須 branch/version scoped |
| LongLineage real-data gate | HCC1395-only report topology units=0、production claim disallowed | 可說該 frozen gate fail-closed，不可泛化為所有資料或生物 NO-GO |
| 66,836 site 解釋 | 未成為 M1 stable assignments；co-occurrence artifact 仍有 79,687 site rows／134,278 pair rows | 「never produce a co-occurrence row」被直接反駁 |

## 8. chr2:18M 完整重跑：能支持什麼、不能支持什麼

輸入、命令、輸出與 stdout 收據見 [`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md)。Fresh output：[`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/chr2_18M_independent_audit_fresh.json`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/chr2_18M_independent_audit_fresh.json)。

| 觀測 | Fresh result | 可支持的最強句子 | 不可推出 |
|---|---|---|---|
| locus scope | chr2:18,066,481–18,110,828；單一 LOH interval | 區域內可做 bounded case audit | 全基因組／跨病人泛化 |
| 六點 truth | 1,2,4,5,6 確認；3 不在 HC、不可評估 | 五個位點有指定 truth support | 六點都同等可靠 |
| HKU rows | 340 records／280 unique／60 duplicate；HP2=232、HP1=4、none=44 | 讀段/標籤 composition 可重算 | cell-clone composition |
| strict α/β | violation=0 | 在觀測 coverage 下與分支 candidate 相容 | sister branches 是唯一解 |
| 1–2 cooccur | HKU 13；Dorado 6 | technical pair 均觀察到同分子 linkage | biological replication |
| 3–5 nested | HKU 4；Dorado 2 | 直接支持 3→5 的局部 nesting candidate | 完整 1→2→3→5 order |
| 4ALT–6 | HKU 10；Dorado 1 | 有局部 joint support | 穩健跨平台 prevalence |
| distant blocks | 無 coverage bridge | 不可直接連成一棵 tree | whole-region ancestry |
| methylation | 只有 3.1/3.2 是 clean replicated；其餘有 matched-normal ASM confound | genetic-label-associated local methyl pattern | acquired clone methylation／independent confirmation |

**整合判定**：此 case 支持至少三個 regional molecular-state candidates；不支持已確認 cellular clones、唯一全域 phylogeny、functional gene silencing 或正確 per-read rescue tag。

## 9. 外部原始研究如何反駁「過度絕對」的教學句

- [PairClone](https://academic.oup.com/jrsssc/article/68/3/705/7058355) 直接利用 short reads 上 proximal mutation pairs；所以「short reads 只有 marginal VAF」是錯的。可防守說法是 long reads 提供更遠、更密集的 joint linkage。
- [Bulk tumour subclonal reconstruction practical guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC7867630/) 說明 single 或 multiple bulk samples 均可在 purity、copy number、multiplicity 與模型假設下做 reconstruction；因此「single bulk 原理上完全不能 reconstruction」過度絕對，但這不表示結果唯一或無偏。
- [2024 DREAM benchmark](https://www.nature.com/articles/s41587-024-02250-y) 以 31 algorithms、51 truth-known simulated tumours、12,061 runs 比較 single-sample reconstruction；這證明問題可被量化 benchmark，也同時顯示 CN、depth、mappability 與算法選擇的重要性。
- [Mazor et al. phyloepigenetic study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4573399/) 在 multi-region／longitudinal glioma 中建立 methylation 與 mutation phylogenies；所以「methylation 不可能用於演化研究」是全領域錯誤。正確限制是：**InterSubMod 目前 single-bulk、read-level、ASM/LOH/CN-confounded 資料只支援 association-only**。
- Same-read pair 方法（PairClone、TreeClone）仍把 clone number、genotype、frequency、tree 當模型推論；這直接支持「同分子共現不是直接看到 cellular lineage」。

完整 source boundary：[`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/external_primary_sources.md`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/external_primary_sources.md)。

## 10. CCU lab-tutorial 現行版本重驗

> 本節比較前次完整稽核 snapshot `46b6f5b` 與目前 live/main `9eb1618`。目前版本共改 57 files、+8,105/−2,774 lines；因此不能沿用舊報告判定而不做 delta audit。

| Current-delta status | 數量 | 解讀 |
|---|---:|---|
| OPEN | 17 | 舊問題或新問題仍存在，尚未實質修正 |
| PARTIAL | 6 | guardrail／部分文字改善，但核心來源、數據或相鄰頁仍未同步 |
| REGRESSED | 6 | 新版新增了超出目前證據的新 claim |
| RESOLVED | 3 | 現行文字與算術／claim ceiling 已一致 |
| **尚未完全解決** | **29/32** | problem-focused tracking denominator；不是全站句子錯誤率 |

### 10.1 已確定的現況

M6/M7 的外部工具判定先查 [`longphase-s.md`](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md) 與 [`longphase-to.md`](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md)，再用 pinned local source 與 binary help 交叉驗證；完整 commit、binary hash 與 stdout 見 command receipt R20。

- Live `index.html` 與 remote main source SHA-256 一致。
- 25 個正式頁面 25/25 HTTP 200 且 hash match；`sr6.html` 是 repository source-only draft、live 404。
- M6 的 `HP1-1-1` 與「LongPhase-TO only supports ONT」沒有被新 commit 修掉。
- M7/glossary 雖已按工具區分 tag type，卻仍錯寫「整個 LongPhase-TO 使用 `HP:Z`」；目前 runtime/code contract 是 TO haplotag `HP:i:1/2`，也可能有 integer substates 11/21/33，而 LongPhase-S somatic 才使用字串 substate。
- SR2c 已明說 0.909 是 reproducibility、不是 accuracy；但可獨立重現的是 technical-pair **VAF Spearman=0.9093**，不是頁面另稱的 linked-region classification/composition agreement；SR2b 又使用「目前最強結果／承載真實訊號」的過強語氣。
- 64.89% 雖被重寫成 64.89%–88.26% 的 parsimony-attributable interval，不再稱 accuracy，但 64.89% 仍沒有 authority numerator/denominator 或可重播 receipt；因此這個 interval 仍不可完整確認。
- SR2c 已明確寫 local node 是 molecular state、不是 cell population、不能證明 global clone tree；SR2b 保留 global clone-tree estimand，但已有 caveat。
- M13 仍宣稱 coverage-mixed synthetic sample 中 tumour DNA fraction 與 cellular purity 恰好相等；若 `f=0.2, κ=3.2`，站內公式給 `p=0.135135`，只有 `κ=2` 與其餘模型條件成立時才有 `f=p`，故是 P0 直接矛盾。

### 10.2 新版新增的六個 regression

1. **把 HP families 各等同一條可識別的實體染色體 copy/homolog。** 在未校正 CN/LOH、沒有 parental origin 時，最多只能稱 phase/haplotag read families。
2. **斷言逐 family likelihood 相乘必然錯，且 misassignment 必然雙重高估。** 是否條件獨立與偏誤方向取決於完整生成模型，需 simulation/calibration，不能先下定論。
3. **宣稱 ξ 跨區差一個數量級、k≥2 時幾乎主導 error floor。** 沒有 ξ decomposition/calibration receipt，且現有式子未涵蓋 mapping、chimera、Hamming-1 與 partial-coverage states。
4. **把 background probability `δ·r` 直接稱 detection limit。** Detection limit 還取決於 n、α、power、multiple testing、state competition 與 detector 定義。
5. **宣稱五參數 optimizer 幾毫秒收斂。** 同頁承認 objective 未知是否 concave、需 multi-start，卻沒有 hardware/input/repetition/failure-rate benchmark。
6. **未驗證 LiLT 式子稱 production likelihood。** 同頁又明示 spec／尚未驗證；校準、simulation、truth-set 與 deployment gate 前只能稱 proposed likelihood。

### 10.3 現行仍需修正的最短清單

- [ ] M6 移除不存在的 `HP1-1-1`；同站 SR2c 已寫「恰好五種 extended values、沒有更長後綴」，目前是站內自相矛盾。
- [ ] M6 把「LongPhase-TO only ONT」改為釘定版本的 ONT/PacBio capability 說明。
- [ ] M7/glossary 移除「整個 LongPhase-TO 是 `HP:Z`」；按實際 binary/release 區分 integer `HP:i` 與 LongPhase-S somatic string tag，並附 runtime receipt。
- [ ] M13 移除 synthetic sample 中無條件 `f=p`；改為只有 diploid `κ=2`、等 DNA contribution 等模型條件成立時相等，並保留一般換算式。
- [ ] SR2b 將「0.909 承載真實訊號」收斂成「technical reproducibility」；若要保留 classification agreement，須補獨立 numerator/denominator、metric code 與 artifact，否則只報可重現的 VAF Spearman=0.9093。
- [ ] SR2c 的 64.89% 應移除或補 numerator、denominator、候選集合定義與 validator；88.26% 可保留但只限 ranked mathematical shapes。
- [ ] 首頁、M9、M12、M13、SR1 的既有 molecule→cell／aspiration→validated-result 問題，凡未在 current delta 明確修正者仍保留 P0 gate。
- [ ] 精確 benchmark 表、case coordinate/count、tool internals 繼續補 source/version/dataset/scope/numerator/denominator/metric definition。

逐項 delta machine table：[`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv)。

## 11. 必須補的資料：依「能改變 claim ceiling 的能力」排序

### P0 — 沒有這些資料，不得宣稱 cellular clone／lineage

1. **Allele-specific CN/LOH + purity + multiplicity + CCF table**  
   每個候選 locus/block 需有 caller、version、segment、confidence、normal contamination、purity/ploidy、multiplicity 與 CCF uncertainty；目前 coverage proxy 不等於 allele-specific CN correction。

2. **獨立 cellular truth**  
   優先次序：single-cell DNA／joint DNA-methylation、single-cell-derived colonies、spatial/multi-region、longitudinal樣本。至少能把 read-level state 對回 cell population，且不能只用同一 bulk reads 派生的 label。

3. **Known-truth simulations／mixtures**  
   模擬或實體 mixtures 必須預先固定 clone genotype、tree、frequency、CN、LOH、coverage、MM/ML error，並報 topology、frequency、abstention、calibration，而不是只報能否生成一棵 tree。

4. **Biological replication**  
   HCC1395 與 Dorado 是同一 cell line technical pair；需不同 biological samples 中重現同類 pattern，預先固定 feature、gate 與方向。

5. **Per-read tag truth 與 rescue benchmark**  
   若要宣稱 methylation rescue/correct tag，需有 truth labels、coverage strata、confusion matrix、abstention、calibration、false-rescue rate 與 held-out validation。

6. **Expression／functional assay**  
   在沒有 RNA、chromatin、perturbation 或功能性證據前，只能描述 methylation state，不得說哪個基因已被 clone 關閉。

### P0 — 公開重現性資料

7. **Portable tiny fixture + manifest**：BAM/CRAM、FASTA、VCF、index、expected outputs、license、exact SHA-256、download step。
8. **Claim registry**：claim ID、文字、source artifact、validator、numerator、denominator、unit、grain、commit、date、status、owner。
9. **Version manifest**：GitHub main、Wiki、Pages deploy、binary、schema、data authority 必須一次性對齊，禁止 current 页面混用不同 commit。
10. **Doc tests in CI**：HTTP/route、code block smoke test、CLI help snapshot、schema header、test count、file count、numeric formula、source existence。

### P1 — 目前四個不可驗證 claim 的直接補件

11. **Tagged-BAM byte registry**：7 個 exact paths、bytes、compression、hash、generation command；才可重算 1.67 TiB 與縮減倍率。
12. **BAM field-level byte census**：才能支持 `SEQ/QUAL >99%`；另需 downstream utility definition 才能支持 `0% use`。
13. **Performance benchmark receipt**：hardware、compiler、commit、input distribution、warm-up、repetitions、median/IQR/tails、thread scaling；才可談 32-core linear 與 <300 ms。
14. **Dated systematic novelty review**：search databases、query、date、inclusion/exclusion、head-to-head baselines；否則移除 first/unique。

### P1/P2 — 降低算法與格式漂移

15. MM/ML grammar fixture matrix，覆蓋 canonical、multi-mod、strand、skip、malformed、boundary cases。
16. Whole-file schema version 與 backward-compatibility test，不只 component-level schema fields。
17. Multi-axis semantics tests，明確說明「repeatable option」與「實際只採第一有效 axis」是否設計一致。
18. `tree.nwk`、`methylation.csv`、deprecated `subclone_structure` 的 machine-readable schema 與 semantic warning。

## 12. 31 個 InterSubMod 公開 artifact 的頁面級裁決

這是「該 artifact 最嚴重 claim」的分類；不是整頁錯誤率。

| 類別 | Artifact IDs | 數量 | 判定 |
|---|---|---:|---|
| CONTRADICTED | A002–A006（除 A001）、W001–W003、W005–W007、P000–P004、P008、P010–P012、P014–P016 | 23 | 至少一項可直接反駁，需修 |
| NEEDS_CORRECTION | A001、W008、P005–P007、P009 | 6 | 沒有必要整頁撤下，但至少一項需降級／重寫 |
| CONFIRMED_WITH_LIMITS | W004、P013 | 2 | 核心可留；須釘 branch/version 與 production ceiling |
| CONFIRMED | — | 0 | 沒有任何現存 artifact 可在完全不加 caveat 下整頁背書 |
| Missing endpoint | M001=`main/README.zh-TW.md` | 1（不納入 31） | HTTP 404，屬發布缺口 |

逐 artifact 的 URL、commit、hash、HTTP 與角色見 `source_scope.tsv`；逐 claim occurrence 見 `claim_inventory.tsv`。

## 13. 建議發布順序與驗收 Gate

1. **先修主入口**：GitHub About、default README EN/ZH、Wiki Home、Pages index。  
   → 驗證：四面 commit/version manifest 一致；所有 cellular claim 經 claim-ceiling linter。
2. **再修可機械反駁項**：denominator、199-column schema、NHD default、test/file counts、quickstart、storage arithmetic。  
   → 驗證：CI 從 artifact 自動產數，禁止手抄。
3. **再修 chr2 教學因果語氣**：LOH、α/β、methylation、rescue、gene silencing。  
   → 驗證：每一節都有 observation／inference／alternative／missing truth 四欄。
4. **補公開 fixture 與 receipt registry**。  
   → 驗證：全新 clone 在無 `/big7`、`/big8` 的環境完成 build＋smoke，exit 0。
5. **取得 P0 科學資料後再升格**。  
   → 驗證：CN/LOH/purity-aware、independent cellular truth、跨 biological sample、pre-registered held-out benchmark 全部通過；未通過前維持 molecular-state candidate 用語。

## 14. 限制與反證條件

- 本次「完整」指全部指定的 public surfaces：GitHub About、default first-hop docs、公開 feature README、獨立 Wiki 全頁、live Pages 全 17 頁，以及目前 CCU lab-tutorial 25 正式 routes＋1 source-only draft；不等於逐句稽核 repository 內每一份歷史/internal Markdown。
- Fresh InterSubMod build 使用本地 build-input-equivalent tracked source；不是 clean remote feature clone。證據足以驗證 current code path，但交付前仍建議在 clean CI runner 重播。
- Wilson interval 可描述 `63,506/71,955` 或 `3/811` 的二項不確定度，但 units 非隨機抽樣的生物樣本，因此不把該 CI 解讀成 population generalization。
- External papers 只用來反駁過度絕對的領域敘述，不把 paper-scope 性能當成本專案 accuracy。
- 本次沒有修改 GitHub remote、Wiki 或 live site；所有修正仍是 checklist，需另開受控文件更新工作。

## 15. 可重現執行收據

完整輸入路徑、命令、輸出路徑、exit code 與 stdout/stderr 片段集中於：

- [`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md)
- [`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/verification_results.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/verification_results.tsv)
- [`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/pre-decision-audit.md`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/pre-decision-audit.md)

最終證據完整性標準：所有精確數字能回到 numerator/denominator/unit/grain；所有 capability 能回到 repo/branch/commit/symbol/runtime；所有 biological claim 不高於其最弱必要證據層。
