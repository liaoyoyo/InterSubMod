<!--
建立時間: 2026-07-10 23:42 Asia/Taipei
目標: 稽核 layered reconstruction 全流程的測試、驗證、數值 provenance 與 adversarial risk
處理範圍: V1-V7、full_v4v5_verification、seeded stress、funnel、7 datasets、CCF/CN、既有 HTML/MD claim、tests/
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md
任務類型: B — Comprehensive validation
服務目標: G4 reproducibility；G5 externally verifiable contribution
狀態: audit_complete_at_moving_snapshot
快照分支: research/subclonal-reconstruction-202606
快照 HEAD: 4fb9e742482b63a660de19a1f1bd07d49d713111
限制: 未執行 full BAM / 7-dataset 長計算；v2 run 於稽核末段仍在執行
-->

# Verification、數值來源與 adversarial risk 稽核

用 SCQA + claim-to-evidence matrix：先給發布判斷，再把每個數字、驗證命題與不能下的結論逐項對回程式、輸出與反例。

> **TL;DR — 舊 HCC1395 HTML 不可當作已完成、可發布的完整驗證報告；funnel 守恆與小型測試可重現，但 recurrence 分類直接矛盾、V5 不是候選樹集合完整性 oracle、4,000 stress 把 capped 未驗案例算作 pass，且主要數據來源未納入 Git。**（影響：高；信心：高）`[O-L1]`

## 0. 稽核判斷與證據語彙

### 0.1 最終判斷

- **[F-L1] 程式正確性（bounded code paths）**：solver 的 8 個 golden cases、v2 的 5 個 synthetic tests、C++ HP-tag 的 12 個 tests 在本快照都通過；這只能證明被覆蓋的小型案例。
- **[F-L1] 內部一致性**：HCC1395 六層 funnel `113,997/113,997` 與 legacy tree class `13,785/13,785` 守恆；7 個 legacy dataset 的 funnel `check_ok=true`。
- **[O-L1] 驗證語意缺陷**：V5 只重算 analytical count，未比對輸出候選樹集合；固定反例在 `n_trees=1,119,744`、只存 32 棵時仍回 `V5=True, overall=True`。
- **[O-L1] 報告數值矛盾**：commit `d2fd4be` 與 `11452e5` 都寫 HCC1395 recurrence `20 CN-gain / 7 LOH / 0 candidate`，其標示來源 `ccf_and_cn_multisample.json` 實為 `10 / 13 / 4`。
- **[O-L1] stress 宣稱過度**：legacy stress code 對 capped case 不執行 V4/V5，卻直接 `stress_pass += 1`；固定 seed 前 100 案有 2 案屬此類。
- **[F-L1] biological validation**：repo 內沒有 single-cell / multi-region truth；ClairS PASS、read co-occurrence、CN annotation 與 read-AF self-consistency 都不能單獨確認 biological subclone tree。
- **[I-L2] 發布 gate**：舊 HTML 應標 `superseded / not validation evidence`；必須等 v2 全 7-dataset run 完成、重新產生 primary-lineage/all-candidate 數字、修正 recurrence、分離 checked/skipped，才能寫工程驗證通過。

### 0.2 Claim tag 定義

| Tag | 定義 |
|---|---|
| `[F-L1]` | 直接由版本化程式或機器輸出讀得的事實；L1 為最直接證據。 |
| `[O-L1]` | 本稽核實際執行命令得到的觀察。 |
| `[I-L2]` | 由兩個以上 L1 事實推出的限制或判斷。 |
| `[U-L5]` | 目前資料或執行狀態不能回答。 |
| L2-L4 | 可重現內部分析、既有專案報告、間接／外部佐證；不可冒充 L1。 |

### 0.3 三種「驗證」必須分開

| 層 | 可回答 | 不能回答 | 本輪狀態 |
|---|---|---|---|
| 程式正確性 | 同一輸入是否符合程式內定義的 V1-V7；schema 與 regression 是否過關 | 模型定義是否符合真實腫瘤演化 | **局部通過** `[O-L1]` |
| 內部數值一致性 | 加總、分母、跨檔欄位、checksum 是否一致 | caller universe 是否是真 somatic truth | **部分通過且有矛盾** `[O-L1]` |
| 生物學驗證 | inferred tree 是否被 orthogonal truth 支持 | 不能由 self-consistency、CN annotation 或同一 bulk data 自我證明 | **未完成** `[F-L1]` |

## 1. 快照、範圍與 moving-target 風險

- **[O-L1] 稽核快照**：branch=`research/subclonal-reconstruction-202606`；HEAD=`4fb9e742482b63a660de19a1f1bd07d49d713111`。
- **[O-L1] working tree 非乾淨**：核心 `tree_enumeration_solver.py`、`layered_tree_reconstruction.py`、`sm_multilocus_combinations.py`、`sm_linkage_genomewide.py`、`ccf_and_cn_multisample.py`、runner 等均有未 commit 變更；`test_layered_reconstruction_v2.py`、`verify_layered_v2.py`、input validation 等為 untracked。
- **[I-L1] 因此本文件同時稽核兩個不同狀態**：已 commit 的 legacy HTML claim，以及當下 working-tree 的 v2 remediation；兩者不可混寫成同一個「已驗證版本」。
- **[O-L1] 稽核末段已有 v2 run 建立**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/`，但 `run_status.tsv` 當時只有 HCC1395 與 HCC1395_DORADO 的 `mlhp START`，沒有 `verification_summary.json`。
- **[U-L5] 本文件未確認 v2 full run 最終是否成功**；最終主報告必須重新讀取該 run 的 completion、checksums、verification summary 與失敗 stage。

| 狀態層 | 內容 | 可重現性判斷 |
|---|---|---|
| HEAD `4fb9e74` | `[F-L1]` legacy committed code/HTML；含被本 audit 否定或降級的舊 claims | 可 checkout，但 ignored JSON 不隨 commit，不能重現完整報告 evidence。 |
| current working tree | `[O-L1]` v2 semantics、runner、verifier、tests；多檔 modified/untracked | 不能單靠 HEAD checkout；必須保留 patch/content。 |
| active v2 run | `[O-L1]` 於啟動時寫 `code.sha256` 與 manifest SHA；C13 當下全數仍 `OK` | hash 可偵測漂移，但 run root **沒有複製 source code**；只有 hash 不能重建未 commit 內容。 |

- **[O-L1] v2 source manifest 本身仍被忽略**：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json` 命中 `.gitignore:10:data/`。
- **[O-L1] prevalidated input lock 仍是 untracked**：`InterSubMod/research/20260710_layered_reconstruction_v2/input_manifest_lock.json` 顯示 `??`。
- **[I-L1] active run 的 `input_manifest.json` copy 與 hashes 能保存「這次用了什麼」，但不能替代 repo-level manifest review，也不能保證未 commit code 可被未來重建**；final evidence package 應複製 source/patch 或綁定 commit。

### 1.1 v2 run 已釘住的設定

| 設定 | 值 | 稽核說明 |
|---|---:|---|
| scope | chr1-22 | `[F-L1]` 不含 chrX/Y；所有「全基因組」字眼都應改為 autosomal scope。 |
| dataset_count | 7 | `[F-L1]` 是 7 個 dataset，不是 7 個獨立生物樣本。 |
| biological_sample_count | 6 | `[F-L1]` HCC1395 與 HCC1395_DORADO 共用 `biological_id=HCC1395`。 |
| VERIFY_EVERY | 1 | `[F-L1]` v2 預設每個 non-capped unit 跑完整 V1-V7。 |
| ANALYSIS_TREE_CAP | 0 | `[F-L1]` v2 working tree 要求分析使用全部候選樹。 |
| DISPLAY_TREE_CAP | 32 | `[F-L1]` UI 仍只顯示前 32，但 schema 應明示 display-only。 |
| MINREAD | 3 | `[F-L1]` full population / partial pattern 至少 3 reads。 |
| MAX_SNV | 8 | `[F-L1]` 大 group 只取 span 最小的連續 8 個 SNV。 |
| TIER_R | 50,000 bp | `[F-L1]` 以相鄰 SNV gap 切 component，**不是 total span ≤50 kb**。 |
| MAPQ_MIN | 20 | `[F-L1]` mapping quality 下限。 |
| BASEQ_MIN | 0 | `[F-L1]` 不額外排除低 base quality；此設定需要 sensitivity。 |

## 2. Claim → source → command → actual output 驗證矩陣

| ID | 被稽核 claim | 來源與證據 | 實際輸出 | 狀態／可下結論 |
|---|---|---|---|---|
| V01 | `[F-L1]` HCC1395 六層 funnel 113,997 守恆 | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/funnel_census_HCC1395.json`；C02 | `funnel=113997/113997 equal=true` | **確認內部一致**；不確認 113,997 為 biological truth。 |
| V02 | `[F-L1]` legacy tree class 13,785 守恆 | 同上；C02 | `tree_sum=13785/13785` | **確認內部一致**；此分母含 HP3 與 root-only legacy units，不等於 v2 primary lineage。 |
| V03 | `[O-L1]` recurrence 27 全是 CN artifact/LOH、真 candidate=0 | HTML commits `d2fd4be`、`11452e5`；`ccf_and_cn_multisample.json`；C01/C02 | HTML=`20/7/0`；JSON=`10/13/4` | **直接矛盾，發布 blocker**。不能保留「真收斂 0」。 |
| V04 | `[F-L1]` 18,103 non-capped units V4/V5 0 fail | `full_v4v5_verification_HCC1395.json`；C02 | total=18,931；skipped=828；run/pass=18,103；fail=0 | **數字可讀，語意只部分確認**；V5 缺候選集合比對，且分母含 5,146 unphased + 1,692 HP3。 |
| V05 | `[O-L1]` 4,000 stress 全部經 V4/V5、0 mismatch | `full_v4v5_verification.py:91-98`；C07 | 前100：checked=98、capped_not_tested=2、fail=0 | **宣稱不成立**；capped 被算 pass。完整 4,000 中實際 skipped 數未保存，現有 JSON 無法回推。 |
| V06 | `[O-L1]` solver golden cases | `tree_enumeration_solver.py`；C04 | 8 cases，`GOLDEN 總結: ALL PASS` | **確認 8 個手工小例**；沒有 negative corruption、budget boundary、large-tree-set oracle。 |
| V07 | `[O-L1]` V5 證明「沒有漏任何候選樹」 | `tree_enumeration_solver.py:349-398`；C08 | `n_trees=1119744, stored=32, trees_complete=false, V5=true, overall=true` | **反例否定**；現行 V5 只驗 analytical count 一致，不驗輸出集合完整。 |
| V08 | `[F-L1]` HCC1395 受 32-tree cap 影響規模 | legacy layered JSON；C09 | units=18,931；over32=371；max=125；sum candidates=66,837 | **確認**；舊 shape/read-AF 只看前32，371 units 可能偏差；v2 cap=0 的實際資源成本仍須 full run。 |
| V09 | `[I-L2]` topology determined=51.1%、真正多拓撲僅2.0% | HTML；untracked `shape_vs_exact_ambiguity.py`；V08 | 來源 JSON 未被保存／追蹤；script 用 solver default cap=32 | **不可確認**；需在 v2 primary units + all candidates 上重算。不能稱「資訊論極限」。 |
| V10 | `[I-L2]` reliable unique=5,384+365=5,749 | quantitative HTML；legacy CCF JSON | 5,384 與 365 來自不同 classification/CN subpopulation；舊 lineage 含 HP3 | **不可確認**；需證明集合互斥、population 相同及 role semantics 相同。 |
| V11 | `[F-L1]` capped 676 全部可解 | commit `d2fd4be` header/§9；later `capped_cluster_conflict.json` | later result=666 pairwise-compatible +10 conflict | **原 claim 被後續結果否定**。而 compatible 只代表程式未見 observed pairwise four-gamete conflict，未產 witness tree。 |
| V12 | `[F-L1]` 7-dataset legacy funnel 守恆 | `funnel_census_7samples.json`；C10 | 7/7 `check_ok=true` | **確認內部守恆**；沒有獨立重讀 BAM/VCF 的第二實作。 |
| V13 | `[F-L1]` 「7 samples」為 7 個獨立 replication | v2 input lock；C11 | datasets=7；biological_samples=6 | **否定**；應寫 7 datasets / 6 biological samples，DORADO 為同細胞株 pipeline replicate。 |
| V14 | `[O-L1]` CN availability 一致 | legacy funnel vs v2 lock | legacy `has_cn=true` 只有 HCC1395；v2 lock available=5/7 | **版本漂移**；報告只能用同一 run manifest 的狀態，禁止跨 artifact 拼接。 |
| V15 | `[F-L1]` input lock 完整且受版本控制 | `InterSubMod/research/20260710_layered_reconstruction_v2/input_manifest_lock.json`；C11/C13 | schema2.0、7/6、all_pass=true、CN 5/7；Git=`??` | **部分確認**；VCF/CN/index 有 SHA，BAM 只有 size/mtime/header/ref dictionary，runner 未比較 prevalidation lock，manifest/lock 皆未被 Git 追蹤。 |
| V16 | `[F-L1]` HCC1395 strict-neutral trustworthy reach=6.1% | `ccf_and_cn_multisample.json` | 597 neutral ambiguous；365 tie-broken；reach=6.1% | **確認為內部 heuristic 指標**；不是 CCF truth accuracy，也不是 independent validation。 |
| V17 | `[O-L1]` upstream HP-tag parser regression | C06 | 12/12 pass | **確認該 C++ 子路徑**；不涵蓋 Python family-role/denominator。 |
| V18 | `[O-L1]` v2 schema remediation tests | C05 | 5/5 pass；另有 2 個 ResourceWarning | **局部確認**；test 未納入 ctest/CI，沒有 full outputs 與 sensitivity。 |
| V19 | `[U-L5]` v2 7-dataset final outputs 完成且全數通過 | active run root；C12 | 稽核時只有兩 dataset `mlhp START` | **尚不可確認**；主報告交付前必須重查。 |
| V20 | `[O-L1]` legacy HTML「provenance-verified、全數機械重算」可由 commit 重現 | C01/C03 | 三個核心 JSON 被 `.gitignore:10:data/` 排除；HTML 沒有 source SHA/run ID | **不可重現**；commit 只保存呈現，不保存精確 evidence package。 |
| V21 | `[F-L1]` regional mutation-state tree 已 biological confirmed | repo input audit / implementation notes | single-cell/multi-region truth unavailable | **未驗證**；claim ceiling 僅為 model-based regional inference。 |
| V22 | `[F-L1]` 基本門檻的 robustness 已確認 | runner params + code search | 只有單組預設；無 MAPQ/BaseQ/MINREAD/MAX_SNV/TIER_R sensitivity matrix | **未驗證**；不得寫「對設定穩健」。 |

## 3. 兩份既有 HTML 的 commit-level source consistency

### 3.1 `d2fd4be` narrative HTML

- **[O-L1] commit 內容**：`d2fd4beef86f90c6a8efc04880654e442ffe0586` 的 HTML SHA-256（以 `git show` stdout 計）為 `f26009fcdec68f3e25055bc505b30cc9770091872b87cb68cdda70da63c9def3`。
- **[O-L1] recurrence 矛盾**：§7/§8 寫 27=`20 CN-gain +7 LOH +0 real`，但同頁列出的 `ccf_and_cn_multisample.json` 是 27=`10 CN-gain +13 LOH +4 candidate`。
- **[O-L1] capped 前後矛盾**：原 commit header 與 §9 寫「676/676 可解」；後續 `0838494` 改成 666 compatible +10 conflict，證明原句不是穩定結論。
- **[I-L2]「51.1% CN-robust、其餘是單一 bulk 資訊論極限」過度推論**：演算法不讀 CN 只能證明計算路徑不依賴 CN；不能證明真實 topology 不受 CN/mapping/homoplasy confounding。
- **[F-L1] current narrative 已經過 `0838494`、`6b0f917` 修改**；audit 必須區分原 commit 與 HEAD，不可只說「d2fd4be 的報告」。
- **[F-L1] 原頁確實有 partial flag**：header 寫 `scope=chr1–22`、`HCC1395 單樣本 ⭐3`；問題不是完全沒標 partial，而是後續用語仍把 model/internal evidence 擴張成「CN-robust」「資訊論極限」，且沒有可 checkout 的 source package。

### 3.2 `11452e5` quantitative HTML

- **[O-L1] commit 內容**：`11452e5161aaf92f042b2cad8ae9805b6d313113` 的原 HTML SHA-256 為 `0ab08bfb13a69c48ce9d2b683da095678ffb91080f395667a786819b15e5061d`；HEAD 後由 `4fb9e74` 加入 funnel SVG。
- **[O-L1] denominator 明顯錯置**：表格段落標「分母=13,785 lineage units」，卻在同段放 `18,103` V4/V5 units；實際 18,931 全 units=`6,132 HP1 +5,961 HP2 +1,692 HP3 +5,146 unphased`。
- **[O-L1] source path 漂移**：HTML 列 `capped_cluster_conflict.py`，實際檔名為 `capped_cluster_conflict_check.py`；列 `full_v4v5_verification.json`，實際為 `full_v4v5_verification_HCC1395.json`。
- **[I-L2] 5,749、51.1%、80.3%、2.0% 均需撤回後重算**：它們混合 legacy role、stored-first-32 shape、read-AF heuristic 與不同 CN subset，現有 commit 無一份 machine-readable join table 能逐 unit 重現。

### 3.3 Evidence package 缺口

| 檔案 | current SHA-256 | Git 狀態 | 影響 |
|---|---|---|---|
| `InterSubMod/.../data/funnel_census_HCC1395.json` | `d81f43c2c612206b663a05d95f893b0afc15e5a022bbb000ef1f56c1df1508f6` | ignored | `[F-L1]` commit 無法釘住數字來源。 |
| `InterSubMod/.../data/ccf_and_cn_multisample.json` | `335d37601cb318ec9562a3c755cddad114b05f62951f1c79bcd3426a0a91f813` | ignored | `[F-L1]` recurrence conflict 無法由 HTML commit 自行揭露。 |
| `InterSubMod/.../data/full_v4v5_verification_HCC1395.json` | `892d43f82ef371492879f58adcc21bac007085065ed4aabdb5f261e426de5fce` | ignored | `[F-L1]` stress skipped count 本來就未記錄。 |

- **[I-L1] 最低可接受修正**：final run root 必須包含 input manifest、input hashes、code hashes、params、per-stage status、per-output hashes、verification JSON/TSV 與 exact report source table；HTML 內嵌 run ID + hashes，而不是只寫「本 session 機械重算」。

## 4. V1-V7 的實際驗證邊界

| Check | 程式真正檢查 | 優點 | 缺口／不能宣稱 |
|---|---|---|---|
| V1 | `[F-L1]` stored trees 是 root-connected arborescence、無 cycle | 可抓 edge 結構錯誤 | 只巡覽 `result["trees"]`；若 output set 被截斷，不看未存候選。 |
| V2 | `[F-L1]` full observation 在 node；partial subcube 與 node set 有交集 | 明確 coverage invariant | coverage 定義與 solver 共用 parse/subcube semantics，非獨立資料 oracle。 |
| V3 | `[F-L1]` full observations 是 nodes subset | 基本一致性 | 與 V2 大幅重疊；不能驗證 reads 是否被正確生成為 genotype。 |
| V4 | `[F-L1]` 以同一 pool/_covers/_is_closed 重搜 `e_min-1` | 可抓部分 early-stop/minimality bug | 稱「independent」不精確：仍共用核心 helper 與 candidate construction；capped skip。 |
| V5 | `[F-L1]` 以同一定義重算 analytical `n_trees` 並比對數字 | 可抓 count formula／report mismatch | **不生成或比較完整 tree set**；C08 已有 `trees_complete=false` 仍 pass 的反例。 |
| V6 | `[F-L1]` classification flags 與 recurrence/refree flags 自洽 | 可抓旗標矛盾 | flags 由同一 solver 結果產生，屬 self-consistency。 |
| V7 | `[F-L1]` determined 必須 n_trees=1、recurrence-free、not capped | 防止最直接 overclaim | n_trees 的 correctness 仍仰賴同一枚舉/分析式。 |

- **[I-L1] 正確用語**：現行 V1-V7 是「implementation invariant suite」；V4/V5 可稱 cross-check，但在沒有獨立實作、tree-set equality 與 negative oracle 前，不能稱 formal independent oracle。
- **[O-L1] `verify_all(light=True)`** 會將 V4/V5 設為 `None`，再對其他非 `None` checks 求 `overall`；legacy consumer 因此可能把 partial check 寫成 overall pass。
- **[F-L1] v2 working tree 已新增 `verification_status={full_pass,partial_pass,fail,not_applicable_capped}` 並只讓 `full_pass` 設 `verify_pass=true`；這是合理修正，但尚未由 completed 7-dataset artifact 驗證。

## 5. Test coverage matrix

| 元件／風險 | 現有測試 | 本輪結果 | 是否進 ctest/CI | 缺口 |
|---|---|---|---|---|
| C++ HP tag demotion / pass-through | `ReadParserHPTagTest.*` 12 tests | `[O-L1]` 12/12 pass | 是，ctest | 不測 Python HP1/2/3 role 與 denominator。 |
| solver basic topology | 8 golden cases | `[O-L1]` all pass | 否，直接執行 script | 案例少；無 corrupted expected output、large set、budget edge。 |
| V1 arborescence | golden + legacy full run | `[F-L1]` stored trees pass | 否 | 未測故意 cycle/disconnected tree 注入。 |
| V2/V3 coverage | golden + legacy full run | `[F-L1]` stored trees pass | 否 | parse/subcube 與 solver 同源；缺第二實作。 |
| V4 minimality | legacy 18,103 + golden | `[F-L1]` reported 0 fail | 否 | 共用 candidate pool/helpers；capped 全跳過。 |
| V5 completeness | legacy 18,103 + golden | `[O-L1]` counterexample shows false sense of completeness | 否 | 缺 explicit canonical tree-set equality。 |
| stress/fuzz | 5 seeds ×800 legacy | `[O-L1]` capped 被算 pass | 否 | 未存 skipped count、generator version、case corpus。 |
| root-only / HP3 / missing CN | v2 synthetic tests | `[O-L1]` 5/5 pass | 否；test untracked | 只有單一 tiny fixture。 |
| analysis cap vs display cap | v2 synthetic test | `[O-L1]` pass | 否 | 未測真實 371 units 與 output size。 |
| funnel conservation | producer self-check + verifier | `[F-L1]` HCC/legacy 7 dataset pass | 否 | 同一輸出內自我加總，缺 upstream VCF position recount in verifier。 |
| input provenance | v2 input validator/lock | `[F-L1]` 7 dataset all_pass | 否 | BAM 無 full hash；CN coverage semantics 未驗；runner 不 enforce lock。 |
| CCF/read-AF ranking | legacy observation scripts | `[F-L1]` 無 unit test | 否 | 無 CN/purity truth、no tie/symmetry/threshold boundary tests。 |
| recurrence CN split | legacy JSON | `[O-L1]` HTML 與 source 矛盾 | 否 | 缺 per-unit join table + regression expected counts。 |
| capped compatibility | pairwise four-gamete script | `[F-L1]` 666/10 | 否 | partial `X` 只忽略；無 completion/witness tree；無 exact subset cross-check。 |
| 7-dataset extension | legacy census、active v2 run | `[U-L5]` v2 final 尚未完成 | 否 | 7 datasets 只有 6 biological samples；無 external truth。 |
| parameter sensitivity | 無 | `[F-L1]` 未測 | 否 | MAPQ20/BQ0/MINREAD3/MAX8/gap50k/posterior thresholds 皆無 robustness。 |
| resource/file lifecycle | v2 tests | `[O-L1]` 2 ResourceWarning | 否 | `open()` 未 context manager；長跑可能累積 descriptor risk。 |

- **[O-L1] CMake/ctest 目前列 234 tests**，但 `rg` 未找到 layered Python solver/funnel/CCF scripts 被 CMake、`tests/` 或 `.github` 引用。
- **[O-L1] `tests/` 的 working-tree diff 是大量 `100644→100755` mode change，沒有對本 Python 流程新增內容**；應避免把 executable-bit drift 誤報成 coverage 增加。

## 6. Risk register

| ID | Risk | Severity | Likelihood | 影響 | 證據 | 必要修正 |
|---|---|---:|---:|---|---|---|
| R01 | `[O-L1]` recurrence 20/7/0 與 source 10/13/4 衝突 | **Critical** | Certain | 改變「沒有真 candidate」核心結論 | V03 | 由 per-unit join 重產；HTML/summary regression 檢查三類加總。 |
| R02 | `[O-L1]` V5 將 count equality 說成 tree-set completeness | **Critical** | Certain | 錯把截斷/漏樹視為完整驗證 | V07 | 更名 `V5_count_consistency`；新增 canonical tree-set oracle，或在不可列舉時明示 proof obligation。 |
| R03 | `[O-L1]` stress capped 未驗卻記 pass | High | Certain | 4,000/4,000 數字被高估 | V05 | 分開 `generated/checked/skipped/failed`；skipped 不得進 pass 分母。 |
| R04 | `[F-L1]` first-32 storage 影響 371 units 的 shape/read-AF | High | Certain | 51.1%、2.0%、ranking 可能偏差 | V08/V09 | v2 全 candidate analysis；完成後逐 unit 比對 digest/count。 |
| R05 | `[F-L1]` old HTML evidence JSON ignored、部分 scripts untracked | High | Certain | commit 不可重現且 source 可無痕漂移 | V20 | 將 compact evidence package 版本化，或 immutable run root + manifest SHA + artifact index。 |
| R06 | `[I-L1]` four-gamete compatible 被寫成「乾淨樹已存在／全可解」 | High | High | 對 partial data 與 capped topology overclaim | V11 | 僅稱 observed pairwise no-conflict；建立 completion/witness tree，並在 exact tractable subset 對照。 |
| R07 | `[F-L1]` read AF 被稱 CCF 或 biological tie-break | High | High | purity/CN confounding 造成錯誤 clone ordering | V16 | 全部改名 `read_af_ordering`; purity/CN corrected 前不得稱 CCF。 |
| R08 | `[O-L1]` 分母混合 13,785 lineage、18,931 all units、primary HP1/2 | High | Certain | 百分比無法解讀、可超出段落分母 | V04/HTML | 每 metric schema 強制 `population_id`, `n_total`, role filters。 |
| R09 | `[F-L1]` 7 datasets 被敘述成 7 independent samples | High | Certain | pseudo-replication、跨樣本信心膨脹 | V13 | 報 7 datasets/6 biological samples；HCC兩流程作 technical/pipeline replicate。 |
| R10 | `[O-L1]` legacy CN availability 1/7 vs v2 5/7 | High | Certain | 同頁拼接不同 policy 產生假一致 | V14 | 僅使用同一 run manifest；所有表帶 schema/run ID。 |
| R11 | `[F-L1]` 無 orthogonal biological truth | High | Certain | 不能升級 clone confirmation/⭐4 | V21 | 保持 regional mutation-state inference claim；另立 single-cell/multi-region validation。 |
| R12 | `[F-L1]` threshold sensitivity 缺失 | High | Certain | 結論可能由基本設定驅動 | V22 | 至少做 one-factor + joint stress：MAPQ、BQ、MINREAD、MAX_SNV、gap、posterior threshold。 |
| R13 | `[O-L1]` 核心修正未 commit、run 來自 moving tree，run root 又未複製 source | High | High | reviewer 無法 checkout 或由 hash 還原同版本 | §1/C13 | 完成後 commit，或保存 source tar/patch bundle + dirty status；執行期間禁止並行修改。 |
| R14 | `[F-L1]` input manifest/lock ignored/untracked，且 input lock 不由 runner enforce | High | High | prevalidated files可能在 run 前漂移；reviewer看不到 manifest evolution | V15/C13 | runner 比較 manifest/lock SHA、fail closed；把 compact manifest/lock納入版本控制。 |
| R15 | `[F-L1]` CN BED「未列即 neutral」語意未被機器驗證 | High | Medium | missing interval 被錯判 diploid neutral | input validator code | manifest 加 coverage contract；檢查 contig coverage/producer semantics。 |
| R16 | `[F-L1]` Python regression tests 不在 CI 且 untracked | Medium | Certain | 後續修改可默默復發 | test matrix | 加 pytest/unittest discovery target，CI fail gate。 |
| R17 | `[O-L1]` ResourceWarning: unclosed read/write files | Medium | Medium | 長計算 descriptor leak、flush risk | C05 | `with open(...)`；warnings-as-errors test。 |
| R18 | `[O-L1]` HTML source filenames 已失效 | Medium | Certain | 使用者無法重跑或找到來源 | V10 | 自動產生 source manifest，不手寫檔名。 |
| R19 | `[F-L1]` BASEQ_MIN=0 且無 sensitivity | Medium | High | 低品質 base 可能改變 rare pattern | params | 顯示設定理由；測 BQ10/20/30。 |
| R20 | `[I-L2]`「algorithm does not consume CN」被當成「CN-robust」 | High | High | 機制獨立被錯寫成 biological robustness | narrative HTML | 改寫為 algorithmically CN-unadjusted；用 stratified CN sensitivity 才可談 robustness。 |
| R21 | `[F-L1]` legacy `clean={neutral,LOH,loss}` artifact 仍可被引用 | High | Medium | 把 LOH/loss 的 read AF 當 CCF | old CCF scripts/JSON | deprecate old output；strict neutral only，並在 artifact index 標 superseded。 |
| R22 | `[I-L2]` ClairS PASS universe 被誤當完整 somatic truth | High | Medium | funnel coverage/accuracy 被高估 | funnel metadata | 報 operational backbone；另做 truth subset/caller sensitivity。 |

## 7. 實際執行命令、輸入、輸出與 exit code

> [F-L1] 以下命令均在 `/big7_disk/liaoyoyo2001/InterSubMod` 執行；除 ctest 的既有 log 與 v2 test 的暫存目錄外，本稽核命令不寫研究輸出。正式稽核文件輸出為 `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/03_verification_risk_audit.md`。

### C01 — 精確讀取兩個 HTML commits

**輸入**：Git objects `d2fd4be`、`11452e5`。

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
git show d2fd4be:docs/methodology/20260710_reconstruction_full_pipeline_narrative_HCC1395.standalone.html \
  | rg -n 'recurrence|20 / 27|7 / 27|51\.1|66\.8|provenance|V4|V5|capped'
git show 11452e5:docs/methodology/20260710_full_pipeline_quantitative_table_HCC1395.standalone.html \
  | rg -n 'recurrence|51\.1|66\.8|18,103|13,785|5,749|V4|V5|capped'
```

**輸出**：stdout；**exit code**：兩命令皆 0。

```text
d2fd4be: recurrence ... artifact(CN-gain) 20 / LOH 7 / 真收斂 0
11452e5: 分母 = 13,785 ... V4/V5 驗證 18,103 ... 5,749
```

### C02 — funnel、tree class、recurrence、V4/V5 機械重算

**輸入**：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/{funnel_census_HCC1395,ccf_and_cn_multisample,full_v4v5_verification_HCC1395}.json`。

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
jq -r '"funnel=" + ((.funnel_sSNV.L2_out_of_scope_chrXY + .funnel_sSNV.L3_positional_singleton + .funnel_sSNV.L4_cap_excluded_densest8 + .funnel_sSNV.L5_read_unsupported + .funnel_sSNV.L6_retained_sSNV)|tostring) + "/" + (.funnel_sSNV.L1_all_pass_universe|tostring), "tree_sum=" + ((.tree_level.determined_mutation_bearing + .tree_level.L2a_root_only_reference_family + .tree_level.ambiguous + .tree_level.capped + .tree_level.recurrence)|tostring) + "/" + (.tree_level.n_lineage_units|tostring)' docs/methodology/_assets/20260627_subclone_4axis_teaching/data/funnel_census_HCC1395.json
jq -r '.samples.HCC1395.recurrence | "recurrence=\(.n_recurrence) candidate=\(.cn_verdict["candidate(real?)"]) loh=\(.cn_verdict.LOH) cn_gain=\(.cn_verdict["artifact(CN-gain)"])"' docs/methodology/_assets/20260627_subclone_4axis_teaching/data/ccf_and_cn_multisample.json
jq -r '.PART_A_full_v4v5 | "v4v5_total=\(.n_units_total) capped_skipped=\(.n_capped_skipped) run=\(.n_v4v5_run) pass=\(.n_v4v5_pass) fail=\(.n_v4v5_fail)"' docs/methodology/_assets/20260627_subclone_4axis_teaching/data/full_v4v5_verification_HCC1395.json
```

**輸出**：stdout；**exit code**：0。

```text
funnel=113997/113997
tree_sum=13785/13785
recurrence=27 candidate=4 loh=13 cn_gain=10
v4v5_total=18931 capped_skipped=828 run=18103 pass=18103 fail=0
```

### C03 — source path、ignore 與 hashes

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
git check-ignore -v docs/methodology/_assets/20260627_subclone_4axis_teaching/data/funnel_census_HCC1395.json docs/methodology/_assets/20260627_subclone_4axis_teaching/data/ccf_and_cn_multisample.json docs/methodology/_assets/20260627_subclone_4axis_teaching/data/full_v4v5_verification_HCC1395.json
sha256sum docs/methodology/_assets/20260627_subclone_4axis_teaching/data/funnel_census_HCC1395.json docs/methodology/_assets/20260627_subclone_4axis_teaching/data/ccf_and_cn_multisample.json docs/methodology/_assets/20260627_subclone_4axis_teaching/data/full_v4v5_verification_HCC1395.json
```

**輸出**：stdout；**exit code**：0。

```text
.gitignore:10:data/ ... three files
d81f43... funnel_census_HCC1395.json
335d37... ccf_and_cn_multisample.json
892d43... full_v4v5_verification_HCC1395.json
```

### C04 — solver golden tests

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py
```

**輸出**：stdout；**exit code**：0。

```text
[PASS] two_siblings ...
[PASS] partial_IDP ...
GOLDEN 總結: ALL PASS ✓
```

### C05 — v2 synthetic tests

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_layered_reconstruction_v2.py
```

**輸出**：stdout/stderr + auto-cleaned `/tmp/tmp*/layered.json`；**exit code**：0。

```text
Ran 5 tests in 0.082s
OK
ResourceWarning: unclosed file ... layered_tree_reconstruction.py:238
ResourceWarning: unclosed file ... layered_tree_reconstruction.py:348
```

### C06 — C++ HP-tag regression subset

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
ctest --test-dir build -R '^ReadParserHPTagTest\.' --output-on-failure
```

**輸出**：stdout + `InterSubMod/build/Testing/Temporary/LastTest.log`；**exit code**：0。

```text
100% tests passed, 0 tests failed out of 12
```

### C07 — stress 前 100 案實際 checked/skipped 分拆

**輸入**：solver、legacy generator 邏輯、seed=20260710；**輸出**：stdout；**exit code**：0。

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 -c 'import random,itertools,sys; sys.path.insert(0,"docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"); import tree_enumeration_solver as S; rng=random.Random(20260710); out={"seed":20260710,"total":100,"checked_v4v5":0,"capped_not_tested":0,"checked_fail":0};
for _ in range(100):
 k=rng.randint(2,6); gs=["".join(b) for b in itertools.product("RA",repeat=k) if "A" in "".join(b)]; rng.shuffle(gs); nf=rng.randint(1,min(len(gs),2**k-1)); fl=gs[:nf]; full={g:rng.randint(3,40) for g in fl}; part=[]
 if rng.random()<0.5 and fl:
  base=rng.choice(fl); pos=rng.sample(range(k),rng.randint(1,max(1,k-1))); part=["".join("X" if i in pos else base[i] for i in range(k))]
 res=S.enumerate_min_trees(full,part,k)
 if res.get("capped"): out["capped_not_tested"]+=1; continue
 out["checked_v4v5"]+=1; out["checked_fail"]+=int(not S.verify_all(res,full,part,k,light=False)["overall"])
print(out)'
```

```text
{'seed': 20260710, 'total': 100, 'checked_v4v5': 98, 'capped_not_tested': 2, 'checked_fail': 0}
```

### C08 — V5 completeness adversarial counterexample

**輸入**：固定 seed=314159 產生第一個 non-capped `n_trees>32` case；**輸出**：stdout；**exit code**：0。

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 -c 'import random,itertools,sys,json; sys.path.insert(0,"docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"); import tree_enumeration_solver as S; rng=random.Random(314159)
for case in range(1,1001):
 k=rng.randint(3,6); gs=["".join(b) for b in itertools.product("RA",repeat=k) if "A" in "".join(b)]; rng.shuffle(gs); nf=rng.randint(1,min(len(gs),2**k-1)); full={g:rng.randint(3,40) for g in gs[:nf]}; part=[]
 if rng.random()<0.5:
  base=rng.choice(list(full)); pos=rng.sample(range(k),rng.randint(1,max(1,k-1))); part=["".join("X" if i in pos else base[i] for i in range(k))]
 res=S.enumerate_min_trees(full,part,k)
 if not res.get("capped") and res.get("n_trees",0)>32:
  ver=S.verify_all(res,full,part,k,light=False); print(json.dumps({"seed":314159,"case":case,"k":k,"n_trees":res["n_trees"],"stored":len(res["trees"]),"trees_complete":res["trees_complete"],"V5":ver["V5"],"overall":ver["overall"]},sort_keys=True)); break'
```

```text
{"V5": [true, "independent=1119744 vs reported=1119744"], "case": 2, "k": 5,
 "n_trees": 1119744, "overall": true, "seed": 314159, "stored": 32,
 "trees_complete": false}
```

### C09 — 真實 HCC1395 tree-cap exposure

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
jq -r '[.detail[].n_trees // 0] as $n | "units=\($n|length) over32=\([$n[]|select(.>32)]|length) max=\($n|max) sum=\($n|add)"' docs/methodology/_assets/20260618_subcluster_pilot/layered_reconstruction_HCC1395.json
```

**輸出**：stdout；**exit code**：0。

```text
units=18931 over32=371 max=125 sum=66837
```

### C10 — legacy 7-dataset funnel

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
jq -r '.samples | to_entries[] | [.key, .value.funnel.check_ok, .value.funnel.L1_all_pass_universe, .value.funnel.L6_retained_sSNV, .value.tree_level.determined_mutation_bearing, .value.tree_level.multi_hp_both_mutation_bearing, .value.has_cn] | @tsv' docs/methodology/_assets/20260627_subclone_4axis_teaching/data/funnel_census_7samples.json
```

**輸出**：stdout；**exit code**：0；七列皆 `check_ok=true`，但只有 HCC1395 的 legacy `has_cn=true`。

### C11 — v2 input lock

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
jq -r '"schema=\(.schema_version) datasets=\(.dataset_count) biological_samples=\(.biological_sample_count) all_pass=\(.all_pass)", "cn_available=" + ([.samples[]|select(.cn.status != "unavailable")]|length|tostring) + "/" + (.dataset_count|tostring)' research/20260710_layered_reconstruction_v2/input_manifest_lock.json
```

**輸出**：stdout；**exit code**：0。

```text
schema=2.0 datasets=7 biological_samples=6 all_pass=true
cn_available=5/7
```

### C12 — active v2 run state（稽核截止）

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
sed -n '1,120p' /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/run_status.tsv
```

**輸出**：stdout；**exit code**：0。

```text
HCC1395            preflight PASS
HCC1395_DORADO     preflight PASS
HCC1395_DORADO     mlhp START
HCC1395            mlhp START
```

### C13 — active run code hash 與 manifest Git 狀態

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
sha256sum -c /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/code.sha256
git check-ignore -v docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json
git status --short -- research/20260710_layered_reconstruction_v2/input_manifest_lock.json docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json
```

**輸出**：stdout；**exit code**：0。

```text
layered_v2_input_manifest.json: OK
sm_linkage_genomewide.py: OK
sm_multilocus_combinations.py: OK
tree_enumeration_solver.py: OK
layered_tree_reconstruction.py: OK
build_region_view.py: OK
verify_layered_v2.py: OK
.gitignore:10:data/  .../layered_v2_input_manifest.json
?? research/20260710_layered_reconstruction_v2/input_manifest_lock.json
```

- **[O-L1] C13 只證明檢查當下 working files 與 run-start hashes 相同**；它不證明未 commit 內容可由 Git 還原，也不取代 run 完成後的 output verification。

## 8. Ambiguity checklist：最終說明書不得留白的定義

| 模糊詞 | 最終文件必須寫出的精確定義 | 目前判斷 |
|---|---|---|
| sample | dataset、biological sample、technical/pipeline replicate 三者分開 | `[F-L1]` 7 datasets / 6 biological samples。 |
| universe | ClairS paired PASS operational backbone，非完整 biological somatic truth | `[F-L1]` manifest/header 可確認。 |
| genome-wide | chr1-22 autosomal scope；chrX/Y excluded | `[F-L1]` 不可只寫「全基因組」。 |
| region | adjacent SNV gap ≤50 kb 的 connected component；total span 可 >50 kb | `[F-L1]` code contract。 |
| retained | 經 scope、singleton、densest-8、read support 後進 solver 的 SNV | `[F-L1]` 必附六層 funnel。 |
| lineage | v2 primary=mutation-bearing HP1/HP2；root-only=reference control；HP3=H3? auxiliary | `[F-L1]` 不可沿用 legacy `is_lineage`。 |
| determined | 在指定 mutation-state model 下唯一最小候選樹 | `[I-L1]` 不等於真實 biological clone determined。 |
| ambiguous | 要分 exact-tree multiplicity 與 topology-shape multiplicity | `[F-L1]` 兩者不可合併百分比。 |
| capped | search budget 未完成；另有 display cap，兩者完全不同 | `[F-L1]` `capped` 不得算 pass/determined。 |
| V1-V7 pass | 必須零 skipped/零 failed；capped 另列 n/a | `[F-L1]` v2 schema 已設計，full artifact 待驗。 |
| stress pass | 只允許 checked cases 作分母；skipped 必另列 | `[O-L1]` legacy 4,000 不符合。 |
| independent oracle | 不共用核心 candidate construction/helper，且比對 canonical output set | `[I-L1]` 現行 V4/V5 不足。 |
| CN neutral | 有 declared segmentation contract 支持；missing source=unavailable | `[F-L1]` SAVANA 未列區 neutral 語意仍待確認。 |
| CCF | 需 purity/CN correction；未修正只能稱 read allele fraction | `[F-L1]` legacy CCF 命名需降級。 |
| tie-broken | posterior heuristic ≥0.6；必附 TEMP/MARGIN 與候選集合 | `[F-L1]` 不是 truth accuracy。 |
| recurrence candidate | CN-neutral 尚未排除的 candidate，不等於已確認 convergent evolution | `[F-L1]` HCC 有 4 個 candidate，需 case review。 |
| compatible | observed pairwise no-four-gamete conflict；是否有完整 completion/witness 需另證 | `[F-L1]` 不可寫「676 全可解」。 |
| robust | 必須指定對哪個 perturbation/threshold robust | `[U-L5]` 目前無 sensitivity matrix。 |
| reproducible | commit/patch + input hash + code hash + params + output hash + command | `[F-L1]` v2 run 正建立；legacy HTML 不符合。 |
| validated | 必須標 code/internal/biological 哪一層 | `[F-L1]` biological 層尚未完成。 |
| all | 每次都附 scope、population、eligible/not-applicable 與 denominator | `[I-L1]` 禁止裸寫「全部通過」。 |

## 9. 最終報告發布前 gates

### 9.1 Engineering validation — 必須全部通過

1. **[U-L5] v2 run completion**：7/7 datasets `complete PASS`，且 `verification_summary.json/tsv`、`output.sha256`、`verification.sha256` 存在。
2. **[U-L5] v2 H4 gate**：所有 non-capped eligible units 的 V1-V7 zero skipped/zero fail；capped 只列 n/a，不列 pass。
3. **[U-L5] v2 H6 gate**：每個 non-capped unit `analysis_candidate_set_complete=true`、`analysis_trees_generated=n_trees`、exact shape 非 null。
4. **[U-L5] recurrence gate**：逐 unit 重算後三類加總與 HTML 完全一致；舊 20/7/0 必撤回。
5. **[U-L5] denominator gate**：所有百分比帶 `population_id`；7 datasets/6 biological samples；primary HP1/2 mutation-bearing。
6. **[U-L5] provenance gate**：final HTML 指向 immutable run ID，列 input/code/output hash 與 dirty/patch state。
7. **[U-L5] test integration gate**：v2 Python tests 納入可重複 test target；修正 ResourceWarning。

### 9.2 Methodological robustness — 發布方法結論前

1. **[U-L5] sensitivity**：MAPQ={10,20,30}、BASEQ={0,10,20}、MINREAD={2,3,5}、MAX_SNV={6,8}、TIER_R 至少兩側 perturbation；報 funnel、class、shape、recurrence、read-AF 的變化。
2. **[U-L5] V5 replacement**：對 tractable cases 以第二實作列出 canonical tree set；對 large/capped cases提供 proof/witness 或誠實 n/a。
3. **[U-L5] capped partial-data check**：pairwise test 之外產生 completion/witness；抽 tractable capped-like cases對 exact solver。
4. **[U-L5] caller/backbone sensitivity**：至少在 truth subset 或第二 caller 評估結構穩定性；ClairS PASS 只能叫 operational backbone。
5. **[U-L5] CN semantics**：SAVANA BED 的 uncovered interval contract 要有來源與機器檢查。

### 9.3 Biological validation — 不完成就維持 claim ceiling

- **[F-L1] 現階段可寫**：read-supported regional mutation-state trees under an explicit model；engineering/internal consistency after gates pass。
- **[F-L1] 現階段不可寫**：confirmed biological clones、true evolutionary history、CN-robust topology、7 independent samples、true convergence=0、single-bulk information-theoretic limit。
- **[U-L5] 升級需求**：single-cell、multi-region、orthogonal phasing/CN truth，或其他不共用同一 bulk read evidence 的驗證。

## 10. 建議的安全用語與禁止用語

| 情境 | 建議用語 | 禁止／需撤回用語 |
|---|---|---|
| V1-V7 | `[I-L1]`「non-capped eligible units 通過 implementation invariant suite」 | 「所有樹已由 independent oracle 完整驗證」 |
| stress | `[I-L1]`「N checked、K skipped capped、0 checked failures」 | 「4,000/4,000 全 pass」 |
| topology | `[I-L2]`「model-defined candidate topology」 | 「真實 topology determined」「資訊論極限」 |
| CN | `[I-L1]`「algorithm does not use CN at L1；CN-stratified robustness 待測」 | 「CN-robust」 |
| read AF | `[I-L1]`「read-AF ordering heuristic」 | 「CCF confirmation」 |
| cross-sample | `[I-L1]`「7 datasets from 6 biological samples」 | 「7 independent samples」 |
| capped | `[I-L1]`「observed pairwise compatible / exact completion not established」 | 「676/676 全可解」 |
| recurrence | `[I-L1]`「10 gain、13 LOH、4 neutral candidates（需重跑 v2確認）」 | 「20 gain、7 LOH、0 real」 |

## 11. Reader questions（主報告必回答）

1. **[U-L5] 最終採用的 immutable run ID 是哪一個？7/7 completion、runtime、peak RSS、output size 各多少？**
2. **[U-L5] v2 primary HP1/2 mutation-bearing 的精確 denominator、determined/ambiguous/capped/recurrence 分別多少？**
3. **[U-L5] 371 個 legacy `n_trees>32` units 在 v2 all-candidate 後，shape classification 改變多少？**
4. **[U-L5] 修正 stress counter 後，4,000 generated cases 中 checked、capped-skipped、failed 各多少？**
5. **[U-L5] 哪個真正獨立的 tree-set implementation/證明支持 V5？若沒有，最終文件是否已更名降級？**
6. **[U-L5] HCC1395 的 4 個 CN-neutral recurrence candidates 是哪些 region/family？人工 read/CN review 結果？**
7. **[U-L5] SAVANA CN BED 未列出的區段為何可視為 neutral？是 producer contract 還是本專案假設？**
8. **[U-L5] MAPQ20、BASEQ0、MINREAD3、MAX8、50kb gap 的 sensitivity 結果是否穩定？**
9. **[U-L5] HCC1395 與 DORADO 如何作 paired pipeline replicate，而非算成兩個 biological n？**
10. **[U-L5] legacy HTML 是否已加 superseded banner、移除 20/7/0、51.1%、2.0%、5,749 等未重驗數字？**
11. **[U-L5] final HTML 能否由一個 command 從 immutable JSON 重新產生，並在 footer 顯示 hashes？**
12. **[U-L5] 哪一份 orthogonal truth 能支持 clone-level claim？若沒有，是否全篇一致維持 regional mutation-state inference？**

## 12. 結論

- **[O-L1] 可保留的結論**：六層 funnel 與 legacy class 加總守恆；8 solver golden、5 v2 synthetic、12 HP-tag C++ tests 通過；v2 runner 已開始建立較完整 provenance package。
- **[O-L1] 必須修正的結論**：recurrence 20/7/0、V5「完整枚舉 oracle」、4,000 全 checked、capped 676 全可解、51.1%/2.0%/5,749 的發布級可信度、7 independent samples。
- **[I-L1] 最終 HTML 應把每個 claim 對到 run-scoped machine-readable source，並把 code correctness、internal consistency、biological validation 分成三個獨立 verdict；任何一層未完成都不得用「全部驗證完成」概括。**
