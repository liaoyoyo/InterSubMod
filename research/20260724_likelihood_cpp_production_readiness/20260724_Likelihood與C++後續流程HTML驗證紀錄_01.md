<!--
建立時間: 2026-07-24
目標: 紀錄 likelihood 教學 HTML 的公式重算、max_sets 數學語意、July raw-all 資料契約、C++ 能力邊界與瀏覽器 QA
處理範圍: HCC1395 July raw-all／chr1-22 compact extraction、既有 exact candidate 與 likelihood 程式、strict W C++；不等同全 7 datasets directed-topology production run
關聯檔案:
  - InterSubMod/research/20260724_likelihood_cpp_production_readiness/20260724_Likelihood與C++後續流程一步步教學_01.standalone.html
  - InterSubMod/research/20260724_likelihood_cpp_production_readiness/scripts/qa_likelihood_cpp_html.py
  - InterSubMod/research/20260724_likelihood_cpp_production_readiness/results/20260724_likelihood_cpp教學HTML_QA_01/qa_receipt.json
-->

# Likelihood 與 C++ 後續流程 HTML 驗證紀錄

用 PREP（Point → Reason → Evidence → Point）：

> **結論：教學 HTML 與其核心數字已通過公式重算、程式來源、receipt、C++ fixture test、61 項 Python focused tests及四種 viewport QA；但 production directed topology 仍是 0/7。**

任務類型：**B — Comprehensive validation**。服務 G3／G4／G5。

## 1. 範圍與 claim ceiling

本次驗證回答四件事：

1. R/A/X＋BQ＋mixture likelihood 的公式與教學數字是否可重算。
2. `max_sets=256` 是否有數學／生物推導，以及超限後應如何處理。
3. July raw-all 與 HCC1395 compact extraction 是否包含後續計算所需資料。
4. 現有 C++ 已完成哪些層、尚缺哪些層。

本次**不宣稱**：

- 全 7 technical datasets 已完成 directed topology；目前是 **0/7**。
- latent state mixture weight \(\pi\) 是真實 cellular clone fraction。
- winning vertex set 唯一即可推出 exact parent edge 唯一。
- 將前 256 個候選排名可代表完整候選 family。

## 2. Step → Verify

| 步驟 | 驗證方式 | 實際結果 |
|---|---|---|
| 1. 核對 likelihood 定義 | 讀取 ranker／optimizer，重算兩個 toy vertex sets | LL 與 \(\pi\) 全部一致；KKT gap < `1e-8` |
| 2. 核對 cap 語意 | 讀取 enumerator complete／stop contract | 達 256 後為 `complete=false`、`MAX_SETS_REACHED` |
| 3. 核對 hard family | 讀取 isolated perfect-family fast-path receipt與原始 manifest | k=m=11、full=0、partial=70；明示非production scope下family=122,281,152 |
| 4. 核對 July 資料 | 讀取 VCF／sidecar／producer、big7 pilot manifest與 schema 1.3 receipts | 22/22 extraction PASS；資料齊備，但big7 authority未升格 |
| 5. 核對 C++ 現況 | source audit＋fixture tests＋HCC1395 parity receipt | strict W層 PASS；candidate／likelihood／directed topology未實作 |
| 6. 驗證 HTML | Playwright 1440／1024／390／320、no-JS、print、互動 | `pass=true`；overflow／console error／外部 request皆0 |

## 3. Likelihood 數字重算

### 輸入

- 程式：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py`
- patterns：`RR×20`、`AR×48`、`AA×30`、`RA×2`，每個固定 call 的 BQ=30。
- 候選 \(V_1=\{RR,AR,AA\}\)、\(V_2=\{RR,RA,AA\}\)。

### 執行命令

```bash
PYTHONPATH=research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 - <<'PY'
from build_m2_patterns_and_rank import fit_quality_aware_mixture
patterns=[('RR',(30,30),20),('AR',(30,30),48),
          ('AA',(30,30),30),('RA',(30,30),2)]
for vertices in [(0,1,3),(0,2,3)]:
    fit=fit_quality_aware_mixture(patterns,vertices)
    print(fit.vertices, tuple(round(x,6) for x in fit.weights),
          round(fit.log_likelihood,6),
          fit.global_log_likelihood_gap_bound)
PY
```

### 實際輸出

```text
(0, 1, 3) (0.207912, 0.480141, 0.311948) -118.913973 9.642334930504148e-09
(0, 2, 3) (0.392125, 0.019683, 0.588192) -427.745135 3.261106940044556e-09
```

驗證：

- \(\Delta LL=308.831162\)，HTML 數字一致。
- 兩個 global LL gap 都小於 `1e-8`。
- 這是在比較 **mutation-state vertex set**；沒有評分相同 vertex set 內的 parent edge。
- \(\pi\) 是固定候選下的 maximum-likelihood latent-state weight估計，不是EM responsibility的條件期望，也不是cellular clone fraction。
- X在現行 conditional-on-observed-mask模型中令emission因子為1；這假設missingness可忽略。若缺值機制與state相關，不可宣稱一般性的無偏邊際化。

## 4. `max_sets=256` 的正式判讀

來源：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`

現行 defaults：

```text
max_vertex_sets = 256
solver_time_limit_seconds = 30.0  # 每一次 MILP
tie_tolerance = 1e-6
certificate_tolerance = 1e-8
```

數學判讀：

- `k=8 → 2^8=256` 只代表超立方體**可能 vertex 數**，不代表候選 vertex sets或 directed trees最多 256。
- 若 active dimension為 \(m\)、mandatory states為 \(f\)、minimum-extra為 \(h^*\)，候選集合的寬鬆組合上界可達：

\[
\binom{2^m-f}{h^*}
\]

- 所以 `max_sets=256` 是工程 guardrail，沒有統計、生物或 clone 數學最適性。
- `256 / 122,281,152 = 0.000209%`；而列舉順序不是 likelihood排名，前256個不是代表性抽樣。
- 每次 MILP 若都用滿30秒，256次就是7,680秒＝128分鐘；只有 cap 仍不能保證單一 unit快速完成。

正確 contract：

```text
F<cap且下一解證明INFEASIBLE，或有外部完成證書
  → candidate_generation_complete=true → 才能排名
達cap／timeout／solver failure → complete=false → ABSTAIN
```

現行enumerator若剛好列出 `F=256`，沒有第257次INFEASIBLE或外部completion certificate，仍保守回報incomplete。

因此 **cap單獨不能保證整批終止或產生完整receipt**。只有再搭配per-unit／global wall-time、memory cap、checkpoint、錯誤處理與fail-closed receipt，才能保證bounded termination並把超限／逾時unit回報為ABSTAIN；即使如此，也不能保證每個區域都有topology。

## 5. 122,281,152 hard family

### 權威輸入

- `InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_stress_panel_v1/authoritative_r3/manifest.json`
- `InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_stress_panel_v1/perfect_family_count_v3/receipt.json`

### 實際值

```text
dataset=HCC1395_DORADO
chrom=chr6
component_basis=PS_HP2
k=11
full patterns=0
partial patterns=70
reduction active groups q=5
ALT comparable pairs=29
forbidden ancestor pairs=0
minimum-extra objective=11
perfect family count=122,281,152
counter elapsed=0.044700817 s
ranking_allowed=false
scope=ISOLATED_PERFECT_EVENT_H_EQ_M_MINUS_D_FAST_PATH_NOT_PRODUCTION
```

這說明：

- 122,281,152只在上述recurrence-free、\(h^*=m-d\) isolated perfect-event fast-path模型中是exact。
- 它不是production M2已完整列舉或排序的候選數；`ranking_allowed=false`。
- `256 / 122,281,152 = 0.000209%`只能作為「cap最多允許materialize數／diagnostic family」的假設性數量比，不是實際覆蓋率或代表性抽樣。
- fast-path計數很快不代表能逐候選 likelihood fit。
- partial patterns對「哪些狀態可解釋 read」有資訊，但對突變先後方向很弱。
- 大量排列都具有相同 minimum-extra objective。

## 6. July raw-all 與 HCC1395 compact extraction

### July sample資料夾

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/
```

已觀察：

| 產物 | 實際值 |
|---|---:|
| normalized ClairS raw-all | 134,122 records |
| LongPhase-S recalibrated PASS | 113,061 records |
| mapped alignment rows | 40,859,727 |
| LongPhase-S tagged rows | 17,700,924 |
| exact duplicate identity conflicts | 0 |
| persisted tagged BAM | `false` |

July資料夾有 VCF、HP／PS sidecar與receipts，但沒有持久化完整 tagged BAM。producer receipt綁定big8 canonical raw tumor BAM；BQ、sequence、CIGAR與MAPQ原則上可從raw BAM讀取，再用sidecar join HP／PS。

### 最新 compact extraction

輸入根：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2/
```

聚合 receipt 實際輸出：

```json
{
  "receipt_count": 22,
  "all_pass_count": 22,
  "exact_join_pass_count": 22,
  "sidecar_missing_sum": 0,
  "n_sSNV_sum": 79687,
  "schema_versions": ["1.3.0"]
}
```

chr6 receipt：

```text
n_sSNV=27,099
canonical eligible alignments=416,944
fixed R/A calls=2,099,726
known PS × HP1/HP2 molecule rows=80,247
sidecar exact matches=496,928
sidecar missing=0
```

但22份schema 1.3 compact extraction綁定的是big7 copy：

```text
/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam
```

July producer綁定的canonical路徑則是：

```text
/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam
```

big7 pilot manifest明示：

```text
task_type=exploratory_pilot
claim_status=PARTIAL
validation_evidence_eligible=false
alignment_transfer_assurance=same_size_plus_fixed_first_middle_last_chunk_sha256_plus_full_bai_sha256
is_full_bam_content_hash=false
```

因此必須分開記錄三層判定，不能把content equivalence與production authority合併：

1. **Gate A — bounded content equivalence：PASS。**22/22 × 6 compact artifacts語意相同；BAM的size、三段固定抽樣chunk與完整BAI相符。但這不是full-BAM hash，也不是完整identity證書。
2. **Gate B — production authority：FAIL／NOT GRANTED。**big7 manifest仍是`PARTIAL`且`validation_evidence_eligible=false`。須用big8 canonical BAM重抽，或另建正式authority contract；Gate A通過不能推導Gate B通過。
3. **Gate C — adapter/software readiness：MISSING。**需建立schema `1.3.0` calls＋strict W membership＋k≤12 partition的composite adapter與passing receipt。現行M2 ranker仍要求schema `1.2.0`同目錄四檔；不可只改版本常數。

結論：**資料種類齊備且bounded content equivalence通過；production authority與adapter仍不足。**

## 7. C++ 能力邊界

### 已完成

- exact `chromosome × PS × HP` strict W。
- fixed R/A endpoint read support。
- k≤12 partition。
- HCC1395 chr1–22 Python／C++ parity：

```text
all_pass=true
components=39,846
observed edge rows=117,760
component mismatch=0
edge mismatch=0
```

### 本次 C++ fixture test

輸入：

```text
InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/cpp/tests/
```

命令：

```bash
bash research/20260723_production_exact_ps_strict_read_linkage/cpp/tests/run_tests.sh
```

輸出片段：

```text
test_strict_endpoint_graph: PASS
strict_endpoint_graph_verify: PASS
run_tests.sh: PASS
```

### 尚未完成

1. R/A/X subcube group constraints 的 C++ adapter。
2. minimum-extra exact candidate solver與完整性證書。
3. BQ-aware emission likelihood。
4. simplex mixture optimizer與KKT gap。
5. winner／tie分類。
6. Hamming-1 directed parent assignment與Topo分類。
7. production fail-closed receipts與Python parity。

source audit顯示目前 `Observation`只有 `node`與`call_code`，沒有BQ；建立的是fixed endpoint無向共現邊與connected components，不是演化 parent edge。

### Python focused regression

為避免修改專案環境，pytest只安裝在 `/tmp/intersubmod_pytest`；測試命令：

```bash
PYTHONPATH=/tmp/intersubmod_pytest:research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q \
research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py \
research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_m2_patterns_and_rank.py
```

實際輸出：

```text
............................................................. [100%]
61 passed in 2.79s
```

## 8. HTML QA

### 輸入

`InterSubMod/research/20260724_likelihood_cpp_production_readiness/20260724_Likelihood與C++後續流程一步步教學_01.standalone.html`

### 執行命令

```bash
python -m py_compile \
research/20260724_likelihood_cpp_production_readiness/scripts/qa_likelihood_cpp_html.py

python \
research/20260724_likelihood_cpp_production_readiness/scripts/qa_likelihood_cpp_html.py
```

### 輸出

`InterSubMod/research/20260724_likelihood_cpp_production_readiness/results/20260724_likelihood_cpp教學HTML_QA_01/`

最終 receipt：

```text
pass=true
HTML SHA256=ec4089b80e4354356e00dbaf51ad822eb23db92f107a23eb536c3a3250f11f77
viewports=1440,1024,390,320
horizontal overflow=0/4
console errors=0
page errors=0
external requests=0
no-JavaScript overflow=false
print details=6/6 visible
```

互動斷言：

- BQ=10 → match=`0.964286`、flip=`0.035714`。
- 切換 \(V_2\) → LL=`−427.745135`。
- cap=8,192 → `0.006699%`。
- 展開／收合 6/6 details。
- paper → night主題切換成功。

## 9. 多視角審查

| 視角 | 審查結果 |
|---|---|
| 數學／cap／C++／資料契約 | 256不是數學上限；cap須fail-closed；July資料種類齊備，但BAM authority與軟體層未完成 |
| likelihood公式／措辭 | 正名為vertex-set likelihood；X邊際化；BQ公式、tie與KKT數字一致 |
| HTML資訊架構／可讀性 | 採研究推論帳本、evidence badges、sticky目錄、摺疊細節 |
| 瀏覽器／視覺 | 四種viewport、no-JS、print與互動均PASS |

## 10. 最終判定

1. HTML 可以用來一步步向教授解釋 likelihood 與現行 production邊界。
2. `max_sets=256` 可暫留為工程預設，但它單獨不保證終止；必須再加per-unit／global wall-time、memory、checkpoint、stop reason與fail-closed ABSTAIN receipt。
3. 若希望救回1.22億candidate family，必須做 compressed exact representation／certified branch-and-bound／BDD或ZDD類方法；調高cap不是解法。
4. C++可完成後續步驟，但**現有binary不能直接產生正式全部結果**。
5. 最快且較低風險的路線：

```text
C++ strict W
→ big8重抽或正式big7↔big8 authority contract
→ schema 1.3 + strict W composite adapter
→ compressed/certified exact candidate backend
→ 先沿用Python likelihood作oracle
→ Python/C++ parity
→ 再決定是否移植likelihood
→ directed topology status + fail-closed receipt
```
