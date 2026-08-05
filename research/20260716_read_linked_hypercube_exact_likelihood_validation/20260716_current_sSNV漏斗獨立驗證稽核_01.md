<!--
建立時間: 2026-07-16 10:05 +08:00
目標: 不共用 current-funnel producer 函式，獨立重算並驗證最新 7-dataset sSNV 漏斗 receipt
處理範圍: Task Type B；7 technical datasets / 6 biological samples；chr1-chr22 自 autosomal 層開始
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_current_funnel_receipt.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_current_funnel_receipt.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/current_funnel_receipt.json
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/independent_verification.json
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
-->

# Current sSNV 漏斗獨立驗證稽核

## TL;DR

**漏斗 receipt 可引用：獨立 verifier 由 7 份 frozen site-ledger 直接重算，322/322 checks PASS，確認
`638,259 → 582,820 → 469,849 → 194,149`；其中 autosomal `469,849 = 50,432 singleton +
225,268 MAX_SNV-excluded + 194,149 retained`。（影響：高；信心：高）**

本 verdict 只證明「最新 sSNV 漏斗數字、資料範圍、分母、來源 SHA 與守恆」可用；它不替代尚未完成的
M2 全 154 dataset×chromosome read-likelihood 結果驗證。

## 1. 問題、單位與範圍

- 服務目標：**G4 多樣本一致性與 reproducibility、G5 外部可驗證**。
- Task Type：**B Comprehensive validation**。
- 計數單位：**dataset-record / dataset-site**，不是 unique genomic locus、cell、clone 或 biological sample。
- 範圍：7 technical datasets、6 biological samples；`HCC1395` 與 `HCC1395_DORADO` 是同一 biological
  sample 的兩份 technical datasets，彙總時依設計各自計數。
- 染色體：`638,259` raw 與 `582,820` LongPhase-S PASS 是全 producer contract；從 `469,849`
  開始才是 **chr1–chr22、biallelic sSNV** scope。
- `PASS` 是 LongPhase-S recalibrated FILTER 狀態，不等同 truth-confirmed somatic mutation。
- `MAX_SNV-excluded` 是 site-ledger 的分支名稱：位置群組受 `MAX_SNV=8` 規則處理後未 retained 的位點；
  不應解讀為該位點不存在或沒有 read evidence。

## 2. 獨立驗證設計（Step → Verify）

1. 直接讀 canonical JSON、7 份 site-ledger summary 與 receipt under test
   → 驗證：實際路徑和 SHA-256 都與 receipt binding 相同。
2. 逐 dataset 檢查 ledger schema、raw-all / recalibrated-PASS contract、7 個 producer checks、四層 record-key
   duplicate excess
   → 驗證：7/7 `pass=true`、所有 producer checks=true、duplicate excess 全為 0。
3. 不呼叫 producer 函式，直接重算五個互斥 branch
   → 驗證：每樣本均滿足：
   - `raw = excluded_by_longphase_filter + out_of_scope_non_autosomal + max_snv_excluded + positional_singleton + retained`
   - `LongPhase-S PASS = raw − excluded_by_longphase_filter`
   - `autosomal = LongPhase-S PASS − out_of_scope_non_autosomal`
   - `autosomal = max_snv_excluded + positional_singleton + retained`
4. 把逐 dataset 重算值與 canonical 的 tree/autosomal/retained 欄逐一對齊
   → 驗證：7/7 dataset、3 個 canonical counts 全相等；aggregate 亦全相等。
5. 把重算值與 current funnel receipt 的 counts、branches、relative ratios、total ratios、source path、source SHA
   逐欄比較
   → 驗證：322/322 checks PASS，0 failures。
6. 故障注入測試
   → 驗證：篡改 aggregate、破壞 source conservation/SHA、製造 dataset duplicate/order drift、篡改 ratio 均會
   fail closed。

獨立性界線：verifier 僅使用 Python standard library，**未 import、未 call**
`build_current_funnel_receipt.py`；測試 fixture 亦直接定義 frozen source counts，而非用 producer builder 生成。

## 3. 最終漏斗結果

| 步驟 | 數量 | 相對前一層 | 占 raw | 意義 |
|---|---:|---:|---:|---|
| ClairS normalized raw-all | 638,259 | 100.000% | 100.000% | 7 datasets 的 producer raw records |
| LongPhase-S recalibrated PASS | 582,820 | 91.314% | 91.314% | 排除 FILTER 未通過 55,439 records |
| chr1–22 biallelic sSNV | 469,849 | 80.616% | 73.614% | 再排除 non-autosomal/out-of-scope 112,971 records |
| retained sSNV | 194,149 | 41.322%（以 autosomal 為分母） | 30.419% | 可進入目前位置群組/樹流程的 retained 位點 |

### 3.1 Autosomal 層的完整互斥分類

| Autosomal 分支 | 數量 | 占 469,849 | 守恆角色 |
|---|---:|---:|---|
| positional singleton | 50,432 | 10.734% | 單一位點位置群組，未進多點樹 |
| MAX_SNV-excluded | 225,268 | 47.945% | `MAX_SNV=8` 規則下未 retained 位點 |
| retained | 194,149 | 41.322% | 進入後續 retained scope |
| **合計** | **469,849** | **100.000%** | **50,432 + 225,268 + 194,149 = 469,849** |

另兩個上游互斥分支是：LongPhase-S FILTER excluded `55,439`、non-autosomal/out-of-scope `112,971`；
五個 raw branches 合計：
`55,439 + 112,971 + 225,268 + 50,432 + 194,149 = 638,259`。

## 4. 每 dataset 重算結果

| Dataset | raw | FILTER excluded | LPS PASS | non-auto | autosomal | singleton | MAX-excluded | retained | retained / autosomal |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| COLO829 | 43,014 | 3,556 | 39,458 | 1,670 | 37,788 | 7,830 | 2,202 | 27,756 | 73.452% |
| H1437 | 83,261 | 3,606 | 79,655 | 2,575 | 77,080 | 6,696 | 32,841 | 37,543 | 48.707% |
| H2009 | 168,638 | 7,043 | 161,595 | 7,130 | 154,465 | 2,853 | 100,448 | 51,164 | 33.123% |
| HCC1395 | 134,122 | 21,061 | 113,061 | 33,374 | 79,687 | 8,279 | 44,306 | 27,102 | 34.011% |
| HCC1395_DORADO | 123,240 | 10,095 | 113,145 | 33,406 | 79,739 | 8,321 | 43,803 | 27,615 | 34.632% |
| HCC1937 | 59,161 | 7,046 | 52,115 | 33,425 | 18,690 | 8,469 | 890 | 9,331 | 49.925% |
| HCC1954 | 26,823 | 3,032 | 23,791 | 1,391 | 22,400 | 7,984 | 778 | 13,638 | 60.884% |
| **合計** | **638,259** | **55,439** | **582,820** | **112,971** | **469,849** | **50,432** | **225,268** | **194,149** | **41.322%** |

每一列都通過上述四條守恆式；不是先彙總後才勉強相等。

## 5. 輸入、命令、輸出與實際結果

### 5.1 輸入

- Canonical summary：
  `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`
- 7 份 ledger root：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`
- Receipt under test：
  `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/current_funnel_receipt.json`

### 5.2 執行命令

工作目錄：`/big7_disk/liaoyoyo2001/InterSubMod`

```bash
/usr/bin/time -v python3 \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_current_funnel_receipt.py \
  --canonical-json research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json \
  --ledger-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 \
  --receipt research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/current_funnel_receipt.json \
  --output research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/independent_verification.json
```

測試命令：

```bash
python3 -m unittest -v \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_current_funnel_receipt.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_current_funnel_receipt.py
```

### 5.3 輸出

- 獨立 verifier：
  `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_current_funnel_receipt.py`
- 故障注入測試：
  `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_current_funnel_receipt.py`
- Machine-readable 結果：
  `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/independent_verification.json`

實際輸出片段：

```text
all_pass: true
n_checks: 322
n_fail: 0
raw_clairs_records: 638259
longphase_s_recalibrated_pass: 582820
autosomal_biallelic_sSNV: 469849
autosomal_partition_identity: 50432 + 225268 + 194149 = 469849
retained_sSNV: 194149
Elapsed: 0.16 s
Maximum resident set size: 16,428 kB
Exit status: 0
```

測試實際結果：`8 tests / 8 PASS / 0 fail`，exit status `0`。其中 verifier 專屬測試 5 個，producer
builder 既有測試 3 個。

## 6. SHA-256 證據鏈

| Artifact | SHA-256 |
|---|---|
| canonical JSON | `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7` |
| current funnel receipt | `a6282ecb73df314782d39a6ae4df410cebca05a42425bd7eebe28f12b1d35d75` |
| independent verifier | `bfa2c25137ec81707d47608936ca0643809602da0123e6b2fc000517ec731c15` |
| verifier tests | `6ee643d7bbeabe0b27c1962d3d4d915c1b1840df9ea4fe80458de30aa4b70f41` |
| independent verification JSON | `a0a098f103980204269d92a1d75ac148b408be3c0ea408349a1a84f7487eb796` |

Site-ledger summaries：

| Dataset | SHA-256 |
|---|---|
| COLO829 | `11b77ed47f8dd06b09e7a3cc0888fc3cc84031141c6b6fb30831517cbdf97612` |
| H1437 | `58e44c81d6fa372ef83a57efa3aff19aeb582200af9e3b0fc99ce369489412c2` |
| H2009 | `6012d222e3f5ea0fa923a5442f4fe87f633650fb121c782977c19a919a5b883e` |
| HCC1395 | `94dc3f49cf28fffb11e1daec7498b5750a9dc97106efccb39d4070e8c232e9bb` |
| HCC1395_DORADO | `a67643ebf997c4587bc5a0ab4b802a8ec12e7ab72955138ce1dde15ddce13a50` |
| HCC1937 | `01a90c2fb9cf4bb1627878f5030c139cc51f7e4dfd8a00d3b1fcb5a401793539` |
| HCC1954 | `44c672b3a1d185394e2e0ba1f6c5023640f645026f3bcb8d0d283fe0d1d54564` |

## 7. Validation verdict 與限制

### Overall Assessment：Ready to share（限漏斗數字）

沒有發現 count、denominator、branch conservation、dataset coverage、canonical reconciliation、ratio 或 source binding
錯誤。這組數字可以放入教授報告，前提是同時保留以下標籤：

1. `7 technical datasets / 6 biological samples`，不能寫成 7 biological samples。
2. chr1–chr22 scope 從 `469,849` 層才開始。
3. `194,149` 是 retained sSNV 位點數，不是 region、tree、candidate、clone 或 subclone 數。
4. `PASS` 是 LongPhase-S filter contract，不是 somatic truth verdict。
5. 本稽核不表示 M2 full extraction/ranking 已完成；最終整份研究報告仍須等待該 terminal receipts 與獨立 verifier。
