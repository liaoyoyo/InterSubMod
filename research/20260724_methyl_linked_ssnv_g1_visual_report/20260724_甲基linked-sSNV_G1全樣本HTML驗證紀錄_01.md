<!--
建立時間：2026-07-24
目標：建立全樣本 methyl-group × linked-sSNV G1 關聯 HTML，回答可驗證位置、甲基差異、read-state 結構與證據限制
處理範圍：Task B comprehensive validation；全基因組；7 technical datasets／6 biological samples
關聯檔案：20260724_methyl_linked_sSNV_G1全樣本關聯圖譜_01.standalone.html、data/report_data.json、results/20260724_G1全樣本HTML_QA收據_01.json
服務目標：G3 read-level epigenetic；G4 多樣本一致性與 reproducibility
-->

# 甲基群 × linked-sSNV G1 全樣本 HTML 驗證紀錄

## TL;DR

全範圍重算確認 **7 個正式 G1 pair／7 個 focal sites**：可支持「在 focal-ALT molecules 內，甲基群 MG 與某個 linked partner sSNV allele 顯著共分離」；**0/7 可確認 mutation order，G2=0，不能升級成因果、cellular clone/subclone 或唯一演化樹**。

## 任務分類與關鍵假設

- Task type：B comprehensive validation，不使用 subset。
- Thread D：是；屬 read-level epigenetic。
- Thread B 撤回範圍：否。
- KDE-corrected：本任務不使用 KDE/VAF；使用 5mC read×CpG matrix、fixed read allele calls 與 exact-PS×HP linkage。
- VCF caller AF：不使用。
- 長計算：只做權威結果重整、矩陣雜湊、HTML 建置與瀏覽器 QA；未重跑 BAM 全基因組 caller。

最重要的語意假設：

1. `MG-1-1` 等是 methyl group，不是 HP tag。
2. G1 的 partner 是「與 MG 共分離的 linked sSNV」；不是被證明造成甲基差異的 causal mutation。
3. W direct edge 是無向 read linkage；不是父子演化箭頭。
4. v8 `EXECUTION_PASS` 只代表 producer execution integrity；仍非 signed Task-B final scientific release。

## Step → Verify

1. 讀取全 407,738 directed pair rows，篩正式 G1  
   → 驗證：正式 G1=`7`、global BY=`10`、exact-family=`147`。
2. 串流全 102,842 stable assignment records，join 7 個 focal sites  
   → 驗證：7/7 assignment 存在；各 matrix SHA 與 assignment identity 相同；826/826 stable-core read IDs 可映射。
3. 讀取 current exact-PS×HP threshold=3 TSV  
   → 驗證：7/7 pair 同 W；10 個 W containers；10/10 direct edge；support 合計 790。
4. 建立 standalone HTML 與 machine-readable JSON  
   → 驗證：7 張完整 locus 卡、7 張 all-core-read heatmap、7 個 focal-first state graph。
5. Playwright 桌機與手機 QA  
   → 驗證：desktop/mobile 都顯示 7 cards、7 index rows、7 canvases；console errors=0；external requests=0；頁面無水平 overflow。

## 全樣本漏斗

| Dataset | 全 sSNV | M1 | M2 evaluable | M2 PASS | 有 partner | Exact focal | Exact pairs | BY | G1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 79,687 | 12,838 | 37 | 4 | 2 | 2 | 2 | 1 | 1 |
| HCC1395_DORADO | 79,739 | 14,789 | 68 | 13 | 1 | 0 | 0 | 0 | 0 |
| COLO829 | 37,788 | 3,579 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| H1437 | 77,080 | 10,187 | 41 | 16 | 5 | 1 | 2 | 0 | 0 |
| H2009 | 154,465 | 54,644 | 1,603 | 816 | 507 | 122 | 138 | 8 | 5 |
| HCC1937 | 18,690 | 1,938 | 86 | 50 | 19 | 2 | 4 | 0 | 0 |
| HCC1954 | 22,400 | 4,867 | 32 | 20 | 14 | 1 | 1 | 1 | 1 |
| **總計** | **469,849** | **102,842** | **1,867** | **919** | **548** | **128** | **147** | **10** | **7** |

主要比例：

- M1：102,842 / 469,849 = 21.8883%
- M2 PASS：919 / 469,849 = 0.1956%
- M2 PASS / evaluable：919 / 1,867 = 49.2234%
- G1 / exact-family pairs：7 / 147 = 4.7619%
- G2 global BY：0 / 58

## 哪些 linked 位點與 MG 關聯最強

下表依 Cramér’s V 排序；V 描述關聯效果量，不表示因果或演化先後。

| 排名 | Dataset | Focal → linked partner | V | ΔALT | BY q | Order |
|---:|---|---|---:|---:|---:|---|
| 1 | H2009 | chr4:2,307,521 T>C → 2,304,921 G>T | 0.964 | 94.7% | 5.54e-13 | 未定 |
| 2 | H2009 | chr12:4,413,414 C>A → 4,414,974 G>C | 0.929 | 87.5% | 2.41e-7 | 未定 |
| 3 | H2009 | chr18:567,920 T>C → 563,687 A>G | 0.894 | 88.9% | 1.35e-14 | 未定 |
| 4 | HCC1954 | chr8:100,071,382 A>G → 100,070,832 C>T | 0.868 | 86.8% | 5.54e-13 | 未定 |
| 5 | H2009 | chr13:28,815,939 G>A → 28,817,498 G>A | 0.847 | 81.7% | 1.61e-7 | 未定 |
| 6 | HCC1395 | chr3:127,986,757 G>A → 127,981,978 C>G | 0.645 | 65.6% | 1.08e-4 | 未定 |
| 7 | H2009 | chr5:44,615,693 G>T → 44,612,223 G>C | 0.614 | 50.9% | 8.97e-9 | 未定 |

H2009 chr12 的 partner truth label 是 FP；其 G1 association 可作 read-level 統計現象，但不能把 partner 自動升級成真 somatic mutation。

## HTML 的四層視覺

每個位置皆包含：

1. genomic locus track：focal 與 partner 的座標、allele、距離；線段明示為無向。
2. MG × partner allele stacked bars：顯示各 MG 的 partner REF/ALT count 與 ALT fraction。
3. 5mC read×CpG heatmap：保留全部 stable focal-ALT core reads；顯示通過 per-group coverage gate 後，群間平均差最大的前 40 CpGs。
4. focal-first observed state graph：RR、RA、AR、AA、X；旁列三種未確認候選形狀與 Endpoint-B 失敗理由。
5. exact-PS×HP W：PS、HP、W k、component span、direct support。
6. supported／not-supported claim cards與 legacy/posthoc audit 欄位。

![HTML 首頁與重點答案](results/01_G1_atlas_top_desktop.png)

![HCC1395 完整 locus 卡與 read×CpG 熱圖](results/02_G1_HCC1395_case_desktop.png)

![手機版首頁](results/03_G1_atlas_mobile.png)

## 輸入、命令、輸出

### 權威輸入

- Pair statistics：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity/methyl_ssnv_pair_results.tsv.gz`
- Whole-scope summary：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity/summary.json`
- Stable assignments：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_stable_multigroup_read_assignments.jsonl.gz`
- Strict linkage root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/`
- Claim contract：`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md`

### 建置命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260724_methyl_linked_ssnv_g1_visual_report/scripts/build_report.py
```

實際輸出片段：

```json
{"status":"PASS","cases":7,"core_reads":826,"strict_w_containers":10}
```

### QA 命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260724_methyl_linked_ssnv_g1_visual_report/scripts/qa_report.py
```

實際輸出片段：

```json
{"status":"PASS","desktop_cases":7,"mobile_cases":7,"console_errors":0,"external_requests":0}
```

### 主要輸出

- [Standalone HTML](20260724_methyl_linked_sSNV_G1全樣本關聯圖譜_01.standalone.html)
- [Machine-readable report data](data/report_data.json)
- [建置收據](results/20260724_G1全樣本HTML建置收據_01.json)
- [瀏覽器 QA 收據](results/20260724_G1全樣本HTML_QA收據_01.json)

## 最終解釋

目前可以嚴格敘述：

> 在 469,849 個全樣本 sSNV 中，7 個 focal-ALT 位點的穩定 5mC methyl groups，與各自一個 linked partner sSNV allele 在相同 molecules 中顯著共分離；7/7 亦有 current exact-PS×HP direct read linkage。

目前不可以敘述：

> 這 7 個 partner mutation 造成甲基差異、代表 7 個 subclones、或已得到 7 棵準確且唯一的真實演化樹。

因此，甲基資訊在這裡的可靠用途是「標記局部遺傳–表觀狀態共分離」，並提供值得後續 single-cell／colony／multi-region 驗證的候選位置。
