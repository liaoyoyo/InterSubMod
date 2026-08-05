<!--
建立時間: 2026-07-24
目標: 將 7/22 exact-PS、k<=12 分區接到 C++ 最簡樹重建與 read-AF 排序，先驗證 HCC1395，再以同一契約跑完整七資料集
處理範圍: Task Type B / comprehensive validation / HCC1395 gate then all seven datasets
關聯檔案: InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/pre-decision-audit.md; InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/implementation-notes.md
-->

# Exact-PS C++ topology + read-AF 全資料集驗證

## 目標

1. 凍結 2026-07-22 exact `PS × HP × read-linked component × k<=12 block` 輸入契約。
2. 以舊 Python `Fraction` read-AF score 作為 oracle，證明或否證固定 vertex set 下的 parent-factorization。
3. C++ 完成最小 vertex-set family、Hamming-1 parent assignment、AF 最佳分數／並列數／唯一性與 topology shape。
4. 先在 HCC1395 chr1–22 做逐 unit exact parity 與 controlled benchmark。
5. HCC gate 全通過後，才執行七資料集 chr1–22，產生壓縮輸出、receipt、時間與數據摘要。

## Claim ceiling

輸出是 read-compatible、recurrence-allowed、unit-cost 最簡 mutation-state
arborescence 的 AF 排序，不是已證實的 cell clone、真實祖先關係、CCF 或
CN/LOH-corrected phylogeny。任何 cap、deadline、輸入缺失或 family 未列完均
必須 `ABSTAIN`，不得產生 winner。

## 2026-07-24 結果

- HCC1395 oracle gate：11,122/11,122 complete units及4,601 read-AF units
  皆0 mismatch；C++ 6.73 s對Python 101.54 s，約15.09x。
- 七樣本chr1–22：technical PASS；98,955 groups、71,955 ranked、
  39,648 unique winners、32,307 tied。
- 10,717 mutation-bearing units命中guard而保守棄權，主要為H2009及
  active k>=6；因此全cohort尚非topology-complete。
- 完整證據與命令：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/20260724_C++_exact_PS拓撲與read_AF全七樣本驗證_01.md`。
