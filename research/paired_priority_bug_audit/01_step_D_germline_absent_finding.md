<!--
build_date: 2026-05-09
agent: paired germline-absent xref pilot (Step D extension)
status: validated
report_class: comparative-mechanism-audit
audience: PI / lab member / 自己未來
parent_audit: InterSubMod/research/paired_priority_bug_audit/00_audit_report.md
inputs:
  - /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/.../longphase_s/HCC1395_tagged.bam
  - InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/vote_dump_{baseline,v3f,v5}_chr19.tsv.gz
outputs:
  - 本檔
  - step_D_germline_absent_xref/cross_tab_summary.txt
  - step_D_germline_absent_xref/paired_chr19_read_hp.tsv
verdict: NEW FINDING — V5 Layer 1.5 在 germline-absent 區域繼承 priority bug 偏移（與 baseline 4.19:1 偏 HP1 完全相同）；V3F 全標 hp=33 反而是更正確的保守處理
last_verified: 2026-05-09
report_template: results-report v1.0
-->

# Step D: Germline-absent 區域 Paired × TO Cross-Reference Audit — V5 Layer 1.5 的盲點

## 0. TL;DR

> 用戶 5/9 進階提問：「paired 中 HP1-1+HP2-1 同位點共現但**沒有 HP1/HP2** 的 reads，TO mode (longphase-to) 是否會錯標都偏向 HP1-1？」
>
> **發現**：是。在 chr19 germline_vote=0 的 5,789 events 內，TO baseline 標 hp=11 (HP1 系列) vs hp=21 (HP2 系列) = **3,312 : 791 = 4.19:1 偏 HP1**（priority bug 偏移）。**V3F 全標 hp=33（純 somatic 方向不定）保守正確；但 V5 Layer 1.5 與 baseline 行為完全相同（同 4.19:1 偏移）— V5 Layer 1.5 在 germline-absent 區域實質上是 priority bug 的 feature 化，非修補**。

---

## 1. 動機

5/9 用戶提問核心：

> 「是否是出現在 paired tag 成 HP1-1 與 HP2-1 同時出現，**沒有 HP1 與 HP2 的位點與 read 部分**，造成 longphase-to tag 錯誤會有發生那種都偏向於 HP1-1 的問題造成標記有問題，還因此標成錯誤的 HP1, HP1-1 或是 HP2 與 HP2-1 這樣的狀況」

對應 V5 audit 已知：T1.2-F1 全基因組 21.7M reads 中 germline_vote=0 / V5 Layer 1.5 觸發 +560,881 reads。**germline 缺席區的 priority bug 是否被修對？V5 Layer 1.5 是修補還是 feature 化？**

## 2. 方法

### 2.1 Cross-reference 設計

對 chr19 同一 read（用 read_name 對齊）：
1. paired tagged BAM (`HP:Z:`) — longphase-s 的 tag 結果（**獨立 ground truth**）
2. T1.2 baseline / V3F / V5 vote dump (`hpResult`) — longphase-to 的 tag 結果

JOIN 後計算每個 (read, position) event 的 (paired_HP × T1.2_hpResult) cross-tab。

### 2.2 Germline-absent events 篩選

從 T1.2 vote dump 篩出 `cnt_HP1 + cnt_HP2 = 0`（germline 缺席）且 `somatic_total > 0` 的 events。對這些 events 做 cross-tab。

### 2.3 資料量

| 資料 | events |
|---|---:|
| paired chr19 reads (HP:Z: tagged primary) | 295,856 reads |
| T1.2 baseline chr19 dump | 549,206 events / 462,435 unique reads |
| paired ∩ T1.2 overlap | ~340K events |
| **germline-absent (cnt_HP1+cnt_HP2=0) events with paired covered** | **5,789** |

## 3. 結果

### 3.1 Germline-absent events cross-tab

| paired HP | events | baseline hp=11 (HP1 系列) | baseline hp=21 (HP2 系列) | V3F hp=33 | V5 hp=11 | V5 hp=21 |
|---|---:|---:|---:|---:|---:|---:|
| HP:Z:1-1 | 2,040 | **1,679** | 318 | 2,040 | 1,679 | 318 |
| HP:Z:2-1 | 1,588 | **1,291** | 295 | 1,588 | 1,291 | 295 |
| HP:Z:3 | 530 | 342 | 178 | 530 | 343 | 177 |
| **加總（hp 系列）** | — | **3,312** | **791** | **2,158**（純 somatic）| 3,313 | 790 |

**關鍵 ratio（不依賴 paired/TO 軸對齊，binary-internal 量化）**：

```
baseline germline-absent  HP1 : HP2 = 3,312 : 791 = 4.19 : 1   ← priority bug 偏移
V3F      germline-absent  hp=33 全部                            ← 保守不選邊
V5       germline-absent  HP1 : HP2 = 3,313 : 790 = 4.19 : 1   ← 與 baseline 一致
```

### 3.2 三版本行為對照

| 版本 | germline-absent events 處理 | 評語 |
|---|---|---|
| **baseline** | 4.19:1 偏 HP1 系列（priority bug 偏移）| 預期：vector ordered check + break early 在 1 票 somatic 觸發 |
| **V3F** | 全部標 hp=33（純 somatic ambiguous，方向不選邊）| **保守正確**：避免錯標方向 |
| **V5** | 4.19:1 偏 HP1（與 baseline 完全相同）| **Layer 1.5 設計繼承 priority bug**！|

### 3.3 對比 chr19 全 events 的偏移

| 區段 | baseline HP1:HP2 |
|---|---|
| chr19 全 events（含 germline）| 1.275:1（接近隨機，因 germline 平衡）|
| chr19 germline-absent only | **4.19:1 偏 HP1** |
| TO mode 全基因組 baseline（全 events 整體）| 17.3:1 偏 HP1（PI §3.3.2）|

→ **germline-absent 區域是 priority bug 偏移的次峰**（介於整體 17.3:1 與 germline-balanced 1.275:1 之間）。

### 3.4 V5 Layer 1.5 的 mechanism 詮釋

V5 `getVote()` Layer 1.5 設計（`d0bcd8c` bundled）：
```cpp
// Layer 1.5: germline 缺席時用 somatic phased votes 決方向
else if (somaticHP1 > 0 || somaticHP2 > 0) {
    germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;
}
```

但 self-phasing 機制下：sub-clone somatic mutations 100% 共現 → 進 phasing graph 後 somatic edges 偏向同一 haplotype → `somaticHP1` vs `somaticHP2` 票數偏向同邊 → Layer 1.5 結果偏向 HP1 系列。

→ **V5 Layer 1.5 在 germline-absent 區域 = priority bug 的 feature 化**：把「baseline 用 somatic vote 蓋過 germline」的 buggy 行為，改成「germline 缺席時才用 somatic vote」的 designed 行為。但**該區域的 4.19:1 偏移本質沒變**。

### 3.5 V3F 為何更正確？

V3F 的設計：germline 缺席時不選邊，標 hp=33（somatic ambiguous）。

```cpp
// V3F (41ff147) two-layer
int germlineResult = 0;  // 0 = undetermined
if (germlineHP1 > 0 || germlineHP2 > 0) {
    germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
} else {
    min = 0; max = 0;  // ← germline 缺席則 germlineResult=0
}
// Layer 2 annotation
if (somaticTotal > 0) {
    hpResult = (germlineResult == 1) ? 11 :
               (germlineResult == 2) ? 21 : 33;  // ← germlineResult=0 → hp=33
}
```

**V3F 的 hp=33** = 「該 read 有 somatic 但 germline 不確定方向」，**避免做不可信的方向判斷**。

→ V3F 正確識別「germline 缺席區域不該強行給方向」，而 V5 試圖補方向反而繼承偏移。

## 4. Caveat 與限制

### 4.1 Paired vs TO HP1/HP2 軸對齊

`paired BAM HP:Z:1` 與 `TO BAM HP:i:1` 在同 chr 內**不保證同 axis**：
- longphase-s 與 longphase-to 各自獨立 phasing
- phase block 內方向一致，跨 phase block 50% 機率 swap
- 所以 cross-tab 中「paired HP:Z:2-1 → TO baseline hp=11」**不能直接解讀為「baseline 把 HP2 reads 錯標到 HP1」**

### 4.2 Binary-internal 量化不依賴軸對齊

但**baseline 自身內部 HP1:HP2 比例**（3,312 vs 791）是 binary-internal 量化，**不依賴 paired 軸**。它告訴我們：

- 在 baseline 自己定義的 HP1/HP2 軸上
- germline-absent events 內
- baseline 標到 HP1 系列（hp=11）的比例是 HP2 系列（hp=21）的 **4.19 倍**
- 如果 sub-clone 結構是隨機（無 priority bug），期望 1:1
- → **4.19:1 是 baseline 自己的 systematic bias 證據**，不是 cross-binary artifact

### 4.3 Sample size

- 5,789 events / 295K paired reads / chr19 only
- 全基因組 germline-absent 區域更大（T1.2-F1 顯示 21.7M reads germline_vote=0）
- 全基因組 cross-ref 估計 ~150K events 量級

## 5. 對 5/8 整合報告 §8.4 的影響

5/8 整合報告 §8.4 把 V5 vs V3F 的 zero-sum 重分配解讀為「Pass 2 reclassify 104K germline het」。本 Step D 揭露**另一層機制**：

| 區段 | V5 vs V3F 行為 |
|---|---|
| germline_vote > 0 區（主流）| Pass 2 reclassify 部分位點為 somatic → reads shift Layer 1 → Layer 1.5（§8.4 量化）|
| **germline_vote = 0 區** | **V3F 標 hp=33 / V5 用 somatic 投票決方向 → 繼承 priority bug 偏移**（本 Step D 量化）|

→ V5 vs V3F 在 germline-absent 區域的差異不只是「tag 數量」，更是**方向偏移 vs 保守 untagged**的設計差異。

## 6. 對 PI 報告 errata E4 的補強

5/9 PI errata E4 已修訂「V5 數值主要功勞 V3F + Layer 1.5（tagging layer）」。本 Step D 補強：

> 在 germline-absent 區域（Layer 1.5 觸發區），V5 行為與 baseline 4.19:1 偏 HP1 完全相同 — Layer 1.5 設計繼承 priority bug 偏移而非修補。V3F 標 hp=33 保守處理在此區域反而**比 V5 更穩健**（避免錯標方向）。

## 7. 後續 follow-up

- **F-paired-D1** 全基因組擴展同 cross-ref（chr19 → 全 chr，估 ~150K germline-absent events）— 確認 4.19:1 偏移是否跨 chr 一致
- **F-paired-D2** Phase block 內 axis-aligned 分析（用 PS tag 對齊 paired vs TO 軸）→ 計算「真正錯標」率
- **F-paired-D3** 評估 V5 Layer 1.5 改用「germline 缺席就標 hp=33（同 V3F）」的影響 — 是否該回歸 V3F 設計？

## 8. 結論

| Q | A |
|---|---|
| paired 同位點 HP1-1+HP2-1 共現但無 HP1/HP2 reads？ | **是** — 5,789 events |
| TO baseline 是否在這些位點偏向 HP1 系列？ | **是** — 4.19:1 偏 HP1（priority bug 在 germline-absent 區域的次峰）|
| V3F 是否修對？ | **是，且更正確** — 全標 hp=33（保守不選邊）|
| V5 Layer 1.5 是否修對？ | **否** — 與 baseline 完全相同 4.19:1 偏移；Layer 1.5 是 priority bug 的 feature 化 |

→ V5 設計可能需要**部分回歸 V3F**：germline-absent 區域改標 hp=33 而非用 somatic 投票決方向。這對 ISM 下游影響的量化是後續 follow-up。
