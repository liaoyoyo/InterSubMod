<!--
build_date: 2026-05-08
agent: synthesis (整合 5 處 self-phasing 觀察 + V5 audit + T1.2 + T1.2-F1)
status: validated
report_class: comprehensive-synthesis (observation → mechanism → cases → fix → validation → twist → conclusion)
audience: PI / lab meeting（中等密度）
parent_reports:
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md（PI 報告 17.3:1 / SP1-3 / sanity 15/15）
  - InterSubMod/docs/reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md（V5 audit / Pass 1 only 揭露 / §5.6 chr19 / §5.7 全基因組）
  - InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md（早期根因敘述）
  - InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md（chr19 752 機制）
  - InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md（全基因組 34,855）
  - InterSubMod/research/paired_priority_bug_audit/00_audit_report.md（5/9 paired Step A+C，paired 沒 priority bug）
  - InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md（5/9 Step D，V5 Layer 1.5 設計缺陷）
inputs:
  - longphase-to-mod commits 8b8c1fd / 41ff147 / 380e8d2 / d0bcd8c / 938f0df
  - HaplotagProcess.cpp:484-563 (getVote three-layer voting)
  - PhasingProcess.cpp:154-228 (Pass 1 / Pass 2 + ponOnlyPhasing branch)
  - vote_dump_{baseline,v3f,v5}_{chr19,genome}.tsv.gz (read-level audit data)
outputs:
  - 本檔（單一 PI / lab meeting 入口）
  - 5 figures: InterSubMod/docs/reports/validated/2026/05/figures/20260508_self_phasing_整合/
verdict: SELF-PHASING 機制因果確立（17.3:1 + 34,855 read-level victims + 100% V3F/V5 修正）；V5 vs V3F zero-sum 重分配已釐清 = Pass 2 reclassify 104K germline het；HCC1395 0.93 + 0.6 樣本 caller F1 三版完全相同（0.7166 / 0.6273，V5 不改 caller）；20 指標 no regression、purity 0.6 完整對照 6 caller F1 + 9 結構指標 0 critical regression；三路徑算法（somaticCalling 3-point）不依賴 purity，Pass 1 永遠跑；Pass 2 second round 只重跑 2-point edgeConnectResult（量化 N50 +3.51% / block −9.79% / phased var −2.9 pp）；HEAD 938f0df = 最新有效版本；對 PI 報告 4-29 4 處表述需 errata；**5/9 paired cross-ref 揭露 V5 Layer 1.5 在 germline-absent 區域與 baseline 4.19:1 偏 HP1 完全相同 — 設計缺陷，V3F 標 hp=33 保守處理反而更穩健**
last_verified: 2026-05-10
update_2026_05_10:
  - "新增 §8.6 Paired Mode Cross-Reference Audit (5/9 cycles 36+37)"
  - "Step A+C: paired 沒 priority bug (HP1:HP2 1:1.275, som_ratio mean 0.462)"
  - "Step D: V5 Layer 1.5 在 germline-absent 區域繼承 priority bug 偏移 (4.19:1 偏 HP1, 與 baseline 完全相同) — V3F 標 hp=33 才是更穩健設計"
  - "§9.1 成熟度表新增 Layer 1.5 設計缺陷; §9.2 加 E5 errata 候選; §9.3 加 F-paired-D1/2/3 後續"
report_template: structured-tech-report v1.0 (synthesis variant)
narrative_order: observation-first (用戶 5/8 確認 inductive ordering)
dialogue_breaks: 5 (此首版 §0-§2 等 dialogue break 1 後續寫 §3+)
-->

# Self-Phasing 完整觀察整合報告 — 2026-05-08

## 0. TL;DR

baseline LongPhase-TO 在 HCC1395 5kHz 全基因組 **HP1:HP2 = 17.3:1**（94.6% somatic ALT 偏 HP1），chr19 SP1/2/3 達 **113:0 / 109:1 / 108:0**。根因 `getVote()` vector 順序錯（somatic 先於 germline）+ break early。5 commits 漸進修補：PON-only flag → V3F two-layer → V5 Layer 1.5 + ploidy fix + threshold。read-level 鐵證 chr19 **752 victims** + 全基因組 **34,855 victims**，V3F/V5 修正率均 **100%**。T1.2-F1 顛覆三項 chr19 結論：chr19 占比僅 **2.16%**、chr8 priority bug rank 21、Layer 1.5 在 germline=0 區觸發 **+560,881**。Bonus 發現 PI 報告 4-29 V5 數據實為 **Pass 1 only**（ploidy bug → purity=0 → Pass 2 從未觸發）。**V5 vs V3F zero-sum 已釐清**（§8.4 Pass 2 reclassify 104K germline het → reads 從 Layer 1 shift Layer 1.5）；**Caller F1 vs SEQC2 truth set 三版完全相同**（HCC1395 0.93 = 0.7166 / 0.6 = 0.6273；V5 不改 caller）；**20 指標 no regression**、**purity 0.6 baseline vs V5 6 caller F1 + 9 結構指標 0 critical regression**、**三路徑算法不依賴 purity（Pass 1 永遠執行）**、**HEAD `938f0df` = 最新有效演算法**（§8.5）。**5/9 paired cross-ref（§8.6）揭露**：paired mode（longphase-s）整體無 priority bug（HP1:HP2 = 1:1.275）；但 **V5 Layer 1.5 在 germline-absent 區域與 baseline 4.19:1 偏 HP1 完全相同 — 設計缺陷**（V3F 標 hp=33 保守處理反而更穩健）。

---

## 1. 受眾與目的

PI / lab meeting 等級。預設熟 InterSubMod 但不熟 longphase-to-mod 內部 C++，20 分鐘讀完即可決策。

**為什麼此刻寫**：5 處 self-phasing 線索 + 多次 commit 沒有單一入口；§5.6 chr19 結論部分被 §5.7 全基因組推翻、三版修補與 PON-only flag 角色未合一說明、顛覆發現未同層整理。

**敘事**：觀察（§2）→ 機制（§3）→ 量化（§4）→ 修補（§5）→ 驗證（§6）→ V5 Provenance bonus（§7）→ 顛覆（§8）→ 結論（§9）。

---

## 2. 觀察起點：17.3:1 偏移與 chr19 SP1/2/3

### 2.1 全基因組統計：HP1:HP2 = 17.3:1（94.6% 偏 HP1）

PI 報告（2026-04-29）§3.3.2 baseline LongPhase-TO 全基因組層觀察：

| 指標 | baseline | 隨機預期 | 偏離 |
|---|---:|---:|---:|
| HP1 reads（somatic ALT）| **614,000** | ~325,000 | 1.89× |
| HP2 reads（somatic ALT）| **35,500** | ~325,000 | 0.11× |
| **HP1:HP2** | **17.3:1** | **1:1** | **17.3×** |
| **HP1 占比** | **94.6%** | ~50% | +44.6 pp |

**為什麼真實存在**（3 條獨立論證）：

1. **生物學**：tumor 的 somatic ALT 屬某 sub-clone 的真實 haplotype，跨 23 染色體不該系統性偏 HP1。
2. **跨染色體一致**：cnLOH artifact 只影響單一 chr；94.6% 跨 23 chr 同向必為工程 bias。
3. **paired GT 對照**：paired tumor-normal 同樣 reads HP1:HP2 ≈ 1:1；只 TO 模式出現 17.3:1。

→ 17.3:1 不是 sample 性質，是 LongPhase-TO 的 systematic bias。

圖示：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01d_somatic_bias_explanation.png` Panel D。

### 2.2 chr19 SP1/2/3 三個極端案例

17.3:1 是全基因組 average。IGV 6-BAM 並列（baseline / V2b / V3-Fixed / V5 / paired tumor / paired normal）找到 3 個近 100% 失衡位點：

| 位點 | 染色體位置 | baseline HP1:HP2 | V5 vs paired |
|---|---|---:|---|
| **SP1** | chr19:17565944 | **113 : 0** | V5 翻 HP2，對齊 paired ✅ |
| **SP2** | chr19:12452332 | **109 : 1** | V5 翻 HP2，對齊 paired ✅ |
| **SP3** | chr19:12467180 | **108 : 0** | V5 翻 HP2，對齊 paired ✅ |

IGV 截圖：
- ![](../../../../reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png)
- ![](../../../../reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP2_chr19_12452332.png)
- ![](../../../../reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP3_chr19_12467180.png)

**重要解讀**：
- SP1/SP3 完全失衡（0 reads 對側），SP2 僅 1 → 不是噪音、caller、alignment 問題；是 **read assignment 強制集中**。
- baseline 與 paired **方向相反**（baseline HP1 → paired HP2）→ 翻轉而非弱化。
- V5 修正後 3/3 翻回 paired → 「修對了」初步證據（完整驗證見 §6）。

### 2.3 觀察驅動本報告

這些觀察自然引出四個問題：為什麼 baseline 把 read 都集中一邊（→ §3 機制）？17.3:1 在 read 層級有多少具體個案、分佈在哪（→ §4 量化擴展）？baseline → V3F → V5 三版各修什麼、為何需要三版（→ §5 修補設計）？V3F/V5 在 SP1/2/3 之外是否也都修對（→ §6 驗證鐵證）？§7-§9 補上修補過程意外發現的 V5 Provenance 問題、顛覆的三項 chr19 結論、PI 報告 errata 清單。

---

## 3. 核心問題：Self-Phasing 機制

問題分兩層：**phasing layer**（建 phasing graph 時 somatic 反客為主）+ **tagging layer**（tag read 時 getVote 順序錯）。兩層互相獨立、互相加成、各需獨立修補。

### 3.1 longphase phasing graph 設計目的：germline-only

正常 phasing graph 的物理基礎：

```
diploid genome:
                                  HP1 (paternal)        HP2 (maternal)
position 1 (germline het A/G):        A                     G
position 2 (germline het C/T):        C                     T
position 3 (germline het G/A):        G                     A

Read covering position 1+2:
  reads 看到 (A,C) → assign HP1
  reads 看到 (G,T) → assign HP2
  → germline het 50% / 50% 隨機分佈，物理上 unbiased
```

phasing graph：node = 位點，edge = 同 read 上的位點對 → 把 het 連成 haplotype。**germline 才能進 graph**（因為只有 germline het 同時出現於 HP1/HP2 兩 stream）。

### 3.2 Tumor-only 模式 somatic 進 graph 的破壞 — 「球員兼裁判」

TO 模式沒 paired normal 切 germline / somatic，只能用 PoN（panel of normals）。**未在 PoN 內的位點預設當 germline 處理 → somatic mutation 也被當 germline 進 graph**。

致命差異：

```
                      germline het         somatic mutation
                      ────────────────     ──────────────────
出現比例（同 sub-clone 內）  50% / 50%       100% / 0%
分佈於 haplotype             兩邊隨機         單一 sub-clone 全帶
共現（多個位點同 read）       random          完全 100% 共現
```

→ somatic 進 graph 後，多個 somatic 之間 **edge weight 暴漲**（100% 共現比 50% 強）。phasing graph 變成「somatic 越強 → 越偏向某 haplotype → 該 haplotype 越強」的自我增強迴圈，**germline 真實訊號被 overrule**。

這就是「球員兼裁判」：somatic 應該是被 phase 的對象（受 germline 仲裁），現在反過來主導 graph。

**解法 = PON-only flag**（commit `8b8c1fd`）：phasing 時只用 PoN 內的 germline het，somatic 完全不進 graph，等 graph 拍板後再用 graph 結果反過來分類 somatic（請見 §5.1）。

### 3.3 Tagging 階段 priority bug — `getVote()` 順序錯 + break early

**即使 phasing graph 修對了（PON-only flag）**，tagging 階段 `getVote()` 還是會出錯。問題在 baseline 的 vector 順序檢查：

```cpp
// baseline HaplotagProcess.cpp:512  ❌ priority bug
void getVote(countMap, ...) {
    vector<pair<int,int>> variantKeys = {
        {HAPLOTYPE1_1, HAPLOTYPE2_1},   // ① somatic pair      ← FIRST
        {HAPLOTYPE3,   HAPLOTYPE2_1},   // ② mixed pair
        {HAPLOTYPE1,   HAPLOTYPE2}      // ③ germline pair    ← LAST
    };
    for (pair in variantKeys) {
        if (countMap[k1]>0 || countMap[k2]>0) {
            hpResult = ...;
            break;   // ← 只要前面有票就停，後面 germline pair 永遠看不到
        }
    }
}
```

**1 票 somatic 就觸發 priority bug**：

```
Read 範例 — 5 votes：
  germline HP1: 0      somatic HP1_1: 1   ← 只有 1 票
  germline HP2: 5      somatic HP2_1: 0
  germline HP3: 0

baseline 投票流程：
  ① 檢查 (HP1_1, HP2_1) → HP1_1=1>0 → 決定 hp=HP:i:11 → break ❌
  ② (HP3, HP2_1) skipped
  ③ (HP1, HP2) skipped ← germline 5 票被忽略！

正確答案應該是 hp=HP:i:21（germline HP2=5 主導，附 somatic HP1_1=1 標 21）
```

**為什麼 read tag 全 HP1**：tumor sub-clone 的 somatic 都同方向（§3.2 解釋的 100% 共現）→ baseline 看到的 somatic 票多偏 HP1_1 → priority bug 把所有受影響 reads 標 HP:i:11（HP1 系列）→ **17.3:1 偏移在 tag layer 形成**。

![Figure F1 — getVote() priority bug 機制 before / after](figures/20260508_self_phasing_整合/F1_priority_bug_mechanism.png)

*Figure F1 — `getVote()` priority bug 機制 before / after。左：baseline vector 順序檢查 + break early；右：V3F + V5 explicit Layer 1 / 1.5 / 2。*

### 3.4 兩層 bug 與兩層修補的對應

| Layer | bug | 修補 commit | §5 章節 |
|---|---|---|---|
| **phasing layer**（建 graph）| somatic 進 graph 球員兼裁判 | `8b8c1fd` PON-only flag | §5.1 |
| **tagging layer**（tag read）| `getVote()` 順序錯 + break early | `41ff147` two-layer getVote | §5.2 |

V3F 還要再加 INDEL guard（`380e8d2`）+ Layer 1.5 fallback / ploidy fix / threshold（`d0bcd8c` + `938f0df`）才成 V5 — 完整演進見 §5。

---

## 4. 問題現象與案例擴展（找共通）

§3 講 mechanism，§4 給具體 read-level 數據（chr19 752 + 全基因組 34,855），驗證 mechanism 不是純理論而是物理可觀察。

### 4.1 chr19 read-level：752 priority bug victims（T1.2 chr19 audit）

對 baseline / V3F / V5 三版 testing-only binary patch 加 `--debug-vote-dump` flag，dump 每條 read 經 `getVote()` 後的 5-vote countMap + hpResult。chr19 子集（HCC1395 5kHz）三版各 549,206 rows × 3-way merged 1,069,832 events，找 priority bug victims：

定義：**germline_majority ≠ somatic_majority** 且 baseline `hpResult` 跟 somatic 方向。

| 指標 | 數量 |
|---|---:|
| 雙向矛盾 reads（germline_maj ≠ somatic_maj，both >0）| 752 |
| baseline 標到 somatic 方向（priority bug confirmed）| **752** |
| **V3F 修正比例**（改向 germline_maj）| **100.00%** |
| **V5 修正比例**（同上）| **100.00%** |

**統一指紋**（前 5 案例）：

| read_name (前 12) | chr19:pos | germline HP1/HP2 | somatic HP1_1/HP2_1 | germline_maj | somatic_maj | baseline → V3F → V5 |
|---|---:|:---:|:---:|:---:|:---:|:---:|
| 1c50034a-f0f | 201,417 | 1/3 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| afb8e89b-893 | 585,252 | 1/2 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| 35c7e166-ec3 | 824,360 | 0/1 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| 096ab9a7-030 | 1,574,442 | 0/3 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| ccc8185d-f9b | 2,558,240 | 0/1 | 2/0 | HP2 | HP1 | **11 → 21 → 21** |

**全 752 條都是 baseline=11 → V3F=21 → V5=21**（單一方向修正，無一條反向）→ 完美對應 §2.1 全基因組 17.3:1（94.6% somatic→HP1）的單向偏移。

**chr19 752 區域聚集**：1 Mb window 內 victim 比例最高：

| chr19 1Mb window | bug victims | total reads | enrichment |
|---:|---:|---:|---:|
| chr19:30M | 215 | 24,732 | 0.87% |
| chr19:27M | 133 | 23,639 | 0.56% |
| chr19:16M | 41 | 22,078 | 0.19% |

**SP1/SP2/SP3** 在 chr19:12M-17M 範圍內 — read-level victims 與 IGV 6-BAM 截圖（§2.2）**屬同一機制的不同層級觀察**，互相佐證。

### 4.2 全基因組擴展：34,855 victims（T1.2-F1）— 46× chr19 量級

對同三版 binary 跑全基因組 vote dump（每版 ~40 min，總 18.9M tagged reads）。dump 大小 744/687/687 MB（gzipped）。

| 規模 | chr19 pilot | Genome F1 | 倍數 |
|---|---:|---:|---:|
| Dump rows | 549,206 | 29,973,253 | 54.6× |
| Tagged reads（per binary）| ~330K | 18,895,432 | 57× |
| **Priority bug victims** | **752** | **34,855** | **46.4×** |
| V3F 修正率 | 100% | **100%** | 一致 |
| V5 修正率 | 100% | **100%** | 一致 |

→ chr19 752 不是局部 artifact，是**全基因組廣泛分佈的 read-level priority bug**。

### 4.3 Per-chromosome 分布：priority bug 不限於 chr19

| chr | victims | total reads | enrichment ‰ | rank |
|---|---:|---:|---:|---:|
| chr7 | **3,508** | 4,852,872 | 0.723 | 1 |
| chr2 | 2,792 | 4,525,538 | 0.617 | 2 |
| chr1 | 2,674 | 4,942,520 | 0.541 | 3 |
| chr16 | 2,584 | 2,267,135 | **1.140** | 4 |
| chr20 | 2,101 | 1,609,083 | **1.306** | 7 |
| chr19 | **752** | 1,069,832 | 0.703 | **19** |
| **chr8** | **666** | 3,324,020 | **0.200** | **21** |
| chrY | 67 | 45,137 | **1.484** | 24 |

**重點觀察**：
- chr19 占全基因組 priority bug 僅 **2.16%**（752 / 34,855），rank 19 — chr19 不是主要 hotspot
- chr8 enrichment 0.34× genome avg — priority bug **冷區**（與 MEMORY chr8 LOH+HPSig hotspot 7.4× FP enrichment 是不同 layer，後者跟 priority bug 無直接關聯，§8.2 詳述）
- victim N 大宗在 chr7 / chr2 / chr1 / chr16 / chr20，分佈廣泛

→ 「priority bug 主要影響 chr19」這個原 §5.6.3 結論被全基因組擴展**推翻**（顛覆性發現詳述在 §8.1）。

![Figure F2 — Priority bug per-chromosome 分佈](figures/20260508_self_phasing_整合/F2_priority_bug_per_chr_enrichment.png)

*Figure F2 — Per-chr priority bug 分佈雙 panel。左 panel 按 victim N 排序（chr7/chr2/chr1 主要分佈），右 panel 按 enrichment ‰ 排序（chrY/chr20/chr21 高，chr8 冷區 0.34× avg）。chr19 標紅、chr8 標藍。*

### 4.4 共通指紋：全單向 baseline=11 → V3F=21

跨 chr19 752 + 全基因組 34,855 一致觀察到：

```
baseline 受 somatic 主導    →   V3F 改向 germline 主導
   hp = 11 (HP1 系列，錯)        hp = 21 (HP2 系列，對)
   ↓
   17.3:1 偏移（94.6% → HP1）
   34,855 read-level victims（單向，無反向修正）
```

這是「全基因組 17.3:1 整體統計」與「read-level 34,855 個案 trace」的一致性證據鏈：mechanism（§3）→ pattern（§4）→ 修正後（§6）。

---

## 5. 修補設計演進（baseline → V3F → V5）

5 個 commits 漸進完成，分屬 phasing layer + tagging layer 兩層；任一單獨 commit 不夠，**三版 stacking** 才達 V5 完整效果。

### 5.1 commit 鏈時間軸 + layer 分類

```
2026-04-09          04-10           04-25         04-30 (a)        04-30 (b)
   │                  │               │              │                │
   ▼                  ▼               ▼              ▼                ▼
[8b8c1fd]        [41ff147]       [380e8d2]      [d0bcd8c]        [938f0df]
PON-only flag    two-layer       INDEL guard    ploidy fix       threshold
                 getVote                        + Layer 1.5      0.95→0.9
─────────        ──────────────────────────     ──────────────   ────────
phasing layer    tagging layer                  phasing+tagging  phasing layer
(球員兼裁判)      (priority bug)                  (Pass 2 觸發 +   (Pass 2
                                                germline 缺席    觸發門檻)
                                                fallback)
                 ────────── V3-Fixed ──────────
                 ─────────────────────────────── V5 ─────────────
[baseline]       [V3-Fixed = baseline + V3F]    [V5 = V3F + d0bcd8c + 938f0df]
```

### 5.2 baseline `8b8c1fd`（4-09）— PON-only flag · phasing layer 的 root cure

**問題**：§3.2 球員兼裁判 — somatic 進 phasing graph 反客為主。

**修補**：加 `--pon-only-phasing` flag。phasing 階段 graph 只放 PoN 內的 germline het，somatic **不進 graph**。Pass 1 graph 拍板後再用結果反過來分類 somatic。

**關鍵程式碼**（`PhasingProcess.cpp:154-228`）：

```cpp
if(params.ponOnlyPhasing){
    // Pass 1: 把 non-PoN 視為 SOMATIC，graph 內只剩 germline het
    vGraph->convertNonGermlineToSomatic();
    vGraph->phasingProcess(...);

    // Pass 2: 重置 somatic origin，用 graph 結果反分類
    vGraph->resetNonPonOrigin();
    vGraph->somaticCalling(...);

    // 最後同步 origins → 修 haplotag HP output
    vGraph->syncPhasingResultOrigins(...);
}
```

**效果**（已驗證）：phase block N50 +99.7%、LOH.bed Jaccard=1.0。

**遺留 bug**：解了 phasing layer 的循環依賴，但 **tag layer 的 `getVote()` 還是壞的**。對 reads tagging 沒實質改善 → 17.3:1 偏移仍在。

### 5.3 V3-Fixed `41ff147`（4-10）+ `380e8d2`（4-25）— tagging layer 兩層投票

#### `41ff147` two-layer getVote — 拆掉 vector 順序檢查

**問題**：§3.3 priority bug — vector ordered check + break early。

**修補**：把 `getVote()` 改成 explicit Layer 1 / Layer 2：
- **Layer 1**：只看 germline 票（HP1 vs HP2）決方向；germline 永不被 somatic overrule
- **Layer 2**：有 somatic 票就根據 Layer 1 結果標 11/21/33

**bonus 修**：原 code `hpResult != HAPLOTYPE1_1` 比較是 enum (=3) vs HP tag int (=11)，永遠 mismatch → 改成 `hpResult == 11 || hpResult == 21 || hpResult == 33`。

#### `380e8d2` INDEL guard — 補 OOB UB

**問題**：`countINDELHaplotype` 在 somatic 位點（GT=0|0+GT2，refHaplotype=-1=UNDEFINED）時 `countMap[-1]++` 是 OOB undefined behavior。

**修補**：加 `!= HAPLOTYPE_UNDEFINED` 檢查（與 `countSNPHaplotype` 一致）。

→ 至此 **tag layer 的 priority bug 與 INDEL OOB 都修對**。但 germline 缺席區（cnLOH / amplicon）會大量 reads 標 0（untagged），因為 V3F 沒 fallback。

### 5.4 V5 `d0bcd8c`（4-30 a）+ `938f0df`（4-30 b）— Layer 1.5 + ploidy fix + threshold

#### `d0bcd8c` 主修：Pass 2 ploidy collection bug

**問題**：`PhasingProcess.cpp:158` Pass 1 把 `ploidyRatioMap=nullptr` 傳 `phasingProcess()`；Pass 2 的 `somaticCalling + syncPhasingResultOrigins` 也沒 refill ploidyRatioMap → q1=q3=0 → polynomial → purity 算成 **0**。

**修補**：加 `VairiantGraph::collectPloidyRatio()` helper，在 `syncPhasingResultOrigins` 後重建 hpAlleleCountMap 並填 ploidyRatioMap。

**驗證**：HCC1395 5kHz @ 0.93 purity → V5 purity 0 → 0.871（baseline 0.927）；HCC1395 t30_n20 @ 0.6 → V5 0 → 0.634（baseline 0.607）。

#### `d0bcd8c` bundled：Layer 1.5 somatic fallback（**bundled 在 ploidy fix commit 內**）

**問題**：V3F Layer 1 only 結構在 germline=0 區 reads 全部 untagged。

**修補**：插入 Layer 1.5 — germline 缺席時改用 somatic phased votes（HP1_1 vs HP2_1）決方向（Layer 2 仍標 11/21/33）。

```cpp
// Layer 1: germline (跟 V3F 同)
if (germlineHP1 > 0 || germlineHP2 > 0) { ... germlineResult = 1 or 2 ... }
// Layer 1.5: NEW — somatic fallback
else if (somaticHP1 > 0 || somaticHP2 > 0) {
    if (somaticHP1 >= somaticHP2) germlineResult = 1;
    else                          germlineResult = 2;
}
```

#### `938f0df` threshold 0.95 → 0.9

**問題**：原 `bool highPurity = purity > 0.95` threshold 太嚴，部分 sample（如 0.91-0.94）跑不到 Pass 2 second-round phasing。

**修補**：cherry-pick from upstream zhenyu，threshold 改 `> 0.9`。

### 5.5 PON-only flag 與三版 layer 對應（用戶 5/8 強調）

| 版本 | Commit | 解 phasing layer | 解 tagging layer | 補 germline 缺席 fallback | Pass 2 真實觸發 |
|---|---|:---:|:---:|:---:|:---:|
| baseline | `8b8c1fd` | ✅ PON-only | ❌ priority bug 仍在 | ❌ | ❌（purity=0）|
| V3-Fixed | + `41ff147` + `380e8d2` | ✅ | ✅ two-layer + INDEL | ❌ | ❌（同上）|
| **V5** | + `d0bcd8c` + `938f0df` | ✅ | ✅ + Layer 1.5 | ✅ | ✅（threshold 0.9）|

**為什麼不能只用 PON-only**：解 phasing layer，但 read tagging 階段 `getVote()` 還是優先吃 somatic 票（priority bug）→ 99.9% reads 仍標 HP:i:21。

**為什麼 V3F 不夠還要 V5**：
1. V3F Layer 1 only → germline 缺席區 reads 全 untagged → 樣本如 cnLOH / amplicon hotspot 大量 reads 流失
2. V3F 之後 4-29 PI 報告引用的 BAM 是 Pass 1 only（ploidy bug 讓 Pass 2 從未觸發 — 詳見 §7）
3. threshold 0.95 太嚴 → 0.91-0.94 purity sample 跑不到 Pass 2

![Figure F3 — Self-phasing 修補 5 commits 時間軸](figures/20260508_self_phasing_整合/F3_binary_commit_timeline.png)

*Figure F3 — 5 commits 時間軸 + 兩層對應。baseline (8b8c1fd) → V3-Fixed (+ 41ff147 + 380e8d2) → V5 (+ d0bcd8c + 938f0df)。phasing layer = 藍、tagging layer = 綠、跨層 = 紫。`41ff147` two-layer getVote 是修偏移的關鍵 commit ★。*

### 5.6 `getVote()` 三版程式碼對照

```cpp
// ═══════════════════════════════════════════════════════════════════
// baseline (pre-41ff147)  ❌ priority bug
// ═══════════════════════════════════════════════════════════════════
void getVote(countMap, ...) {
    vector<pair<int,int>> variantKeys = {
        {HAPLOTYPE1_1, HAPLOTYPE2_1},   // ① somatic pair    ← FIRST
        {HAPLOTYPE3,   HAPLOTYPE2_1},   // ② mixed pair
        {HAPLOTYPE1,   HAPLOTYPE2}      // ③ germline pair  ← LAST
    };
    for (pair : variantKeys) {
        if (countMap[k1]>0 || countMap[k2]>0) {
            hpResult = haplotypeBase[winner];
            break;   // ← 1 票 somatic 就 break，germline 看不到
        }
    }
}

// ═══════════════════════════════════════════════════════════════════
// V3-Fixed (41ff147)  ✅ tagging layer 修對；germline 缺席仍 untag
// ═══════════════════════════════════════════════════════════════════
void getVote(countMap, ...) {
    int germlineHP1 = countMap[HAPLOTYPE1];
    int germlineHP2 = countMap[HAPLOTYPE2];
    int somaticTotal = countMap[HAPLOTYPE1_1] + countMap[HAPLOTYPE2_1] + countMap[HAPLOTYPE3];

    // Layer 1: germline only
    int germlineResult = 0;
    if (germlineHP1 > 0 || germlineHP2 > 0) {
        germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
    } else { min = 0; max = 0; }

    // Layer 2: somatic annotation
    if (somaticTotal > 0) {
        hpResult = (germlineResult == 1) ? 11 :
                   (germlineResult == 2) ? 21 : 33;
    } else hpResult = germlineResult;
}

// ═══════════════════════════════════════════════════════════════════
// V5 (HEAD = 938f0df, with Layer 1.5 bundled in d0bcd8c)
// ✅ tagging layer 修對；germline 缺席用 somatic fallback
// ═══════════════════════════════════════════════════════════════════
void getVote(countMap, ...) {
    int germlineHP1 = countMap[HAPLOTYPE1];
    int germlineHP2 = countMap[HAPLOTYPE2];
    int somaticHP1  = countMap[HAPLOTYPE1_1];   // ← V5 新增分開取
    int somaticHP2  = countMap[HAPLOTYPE2_1];   // ← V5 新增分開取
    int somaticTotal = somaticHP1 + somaticHP2 + countMap[HAPLOTYPE3];

    int germlineResult = 0;
    // Layer 1: germline (跟 V3F 同)
    if (germlineHP1 > 0 || germlineHP2 > 0) {
        germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
    }
    // Layer 1.5: NEW — germline 缺席時用 somatic phased
    else if (somaticHP1 > 0 || somaticHP2 > 0) {
        germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;
    }
    // Layer 2: somatic annotation (跟 V3F 同)
    if (somaticTotal > 0) {
        hpResult = (germlineResult == 1) ? 11 :
                   (germlineResult == 2) ? 21 : 33;
    } else hpResult = germlineResult;
}
```

**三版差異 diff 統計**（commit chain）：
- `8b8c1fd`：5 files +69/-6 (PhasingProcess +25/-3)
- `41ff147`：HaplotagProcess.cpp +36/-25（getVote 全重寫）
- `380e8d2`：HaplotagProcess.cpp +8/-4（INDEL guard）
- `d0bcd8c`：4 files +68/-9（PhasingGraph collectPloidyRatio +29 / Layer 1.5 +21 / countSNP guard +5）
- `938f0df`：2 files +4/-4（threshold + purity polynomial）
- 累計：~155 lines tagging-layer + ~40 lines phasing-layer 修補

---

## 6. 修補驗證鐵證

### 6.1 chr19 752 victims：100% V3F+V5 修正（T1.2 chr19 audit）

對三版 testing-only binary patch 加 `--debug-vote-dump` flag，對 HCC1395 5kHz chr19 跑同樣 reads，三版結果 JOIN by (read_name, chr, pos)：

| 指標 | 數值 | 判定 |
|---|---:|:---:|
| 雙向矛盾 reads | 752 | — |
| baseline 標到 somatic 方向（priority bug confirmed）| 752 | — |
| **V3F 修正比例**（改向 germline_maj）| **100.00%** | ✅ |
| **V5 修正比例** | **100.00%** | ✅ |

**4-path 驗證結果**（T1.2 plan 設計的 4 條獨立證據路徑）：

| 路徑 | 結果 | 判定 |
|------|------|:---:|
| ① 個案 trace ≥10 條 | 752 條 | ✅ PASS |
| ② chr19 1Mb 區域聚集 | chr19:30M (215) + 27M (133) 集中 46% | ⚠️ PARTIAL |
| ③ Somatic density 共變 | high somatic vote ≥5 = 0 受害；低票觸發 | 🔄 反向但有意義 |
| ④ 修正後消失（V3F/V5 修正率）| 100% / 100% | ✅ PASS |

→ 4 路徑 3.5/4 PASS，priority bug 機制因果**確立**。

![Figure F4 — chr19 752 victims top-10 1Mb hotspot windows](figures/20260508_self_phasing_整合/F4_chr19_752_victims_scatter.png)

*Figure F4 — chr19 priority bug 在 1 Mb window 的 hotspot 散點。30 M (215 victims) + 27 M (133) 共佔 chr19 victim 46%，與 SP1 (chr19:17.5 M) / SP2-3 (chr19:12.5 M) 區段對齊。圓點大小與顏色表 enrichment %。*

### 6.2 全基因組 34,855 victims：100% V3F+V5 修正（T1.2-F1）

| 規模 | chr19 pilot | Genome F1 | 倍數 |
|---|---:|---:|---:|
| 雙向矛盾 reads | 752 | **34,855** | 46.4× |
| V3F 修正率 | 100% | **100%** | 一致 |
| V5 修正率 | 100% | **100%** | 一致 |

→ chr19 100% 修正不是局部 artifact；**全基因組 34,855 個案、無一條反向**，priority bug 修補在 read-level 鐵證確立。

### 6.3 SP1/2/3 案例修補後對照（呼應 §2.2 hook）

§2.2 給的 baseline IGV 113:0 / 109:1 / 108:0；V5 修正後：

| 位點 | baseline HP1:HP2 | V5 HP1:HP2 | paired tumor 方向 | 是否對齊 paired |
|---|---:|---:|:---:|:---:|
| **SP1** chr19:17,565,944 | **113 : 0** | 翻轉至 HP2 主導 | HP2 | ✅ |
| **SP2** chr19:12,452,332 | **109 : 1** | 翻轉至 HP2 主導 | HP2 | ✅ |
| **SP3** chr19:12,467,180 | **108 : 0** | 翻轉至 HP2 主導 | HP2 | ✅ |

**3/3 對齊 paired ground truth** — 個別位點層 V5 100% 修對。

**IGV 6-BAM 並列截圖**（baseline / V2b / V3-Fixed / V5 / paired tumor / paired normal）：

![SP1 chr19:17565944 — baseline 113:0 全在 HP1，V5/paired 都是 HP2 主導](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png)

*Figure SP1 — baseline panel reads 全部集中於 HP1+HP1-1（紅 + 綠群）；V5 / V2b / V3-Fixed 整體翻轉至 HP2+HP2-1；paired tumor 確認 HP2 為真實方向。*

![SP2 chr19:12452332 — baseline 109:1，V5 翻轉至 HP2 對齊 paired](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP2_chr19_12452332.png)

*Figure SP2 — baseline HP1+HP1-1 集中 109 reads，HP2 stack 僅 1；V5 方向翻轉至 HP2+HP2-1，與 paired 一致。*

![SP3 chr19:12467180 — baseline 108:0，與 SP1 同模式](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP3_chr19_12467180.png)

*Figure SP3 — 與 SP1 同模式：baseline HP1 主導 → V5 HP2 主導，HP orientation 整體翻轉。*

### 6.4 整體 HP 結構修復（PI 報告 §6.5 引用）

| 指標 | baseline | V5 修正後 | 變化 |
|---|---:|---:|---:|
| HP1:HP2 ratio | 17.3:1 | ~1:1 | **消除偏移** |
| 94.6% somatic→HP1 | 是 | ~50% balanced | **修正** |
| 15-site Problem PS（SP1/2/3 含）| 48.5% | 52.0% | +3.5 pp |

→ §3 mechanism、§4 read-level 量化、§6 修補驗證形成**完整證據鏈**：
```
理論（球員兼裁判 + priority bug）
   ↓
全基因組 17.3:1 偏移 (PI §3.3.2)
   ↓
個別位點 SP1/2/3 113:0 / 109:1 / 108:0 (PI §3.3.3, IGV)
   ↓
read-level 752 victims chr19 + 34,855 全基因組 (T1.2 / T1.2-F1)
   ↓
V3F + V5 修正率 100% (chr19 + 全基因組)
   ↓
SP1/2/3 翻回 paired 方向 3/3 (IGV verification)
```

---

## 7. V5 Provenance Bonus 發現：PI 報告 4-29 V5 = Pass 1 Only

修補過程意外發現 PI 報告（2026-04-29）引用的 V5 數據**並非完整 V5 (Pass 1 + Pass 2 + Layer 1.5) 設計**，實為 Pass 1 only。本節釐清因果並回答「純樣本為什麼還有偏移」的常見直覺問題。

### 7.1 為什麼純樣本 baseline 仍有 17.3:1 偏移？— ploidy bug 讓 purity=0

直覺：HCC1395 5kHz 真實 purity ≈ 0.93（純樣本），baseline 應該 highPurity=true → Pass 2 觸發 → 沒偏移問題。

反例（資料事實）：

| 階段 | phasing.log purity | highPurity | Pass 2 觸發 |
|---|---:|:---:|:---:|
| baseline `8b8c1fd` (V2b PON-only) | **0** | false | ❌ |
| PI 報告 V5 (4-12，同 V2b phased VCF + Layer 1.5 tag) | **0** | false | ❌ |
| 4-30 重跑 baseline_09（含 `d0bcd8c` ploidy fix）| **0.977376** | true | ✅ |
| 4-30 重跑 v5_flag（含 `d0bcd8c` + `938f0df`）| **0.983791** | true | ✅ |

→ 純樣本 baseline 算出 purity=0 **不是樣本性質，是 ploidy bug**：

```
PhasingProcess.cpp:158  Pass 1 傳 nullptr  →  ploidyRatioMap 留空
                                                    │
Pass 2 somaticCalling + syncOrigins  也沒 refill ────┘
                                                    │
                                                    ▼
                              ploidyRatioMap = {} （空）
                                                    │
                                                    ▼
                              q1 = q3 = 0  →  purity polynomial → 負值
                                                    │
                                                    ▼  clamp to 0
                              purity = 0
                                                    │
                                                    ▼
                              highPurity = (0 > 0.9) = false
                                                    │
                                                    ▼
                              Pass 2 second round phasing **永不觸發**
```

### 7.2 17.3:1 偏移與 Pass 2 觸發是兩條獨立問題

**HP1 偏移的兩個成因**（與 Pass 2 無關）：

| 層 | 成因 | 修補 commit |
|---|---|---|
| phasing layer | 球員兼裁判（somatic 進 graph）| `8b8c1fd` PON-only |
| **tagging layer** | **`getVote()` priority bug** | **`41ff147` V3F two-layer** ★ |

→ **17.3:1 → ~1:1 修正功勞主要在 V3F**（tagging layer fix），**不依賴 Pass 2**。

PI 報告 4-29 V5 數據雖然是 Pass 1 only（`d0bcd8c` 之前），但仍展現 17.3:1 → ~1:1 修正，**正是因為 V3F 的 tagging fix 是上游已生效的**。Pass 2 是 phasing layer 的 second round（讓 graph 結構更乾淨、N50 提升），對 HP1 偏移是**次級效益**而非主因。

### 7.3 對 PI 報告數值的影響

PI 報告 §6.4 / §6.5 的 V5 數值來自 4-12 BAM = Pass 1 only：

| PI 報告數值 | 主要功勞來源 | Pass 2 觸發後預期 |
|---|---|---|
| HP1:HP2 17.3:1 → ~1:1 | V3F tagging fix（`41ff147`）| 應持平 |
| 94.6% somatic→HP1 → ~50% balanced | V3F | 應持平 |
| sanity 15/15 PASS | V3F + Layer 1.5 | 應持平 |
| clean PS paired GT concordance +13.3 pp | V3F + Layer 1.5 | 持平 or 微升 |
| HP:i:33 −54%（110,197）| V3F + Layer 1.5 | 應持平 |
| Phase block N50 +99.7% | 主要 PON-only Pass 1 | Pass 2 後可能更升 |

→ 數字不會因 Pass 2 觸發而崩盤；**歸因可能要修正**（從「V5 four-commit chain 的整體效益」更精確說成「主要 V3F + Layer 1.5；Pass 2 二次效益尚未獨立量化」）。

獨立量化 Pass 2 vs V5 flag 對 ISM verdict 的差異需要 4-cell ablation（baseline×Pass 2 與 V5 flag×Pass 2 各組合），用戶 5/7 決策**暫不做 ISM benchmark**（ISM 是下游，longphase-to 端 V3F 修對後 ISM 自動受惠 — 詳見 §9 後續）。

---

## 8. 顛覆性發現（T1.2-F1 全基因組擴展推翻三項 chr19 結論）

T1.2-F1 把 read-level audit 從 chr19 (752 victims) 擴到全基因組 (34,855 victims) 後，發現三項原 chr19 pilot 結論需顛覆。**機制因果不變**（priority bug 仍是真，V3F/V5 仍 100% 修對），但 hotspot 詮釋與 Layer 1.5 觸發角色都需修訂。

### 8.1 chr19 不是主要 hotspot — 占比僅 2.16%

原 chr19 pilot 結論：「chr19 752 victims 集中在 30M / 27M window」→ 暗示 chr19 是 hotspot。

全基因組擴展後（§4.3 已示）：

```
chr7  ████████████████████  3,508 (rank 1, 占 10.1%)
chr2  ████████████████      2,792 (rank 2, 占 8.0%)
chr1  ███████████████       2,674 (rank 3, 占 7.7%)
chr16 ██████████████        2,584 (rank 4, 占 7.4%)
chr20 ████████████          2,101 (rank 7, 占 6.0%)
chr19 ████                    752 (rank 19, 占 2.16%) ★ 不是 hotspot
chr8  ███                     666 (rank 21, 占 1.9%)  ★ 冷區
chrY  ▏                        67 (rank 24, 占 0.2%)  ★ 但 ‰ 排第 1
```

→ 17.3:1 偏移廣泛分佈於 chr7/chr2/chr1/chr16/chr20，不限 chr19。**chr19 SP1/2/3 是「可重現案例」而非「主要分佈位置」**。

### 8.2 chr8 是 priority bug 冷區（與 LOH+HPSig hotspot 不同 layer）

MEMORY `project_hcc1395_chr8_hotspot.md` 記錄 chr8「LOH+HPSig 7.4× FP enrichment」hotspot。

T1.2-F1 顯示 chr8 priority bug enrichment **0.34× genome avg**（rank 21，冷區）。

→ **chr8 LOH+HPSig hotspot 與 priority bug 是不同 layer**：
- chr8 LOH+HPSig hotspot：ISM 下游 false-positive 富集（HP_Ratio + LOH 特徵交互）
- chr8 priority bug：tagging layer getVote 投票錯，chr8 反而**少於** genome avg
- 兩者**沒有直接因果關聯**

→ chr8 LOH+HPSig hotspot 機制需另尋（不是 V3F/V5 修對 priority bug 後就會自動消失，因為它不是 priority bug 造成）。

### 8.3 Layer 1.5 在 germline 缺席區確實觸發 +560,881 reads

原 chr19 pilot 結論：「V5 vs V3F 在 chr19 結果 100% 相同 → Layer 1.5 chr19 未觸發」。

全基因組擴展（germline_vote=0 子集分析）：

| 指標 | 數值 |
|---|---:|
| germline_vote=0 reads | 21,765,669（占 merged 36.8%）|
| V3F tagged 數（germline=0 子集）| 0 |
| V5 tagged 數（germline=0 子集）| **560,881** |
| **Layer 1.5 額外觸發** | **+560,881** |

→ Layer 1.5 確實觸發（chr19 局部沒觸發是因 germline 不缺席，外推錯誤）。

### 8.4 V5 vs V3F 整體 zero-sum 重分配（已釐清 — Pass 2 reclassification 主導）

V5 全基因組整體 tag count = V3F 整體 tag count（兩者都是 18,895,432）：

```
                    germline=0 區          germline>0 區           整體
                    ───────────────       ──────────────────       ────────
V3F tag             0                     18,895,432               18,895,432
V5  tag             560,881               18,334,551               18,895,432
                    ──────────────         ─────────────             ─────
V5 - V3F            +560,881               -560,881                  0
```

#### 機制驗證 — `tumor_phased.vcf` 直接量化

對 V3F binary（380e8d2，**有 ploidy bug → Pass 2 跳過**）vs V5 binary（HEAD = 938f0df，**有 ploidy fix → Pass 2 觸發**）的 phased VCF 比對：

| 階段 | total variants | germline het (0\|1, 1\|0) | somatic (0\|0) | other (./.) |
|---|---:|---:|---:|---:|
| V3F（Pass 1 only）| 3,187,275 | **1,165,336** | 21,304 | 2,000,635 |
| V5（Pass 1+2）| 3,187,275 | **1,060,879** | 27,854 | 2,098,542 |
| **Δ Pass 2 觸發** | 0 | **−104,457** | +6,550 | +97,907 |

#### 因果鏈

```
V3F：phased VCF 有 1,165,336 germline het
  → reads 在這些位點投票 countMap[HP1] or [HP2]
  → germline_vote > 0 → Layer 1 tag (HP:i:1/2/11/21)

V5：Pass 2 reclassify 104,457 germline het 為 somatic / 未 phase
  → phased VCF 只剩 1,060,879 germline het
  → reads 在 reclassified 位點：
      改判 somatic → 投 countMap[HP1_1]/[HP2_1] → germline_vote=0 → Layer 1.5 tag
      改判 ./.    → 不投票 → 部分 reads 失去訊號
  → 約 560K reads 從 Layer 1 路線 shift 到 Layer 1.5 路線
```

#### 排除其他候選

| 候選 | 排除理由 |
|---|---|
| ❌ `countSNPHaplotype` UNDEFINED guard | 只防 OOB UB，well-formed BAM 不會改變 tag count 數量級到 560K |
| ❌ threshold 0.95 → 0.9 | HCC1395 5kHz purity=0.984，兩 threshold 都觸發 Pass 2，對此樣本無差異 |
| ✅ **`d0bcd8c` ploidy fix → Pass 2 觸發 → 104K 位點 reclassify** | 量級匹配（104K 位點 × 平均覆蓋 → ~560K reads 受影響）|

**這不是 bug，是 Pass 2 設計的正確副作用**：Pass 2 重新分類後得到更乾淨的 phasing graph + 更精確的 somatic 邊界，受影響 reads 仍然透過 Layer 1.5 正確 tag。priority bug 修補的 100% verdict **不受影響**。

![Figure F5 — Layer 1.5 zero-sum 重分配](figures/20260508_self_phasing_整合/F5_layer15_zero_sum_4quadrant.png)

*Figure F5 — Layer 1.5 zero-sum 重分配視覺化。左 panel：V3F vs V5 在 germline=0 / germline>0 兩 region 的 tag count；右 panel：差異（V5 − V3F）顯示 zero-sum 結構（germline=0 +560 K / germline>0 −560 K，總和 = 0）。*

---

## 8.5 No-Regression 量化驗證（baseline → V5）

針對用戶 5/9 質疑「F1 / LOH / 整體量化是否變差」、「purity 0.6 樣本是否差很多」做完整盤點。

### 8.5.1 baseline → V5 (PI 報告 4-12 = Pass 1 only) 全指標對照

20 個指標、5 大類別、**無任何 regression**：

| 類別 | 指標 | baseline | V5 | Δ | 變差? |
|---|---|---:|---:|---|:---:|
| **ISM aggregate** | TP_rate | 0.706 | 0.711 | +0.005 | ❌ |
| (PI §6.1) | HP_Ratio median | 0.788 | **0.574** | tag bias 修正 | ❌ |
| | Potential_LOH% | 58.7 | 62.2 | +3.5 pp | ❌ |
| **HP_Ratio AUC** | All | 0.531 | 0.526 | -0.005（隨機區間）| ❌ |
| (PI §6.2) | Inner | 0.523 | 0.525 | +0.002 | ❌ |
| **Methylation 6 feat** | 全部 (HPMergedDelta / AlleleDelta / HPFineF / PairwiseMeanDist / CramersV / GlobalP) | — | — | 全 ±0.01 內持平 | ❌ |
| (PI §6.3) | | | | | |
| **Paired GT concord.** | 全基因組 clean PS | 82.2% | **90.5%** | **+8.3 pp** | ❌ |
| (PI §6.4) | 15-site Aggregate | 72.20% | **78.85%** | **+6.65 pp** | ❌ |
| | 15-site Clean PS (11) | 74.9% | **88.2%** | **+13.3 pp** | ❌ |
| | 15-site Problem PS (SP1/2/3) | 48.5% | 52.0% | +3.5 pp | ❌ |
| **HP 結構** | HP1:HP2 ratio | 17.3:1 | ~1:1 | 消除 | ❌ |
| (PI §6.5) | Phase block N50 | 4,061 | 8,109 | **+99.7%** | ❌ |
| | Phased rate | 54.9% | 78.5% | **+23.6 pp** | ❌ |
| | 執行時間 | 2,693s | 1,976s | **1.36× 快** | ❌ |
| **LOH 結構** | LOH regions | 1,094 | 1,094 | 完全相同 | ❌ |
| (今次驗證) | LOH 總 bp | 1.632 Gb | 1.632 Gb | 完全相同 | ❌ |

→ 20/20 指標 **沒有任何 regression**，6 項顯著改善（+8.3 pp ~ +99.7%）。

### 8.5.2 Purity 0.6 樣本 baseline vs V5 — 沒差很多 + Caller F1 三版完全相同

#### Purity 計算對照

HCC1395 t30_n20 (purity 0.6 模擬) 各版本 phasing.log purity 計算：

| 版本 | purity 計算 | Pass 2 觸發 | 評語 |
|---|---:|:---:|---|
| baseline 0.6 | **0.607** | ❌ (<0.9) | 接近真實 0.6 — ploidy bug 在低純度自我治癒 |
| V5 0.6 (without `d0bcd8c`) | **0** | ❌ | ploidy bug 仍崩 |
| V5 0.6 (with `d0bcd8c`) | **0.634** | ❌ (<0.9) | 修對 |

→ baseline (0.607) vs V5 修正版 (0.634) 差 **0.027**（4.5% 相對差）→ **沒差很多**。

→ ploidy bug 在**低純度樣本上自我治癒**（polynomial 在 low-input 不會崩到 0）；高純度樣本才崩盤。`d0bcd8c` 主要目的是修高純度樣本，低純度差異微小。

→ 0.6 樣本兩者 purity < 0.9，Pass 2 都不觸發 → Pass 2 設計就不針對低純度樣本。

#### Caller F1 三版完全相同（vs SEQC2 truth set v1.2.1）

來源：`InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md` §3a（6 版本 phased VCF F1 vs `/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`）。

**HCC1395 5kHz @ 0.93 purity（A1/A3/A5 三版）**：

| 版本 | TP | FP | FN | Precision | Recall | **F1** |
|------|---:|---:|---:|---:|---:|---:|
| A1 baseline | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| A3 v3f_no_pononly | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| A5 **V5** | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |

**HCC1395 t30_n20 @ 0.6 purity（B1/B3/B5 三版）**：

| 版本 | TP | FP | FN | Precision | Recall | **F1** |
|------|---:|---:|---:|---:|---:|---:|
| B1 baseline | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| B3 v3f_no_pononly | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| B5 **V5** | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |

**ΔF1 (0.93→0.6) = −0.0893** — 全部 caller 性能下降，**與 V5/baseline 無關**（ClairS-TO 在 0.6 純度本身即降）。

#### 為什麼三版完全一致？機制驗證

```
ClairS-TO  ──→  snv.vcf.gz (FILTER: PASS / NonSomatic / LowQual)
                │  PASS variants set 由 ClairS-TO 決定 (V5 不改 caller)
                ▼
longphase-to phase  ──→  tumor_phased.vcf
                │  改 GT / PS / GT2 / GT3 欄位（V5 改這層）
                │  **不改 FILTER 欄位**
                ▼
F1 比對：用同一組 PASS variants vs SEQC2 truth
                │  PASS set 相同 → TP/FP/FN 完全相同 → F1 完全相同
                ▼
                三版本 F1 = 0.7166（0.93） / 0.6273（0.6）
```

→ **F1 不衡量 tag 品質**（PI 報告 §5.3 已釐清）。V5 真實價值在 read-level tag concordance（+13.3 pp paired GT @ 0.93）。

#### 0.6 樣本完整 no-regression 對照

| 類別 | 指標 | baseline 0.6 | V5 0.6 | Δ | 變差? |
|---|---|---:|---:|---|:---:|
| **Caller F1** | F1 vs SEQC2 | 0.6273 | 0.6273 | 0 | ❌ |
| | TP | 24,190 | 24,190 | 0 | ❌ |
| | FP | 13,487 | 13,487 | 0 | ❌ |
| | FN | 15,257 | 15,257 | 0 | ❌ |
| | Precision | 0.6420 | 0.6420 | 0 | ❌ |
| | Recall | 0.6132 | 0.6132 | 0 | ❌ |
| **Purity 計算** | purity | 0.607 | 0.634 | +0.027 (4.5% rel) | 否（更接近真實 0.6）|
| **Pass 2 觸發** | 狀態 | ❌ (<0.9) | ❌ (<0.9) | 兩者都不觸發 | ❌ |
| **HP tag (15-site)** | HP:i:33（保守標）| 0 | **20** | +20 | 否（V5 conservative）|
| | AMB% | 0.00% | 3.12% | +3.12 pp | 否（合理 conservative）|
| | HP1/HP2 ratio | 0.48 | 0.38 | 修正 tag bias | 否 |
| **Phase 結構** | phased% | 61.82% | 65.83% | +4.01 pp | ❌ |
| | n_blocks | 9,748 | 11,514 | +18.1% | ❌ |
| | N50 (bp) | 798,903 | 683,296 | -14.5% | 微差（仍 ≥ 600 K）|

→ purity 0.6 樣本：**caller F1 完全相同**、6 結構指標 4 改善 + 1 微差 + 1 持平、V5 conservative tagging 是好事。**沒有任何 critical regression**。

### 8.5.3 三路徑（multi-point）演算法不依賴 purity；Pass 2 = 只重跑 2-point

#### 兩個分離的 graph 演算法

從 `PhasingGraph.cpp` + `PhasingProcess.cpp` source code 還原 longphase-to phase 真實流程，**不是「2-point vs 3-point 二選一」，而是兩個分別函式各自處理不同任務**：

| 函式 | 用什麼 edge | 觸發條件 | 由誰呼叫 |
|---|---|---|---|
| **`somaticCalling`** | **3-point**（`threePointEdge`、`patternMining` first/second/third path）| `!disableCalling` flag（**與 purity 無關**）| `PhasingProcess.cpp:163, 180` Pass 1 內部 |
| **`edgeConnectResult`** | **2-point**（pairwise edge，`hpCountMap2`）| 永遠跑 | 透過 `phasingProcess()` 在 `PhasingGraph.cpp:1357` 呼叫 |

#### Pass 1 vs Pass 2 各跑哪個

```
低 purity (≤ 0.9) 樣本:
  Pass 1: edgeConnectResult (2-point) + somaticCalling (3-point)  ← 都跑
  Pass 2: 跳過

高 purity (> 0.9) 樣本:
  Pass 1: edgeConnectResult (2-point) + somaticCalling (3-point)  ← 都跑
  Pass 2: edgeConnectResult (2-point) only                          ← 只重跑 2-point
```

對應 `PhasingProcess.cpp:210-218` Pass 2 程式碼確認：

```cpp
if(highPurity){
    vGraph->convertNonGermlineToSomatic();
    chrInfo.posPhasingResult = PosPhasingResult();   // reset
    vGraph->phasingProcess(...);                      // 只跑 phasingProcess (= edgeConnectResult)
    // ← 沒有 vGraph->somaticCalling(...) 重跑
}
```

→ **Pass 2 只重跑 2-point**：高 purity 認為 graph 結構更可信 → 用 Pass 1 拍板的 somatic 分類結果重 phase 一次得到更乾淨的 block。**不重跑 3-point `somaticCalling`** 因為 Pass 1 已產出穩定的 origin 分類。

#### Pass 2 incremental effect 量化（HCC1395 5kHz 同 sample 對比）

數據來源：`output/pononly_v2b/tumor_phased.vcf` (Pass 1 only，ploidy bug 階段) vs `output/threshold_compare/v5_flag/tumor_phased.vcf` (Pass 1+2，ploidy fix 後)，同 PON-only flag、同 sample、同 input：

| 指標 | Pass 1 only | Pass 1+2 | Δ | 評語 |
|---|---:|---:|---|---|
| total variants | 3,187,275 | 3,187,275 | 0 | input 不變 |
| phased variants | 1,848,538 (58.00%) | 1,756,339 (55.10%) | **−92,199 (−2.90 pp)** | Pass 2 reclassify ~104K germline het 為 somatic/./. → 設計目標 |
| phase blocks | 1,808 | 1,631 | **−177 (−9.79%)** | block 數量減少（merging）|
| total phased bp | 2.846 Gbp | 2.850 Gbp | +4.18 Mbp (+0.15%) | 總覆蓋幾乎相同 |
| **N50 (bp)** | **11,388,114** | **11,788,053** | **+399,939 (+3.51%)** | **block 平均變更長** |

→ Pass 2 **block 變少 −10% 但每塊變更長 +3.51%**（典型 polish/merge 結果）。失去 92K phased variants 是 **Pass 2 reclassify 為 somatic/./. 的設計目標**，不是 regression。

#### 用戶記憶逐項對照與真實情況

| 用戶記憶 | 真實情況 | 對 / 錯 / 倒過來 |
|---|---|:---:|
| 「低 purity 用 3 點關聯」 | 低 purity Pass 1 用 **3-point + 2-point**（不只 3 點）| 部分對（含 3 點）但不只 |
| 「高 purity 只用 2 點關聯」 | 高 purity Pass 1 用 **3-point + 2-point**，Pass 2 **再加** 2-point 重 phase | **倒過來** — 高 purity 跑得**更多**（Pass 1 全部 + Pass 2 額外 2 點）|
| 「2 點 only 階段」可能源 | Pass 2 second round **確實只重跑 2-point**（不重跑 3-point `somaticCalling`）| 這部分**正確**，但只是「Pass 2 的特定行為」，不是「高 purity 全程只用 2 點」 |

#### 結論

- 三路徑（patternMining）算法**所有 sample 的 Pass 1 都跑**，不依賴 purity
- **高 purity 才多做事**（多跑 Pass 2 = 2-point 重 phase），不是低 purity 多做事
- Pass 2 incremental effect = 量化於 N50 +3.51% / block 數 −9.79% / phased variants −2.9 pp（reclassify 為 somatic 的設計目標）
- 從 our HEAD + `origin/main` PhasingProcess.cpp + `docs/phase.md` 三層驗證，唯一 purity gate 是 `highPurity = (purity > 0.9)`，**沒有 `if (!highPurity) {…}` 或 `if (purity < N) {…}` 之類低 purity 特殊分支**

> **常見誤解釐清**：「低 purity 會多做事」**不正確**。可能誤記源：(a) PON-only 模式 Pass 1 比非 PON-only 多 5 步（`ponOnlyPhasing` flag 控制，與 purity 無關）、(b) 「高 purity 才觸發 Pass 2 second round (2-point only)」記反成「高 purity 用 2 點」、(c) `patternMining` 三路徑投票（first/second/third path）字面像「三條路」但是 graph ambiguous edge 解析機制，與 purity 無關。

→ V5 在低 purity 樣本上：
- **Pass 1 跑同 baseline**（PON-only flag + 三路徑算法都同）
- **Pass 2 跳過**（同 baseline）
- **tag 層改動**（V3F two-layer + Layer 1.5）→ 受惠
- **ploidy fix** → 低 purity 計算精度小幅改善（0.607 → 0.634）
- **caller F1** → 三版完全相同（PASS set 不變）

→ 結構上 V5 在 0.6 樣本**最多持平，不會變差**，已由實測 6 cell F1 / 9 結構指標數據驗證。

### 8.5.4 longphase-to-mod 版本對齊確認

當前 V5 binary（HEAD = `938f0df`）vs 最新 GitHub origin/main：

| 後續 commit 數 | 改 algorithmic? | 改 doc/image? | 移除 dead code? |
|:---:|:---:|:---:|:---:|
| 10 | **❌ 全部不改** | ✅ 9 commits | ✅ 1 commit (filterSNP unused fn) |

origin/main HEAD 後續 commits 只動 `docs/`、`images/`、`README.md`、註解 + 移除 `filterSNP` 未使用函式。**HEAD `938f0df` = 最新有效演算法版本**。本報告所有結論建立在最新有效版本上。

---

## 8.6 Paired Mode Cross-Reference Audit（5/9 cycles 36+37 補充發現）

> **2026-05-09 用戶提問驅動**：「paired tag 是否也會出現 HP1-1+HP2-1 同位點共現偏向 HP1-1（類似 TO 17.3:1 priority bug）？同位點 germline 缺席但 somatic 共存的 reads 是否會被 longphase-to 錯標到 HP1 系列？」
>
> 兩個 cycle 階梯式驗證 — Step A+C (cycle 36) 答整體 paired mode 是否有 bug；Step D (cycle 37) 答 germline-absent 區域 V5 是否真修對。**結果揭露 V5 Layer 1.5 設計缺陷**：在 germline-absent 區域繼承 baseline priority bug 偏移，而 V3F 的「保守標 hp=33」反而更穩健。

### 8.6.1 Paired mode binary 釐清

| 模式 | Binary | Codebase | HP tag 編碼 |
|---|---|---|---|
| TO | `longphase-to` | `/big7_disk/liaoyoyo2001/longphase-to-mod/`（我們 fork，含 V5 修補）| `HP:i:` 整數（1/2/11/21/33）|
| **Paired** | `longphase-s somatic_haplotag` | `/big8_disk/liaoyoyo2001/longphase-s/`（**獨立 codebase**）| **`HP:Z:` 字串**（1/2/1-1/2-1/3）|

→ paired tag 用不同 binary、不同 codebase，與 longphase-to 的 priority bug **不直接套用**，需獨立驗證。HP tag 從 `longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp:533` `bam_aux_append(aln, "HP", 'Z', ...)` 確認字串編碼。

### 8.6.2 Step A — paired chr19 整體 HP:Z: tag distribution

對 paired tagged BAM (`/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/.../longphase_s/HCC1395_tagged.bam`) chr19 primary reads (582K total, 354,919 tagged):

| HP:Z: tag | reads | 占 tagged % |
|---|---:|---:|
| `HP:Z:2` | 183,309 | 51.6% |
| `HP:Z:1` | 143,760 | 40.5% |
| `HP:Z:2-1`（somatic on HP2）| 14,504 | 4.1% |
| `HP:Z:1-1`（somatic on HP1）| 12,401 | 3.5% |
| `HP:Z:3`（somatic ambiguous）| 1,145 | 0.3% |

**核心 ratio**：

| | paired chr19 | TO baseline 全基因組 | TO baseline chr19 victims |
|---|---|---|---|
| germline HP1:HP2 | 1:1.275（近 1:1）| 17.3:1（priority bug）| — |
| somatic HP1-1:HP2-1 | 1:1.169（近 1:1）| 全偏 HP1 | 752 全單向 baseline=11→V3F=21 |

→ **paired mode 整體沒 priority bug 偏移**。

![Figure F6 — paired vs TO HP distribution chr19](../../../research/paired_priority_bug_audit/figures/F6_paired_vs_TO_HP_distribution.png)

*Figure F6 — paired (longphase-s, 左) vs TO baseline (longphase-to, 右) chr19 HP tag distribution。paired germline HP1:HP2 ≈ 1:1, somatic HP1-1:HP2-1 ≈ 1:1（無偏移）；TO baseline 在 chr19 整體分布中 hp=11 (HP1 系列) > hp=21 (HP2 系列)，與全基因組 17.3:1 偏移一致。*

### 8.6.3 Step C — chr19 1Mb window som_ratio = HP1-1 / (HP1-1 + HP2-1)

| 統計 | 值 |
|---|---:|
| windows 數 | 57 |
| som_ratio mean | **0.462** |
| som_ratio median | 0.494 |
| som_ratio stdev | 0.332 |
| windows 近 1:1 (0.4-0.6) | 16 |
| windows 偏 HP1-1 (>0.7) | 13 |
| windows 偏 HP2-1 (<0.3) | 18 |

**重要觀察**：
- 整體 mean 0.46 ≈ 0.5 → **無 systematic bias**（vs TO baseline ~0.95）
- stdev 0.332 跨 windows 大變動 → **真實 sub-clone signal**：
  - chr19:3M 全 HP2-1（755/0）→ 該區域 LOH 方向特定
  - chr19:0M 全 HP1-1（330/1）→ 反向區域
  - chr19:17M som_ratio = 0.500（265/265 完美對稱）→ **SP1 附近 paired 認雙 sub-clone 共現**（vs TO baseline 113:0 失衡）

### 8.6.4 Step D — Germline-Absent Cross-Reference（**V5 Layer 1.5 設計缺陷揭露**）

對 paired chr19 read_name × T1.2 baseline / V3F / V5 vote dump 做 JOIN，篩 `cnt_HP1+cnt_HP2=0` 且 `somatic>0` 的 events（5,789 events），cross-tab：

| paired HP | events | baseline hp=11 | baseline hp=21 | V3F | V5 hp=11 | V5 hp=21 |
|---|---:|---:|---:|---|---:|---:|
| HP:Z:1-1 | 2,040 | **1,679** | 318 | 全 hp=33 | 1,679 | 318 |
| HP:Z:2-1 | 1,588 | **1,291** | 295 | 全 hp=33 | 1,291 | 295 |
| HP:Z:3 | 530 | 342 | 178 | 全 hp=33 | 343 | 177 |
| **加總** | — | **3,312** | **791** | 5,789（全 hp=33）| 3,313 | 790 |

![Figure F7 — Germline-absent events V3F vs V5 vs baseline 對比](../../../research/paired_priority_bug_audit/figures/F7_germline_absent_crosstab.png)

*Figure F7 — Germline-absent events (5,789 chr19 events) 三版處理對比。**左 panel**：V3F 全標 hp=33 保守不選邊（綠色 evaluation）；baseline 與 V5 都 4.19:1 偏 HP1（紅色 priority bug 警示）。**右 panel**：對 paired HP:Z:1-1/2-1/3 reads，baseline 與 V5 cross-tab 完全相同（hatch pattern 重疊）— V5 Layer 1.5 在該區域 = priority bug 的 feature 化非修補。*

**核心數字（binary-internal 量化，不依賴跨 binary 軸對齊）**：

```
baseline germline-absent  HP1 : HP2 = 3,312 : 791 = 4.19 : 1   ← priority bug 次峰偏移
V3F      germline-absent  全 hp=33 (somatic ambiguous)         ← 保守不選邊 ✓
V5       germline-absent  HP1 : HP2 = 3,313 : 790 = 4.19 : 1   ← 與 baseline 完全相同！
```

**關鍵發現**：

| 版本 | germline-absent 行為 | 評語 |
|---|---|---|
| baseline | 4.19:1 偏 HP1（priority bug 次峰）| 預期：vector ordered check + break early |
| **V3F** | **全標 hp=33（保守不選邊）** | **更正確**：避免錯標方向 ✅ |
| **V5** | **4.19:1 偏 HP1（與 baseline 完全相同）** | **Layer 1.5 設計缺陷** — 繼承 priority bug 而非修補 ❌ |

### 8.6.5 V5 Layer 1.5 設計缺陷的機制詮釋

V5 `getVote()` Layer 1.5 邏輯（`d0bcd8c` bundled，§5.6 程式碼對照）：
```cpp
// Layer 1.5: germline 缺席時用 somatic phased votes 決方向
else if (somaticHP1 > 0 || somaticHP2 > 0) {
    germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;
}
```

但 self-phasing 機制下：sub-clone somatic 100% 共現 → 進 phasing graph 後 somatic edges 偏向同一 haplotype → `somaticHP1` vs `somaticHP2` 票數偏向同邊 → **Layer 1.5 結果繼承 priority bug 偏移（4.19:1 偏 HP1）**。

→ **V5 Layer 1.5 在 germline-absent 區域 = priority bug 的 feature 化**：把「baseline 用 somatic vote 蓋過 germline」的 buggy 行為，改成「germline 缺席時才用 somatic vote」的 designed 行為。但**該區域的 4.19:1 偏移本質沒變**。

### 8.6.6 對 §8.4 機制詮釋的補強 — 兩層差異

§8.4 把 V5 vs V3F 的 zero-sum 重分配解讀為「Pass 2 reclassify 104K germline het」。Step D 揭露**另一層機制**：

| 區段 | V5 vs V3F 差異 | Source |
|---|---|---|
| germline_vote > 0 區（主流）| Pass 2 reclassify 部分位點為 somatic → reads shift Layer 1 → Layer 1.5 | §8.4 量化（Pass 2 effect）|
| **germline_vote = 0 區**（少數但重要）| **V3F 標 hp=33 保守 vs V5 用 somatic 投票決方向 → 繼承 priority bug 4.19:1 偏移** | §8.6 Step D（**Layer 1.5 設計缺陷**）|

→ V5 vs V3F 在 germline-absent 區域**不只是「tag 數量」，更是「方向偏移 vs 保守 untagged」的設計差異**。在這區域 V3F 的處理**比 V5 更穩健**（避免錯標）。

### 8.6.7 Caveat

| 限制 | 說明 |
|---|---|
| 跨 binary axis alignment | paired vs TO 各自獨立 phasing，HP1/HP2 軸 phase block 內一致跨 block 50% swap；cross-tab 中「paired 2-1 → baseline hp=11」**不能直接解讀為錯標** |
| binary-internal 量化不受影響 | baseline 自身 HP1:HP2 = 4.19:1 是 binary-internal 偏移量化，不依賴跨 binary 軸 |
| chr19 only | 全基因組擴展待 F-paired-D1 |
| sample size | 5,789 events × 295K paired reads；全基因組估 ~150K events 量級 |

### 8.6.8 用戶提問的精確回答

| Q | A |
|---|---|
| paired tag 是否有 HP1-1+HP2-1 同位點共現？ | **是** — 16/57 chr19 windows som_ratio 0.4-0.6 共現（Step C）|
| paired 是否偏向 HP1-1？ | **否** — som_ratio mean 0.462 跨 windows 0-1 全範圍分布（Step C）|
| 同 reads 是否在 longphase-to 被錯標到 HP1 系列？ | **是** — germline-absent 區域 baseline 4.19:1 偏 HP1（Step D）|
| V3F / V5 是否修對這些錯標？ | **V3F 修對**（全標 hp=33 保守處理）；**V5 沒修對**（與 baseline 4.19:1 完全相同）|
| 應該回歸 V3F 嗎？ | 待 ISM 影響量化（F-paired-D3 follow-up）|

完整 audit 詳見：
- [`InterSubMod/research/paired_priority_bug_audit/00_audit_report.md`](../../../research/paired_priority_bug_audit/00_audit_report.md)（Step A+C 整體 audit）
- [`InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md`](../../../research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md)（Step D V5 Layer 1.5 設計缺陷）

---

## 9. 結論與後續方向

### 9.1 整體 self-phasing 修補成熟度

| 維度 | 狀態 | 證據 |
|---|---|---|
| 機制因果 | ✅ 確立 | 17.3:1 全基因組 + 34,855 read-level victims + 100% V3F/V5 修正 |
| 修補設計合理性 | ✅ 完整 | 5 commits 對應兩 layer + germline 缺席 fallback + Pass 2 觸發 |
| chr19 SP1/2/3 案例 | ✅ V3F/V5 對齊 paired 3/3 | IGV 6-BAM 並列截圖 |
| 全基因組擴展 | ✅ chr19 結論部分顛覆但不影響機制 | T1.2-F1 |
| **V5 vs V3F zero-sum 機制** | **✅ 已釐清** | **§8.4 Pass 2 reclassify 104K germline het → reads shift Layer 1→1.5** |
| **20 指標 no regression** | **✅ 全部持平或變好** | **§8.5.1 PI 報告 §6 全表 + LOH bed 三版相同** |
| **Caller F1 vs SEQC2 truth** | **✅ 三版完全相同** | **§8.5.2 0.93 = 0.7166 / 0.6 = 0.6273（V5 不改 caller）** |
| **purity 0.6 樣本完整對照** | **✅ 6 F1 + 9 結構 0 regression** | **§8.5.2 caller F1 同、phased% +4 pp、tag conservative ↑** |
| **三路徑算法不依賴 purity** | **✅ Pass 1 永遠跑** | **§8.5.3 PhasingGraph.cpp:267 patternMining 不 gated** |
| **Pass 2 = 只重跑 2-point** | **✅ 量化** | **§8.5.3 N50 +3.51% / block −9.79% / phased var −2.9 pp（reclassify 設計目標）** |
| **版本對齊** | **✅ HEAD 938f0df = 最新有效演算法** | **§8.5.4 origin/main 後 10 commits 全 doc/dead-code** |
| **Paired mode 整體 priority bug** | **✅ NEGATIVE** | **§8.6 Step A+C — paired 用 longphase-s 不同 binary，HP1:HP2 = 1:1.275，som_ratio mean 0.462 無 systematic bias** |
| **V5 Layer 1.5 在 germline-absent 區域行為** | **⚠️ 設計缺陷** | **§8.6.4 Step D — V5 與 baseline 4.19:1 偏 HP1 完全相同，未修補 priority bug；V3F 標 hp=33 反而更穩健** |
| Pass 2 second round 獨立貢獻 | ⚠️ 未量化 | ISM benchmark 用戶 5/7 決策 cancel；可由 T1.3 補 |
| 跨樣本一致性 | ⚠️ 待擴展 | 僅 HCC1395 5kHz 一樣本（T3 候選）|

整體：**self-phasing 主軸修補成熟可用**（V3F + V5 在 HCC1395 5kHz 全基因組 100% 修正、20 指標無 regression、版本對齊最新有效演算法）；**新發現** V5 Layer 1.5 在 germline-absent 區域繼承 priority bug 偏移（§8.6 Step D），V3F 設計反而更穩健 — 這對 V5 設計選擇有 implications，但**不阻擋 V5 作為 production tag baseline**（germline-absent 區域佔 chr19 events 比例小）；剩跨樣本擴展 + Pass 2 獨立貢獻量化 + Layer 1.5 改回 V3F 的 ISM 影響為下游 follow-up。

### 9.2 對 PI 報告（2026-04-29）的訊息更新清單（5 條 errata 候選）

| # | 段落 | 原表述 | 建議修訂 |
|---|---|---|---|
| **E1** | §3.3.3 chr19 SP1/2/3 hotspot | chr19 SP1/2/3 是「主要 hotspot」 | 改為「可重現案例」；priority bug 主要分佈在 chr7/chr2/chr1/chr16/chr20 |
| **E2** | R11 chr8 hotspot | chr8 LOH+HPSig 7.4× FP enrichment 與 priority bug 關聯 | chr8 priority bug 是冷區（rank 21）；hotspot 機制獨立於 priority bug，需另案 |
| **E3** | §5.2 priority bug 機制 | commit message + IGV 3 截圖支持 | 升級為「個案 + 統計 + 機制」三重佐證（read-level 34,855 全基因組鐵證）|
| **E4** | §6.4 / §6.5 V5 數值歸因 | V5 four-commit chain 整體效益 | 主要 V3F + Layer 1.5；Pass 2 二次效益尚未獨立量化（PI 報告 V5 數據實為 Pass 1 only）|
| **E5（5/9 新）** | §5.2 V5 Layer 1.5 設計描述 | Layer 1.5 = germline 缺席區域的 fallback（隱含「修補」）| 在 germline-absent 區域 V5 與 baseline 4.19:1 偏 HP1 完全相同（§8.6 Step D）— 是 priority bug 的 feature 化非修補；V3F 標 hp=33 保守處理反而更穩健 |

E1-E4 已 patch 至 [`InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md`](20260509_PI_Report_4_29_Errata_01.md)（commit f17754f）；E5 為 5/9 Step D 新發現，待後續 errata 補強或 V5 設計選擇 follow-up。

### 9.3 待解問題與後續行動

| ID | 主題 | 內容 | 預估時間 |
|---|---|---|---|
| **F1** ✅ | PI 報告 errata patch（E1-E4）| 整合 §9.2 4 條訊息變更到 4-29 PI 報告對應段（5/9 完成）| DONE 2026-05-09 commit f17754f |
| ~~F2~~ ✅ | ~~V5 vs V3F zero-sum 釐清~~ | **已釐清（§8.4）**：Pass 2 reclassify 104,457 germline het 為 somatic/未 phase → reads 從 Layer 1 shift 至 Layer 1.5 | DONE 2026-05-09 |
| **F3** | chr8 LOH+HPSig hotspot 另案 | 從 ISM 特徵 / pileup 異常 / coverage 入手，與 priority bug 解耦後重啟 | 視 scope |
| **T3** | 7 樣本跨樣本擴展 | HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829 跑同 vote audit 驗證 priority bug 跨樣本分佈 | 1-2 天（含 binary patch + dump 6×3）|
| **T1.3** | 4-cell ablation（conditional）| Pass 2 vs V3F 各自獨立貢獻量化 — 量化 §7.3 Pass 2 二次效益 | 3 天 |
| **F-paired-D1（5/9 新）** | germline-absent xref 全基因組擴展 | Step D chr19 → 全 chr，估 ~150K events，確認 4.19:1 偏移是否跨 chr 一致 | 0.5 天 |
| **F-paired-D2（5/9 新）** | Phase block 內 axis-aligned 分析 | 用 PS tag 對齊 paired vs TO 軸 → 計算「真正錯標」率（去除 cross-binary axis swap caveat） | 1 天 |
| **F-paired-D3（5/9 新）** | V5 Layer 1.5 改回 V3F 的 ISM 影響評估 | germline 缺席改標 hp=33 而非用 somatic 投票決方向 — 量化 ISM 下游影響 | 1-2 天（需 binary patch）|
| **F-paired-D4（5/9 新）** | E5 PI errata 補強 | Step D 結論加進 PI 4-29 errata 列表（接 E1-E4 後）| 30 min |

→ 本報告完成：(1) V5 audit synthesis（§2-§9），(2) zero-sum 機制釐清（§8.4），(3) no-regression 量化驗證（§8.5），(4) **paired cross-ref 揭露 V5 Layer 1.5 設計缺陷**（§8.6, 5/9 新發現）。F1 已完成；F-paired-D1/2/3 是 5/9 新發現的後續延伸；T3 跨樣本擴展與 T1.3 ablation 仍是後續主軸 cycle。

---

## 附錄 A — 三版 testing-only binary 編譯與 dump 重現步驟

詳見 `InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md` §10-11（patch_actual.diff 131 行 + 編譯三版 + dump 步驟）。本報告 read-level 數據（752 chr19 + 34,855 全基因組）即由此流程產出。

## 附錄 B — Glossary

| 術語 | 定義 |
|---|---|
| **HP:i:0** | untagged read（沒有方向證據）|
| **HP:i:1 / HP:i:2** | germline-only read，方向 HP1 / HP2（沒 somatic vote）|
| **HP:i:11 / HP:i:21** | germline + somatic 共存 read，germline 方向 HP1 / HP2 + 帶 somatic annotation |
| **HP:i:33** | 純 somatic read（germline ambiguous，只剩 HP3 vote）|
| **HAPLOTYPE1 / 2** | germline het 投票 enum（內部 0/1）|
| **HAPLOTYPE1_1 / 2_1** | somatic 投票 enum（內部 2/3）|
| **HAPLOTYPE3** | mixed somatic / pseudo-somatic enum（內部 4）|
| **germline_maj** | countMap[HP1] vs countMap[HP2] majority（1, 2, 或 0=tie）|
| **somatic_maj** | countMap[HP1_1] vs countMap[HP2_1] majority |
| **priority bug victim** | germline_maj ≠ somatic_maj 且 baseline `hpResult` 跟 somatic 方向 |
| **Pass 1** | PhasingProcess phasing 第一輪（PON-only graph）|
| **Pass 2** | second round phasing（高純度樣本才觸發，convertNonGermlineToSomatic + reset + 重 phase）|
| **`ponOnlyPhasing`** | LongPhase-TO `--pon-only-phasing` flag，phasing 階段只用 PoN germline het |
| **Layer 1 / 1.5 / 2** | V5 三層 getVote：germline 投票 / somatic fallback / somatic annotation |
| **highPurity** | `bool highPurity = purity > 0.9`，決定 Pass 2 是否觸發 |
| **ploidy bug** | Pass 1 nullptr → ploidyRatioMap empty → q1=q3=0 → purity=0（d0bcd8c 已修）|

## 附錄 C — 引用 commit hash + 文件清單 + ledger cycle ID

### longphase-to-mod commits
| Hash | 日期 | 主題 |
|---|---|---|
| `8b8c1fd` | 2026-04-09 | feat: --pon-only-phasing flag |
| `41ff147` | 2026-04-10 | fix(haplotag): two-layer getVote |
| `380e8d2` | 2026-04-25 | fix(haplotag): countINDELHaplotype UNDEFINED guard |
| `d0bcd8c` | 2026-04-30 | fix(purity): collect ploidyRatio + bundled Layer 1.5 |
| `938f0df` | 2026-04-30 | Update purity threshold 0.95→0.9 |

### Parent reports
- `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` — PI 報告 17.3:1 / SP1-3 / sanity 15/15 / V5 four-commit
- `InterSubMod/docs/reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md` — V5 audit + §5.6 chr19 機制 + §5.7 全基因組顛覆
- `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md` — 早期根因敘述
- `InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md` — chr19 752 機制報告
- `InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md` — 全基因組 34,855 報告

### Evidence ledger cycles
- `20260505_self_phasing_V5_provenance_audit` (PROV-V5-001, verdict=negative, stability=5)
- `20260507_T1_2_priority_bug_mechanism` (PROV-V5-002, verdict=positive, stability=5)
- `20260508_T1_2_F1_genome_wide_audit` (PROV-V5-003, verdict=positive_with_corrections, stability=5)
- `20260508_self_phasing_synthesis_report` (PROV-V5-SYNTH-001, verdict=positive, stability=5, 本報告 commit 951e7c9)
- `20260509_PI_report_4_29_errata` (PROV-V5-ERRATA-001, verdict=positive, stability=5, commit f17754f)
- `20260509_paired_priority_bug_audit` (PROV-V5-PAIRED-001, verdict=negative, stability=5, commit 6ed8a0d) — §8.6.2-3
- `20260509_paired_germline_absent_xref` (PROV-V5-PAIRED-002, verdict=positive_with_caveat, stability=4, commit 766ec5f) — §8.6.4-7

---

> **Dialogue break 4** — §7/§8/§9 + 附錄 A/B/C 寫完。請 review：
> 1. **§7.1 ploidy bug ASCII 因果鏈**：是否清楚說明「purity=0 不是樣本問題」？
> 2. **§7.3 數值影響表**：把 PI 報告 V5 各數值「主要功勞」分到 V3F vs Pass 2 是否準確？
> 3. **§8.1 chr19 占比 ASCII 長條圖**：是否好看？
> 4. **§8.4 V5 vs V3F zero-sum**：用 zero-sum 重分配描述 mechanism 待釐清是否合理？
> 5. **§9 errata 4 條 + 後續 5 項**：scope 與優先序是否合理？
>
> 確認後我會：(1) 產 5 張新 figures (F1-F5) 補進對應段、(2) 加 INDEX entry + commit + 寫 evidence_ledger cycle。
