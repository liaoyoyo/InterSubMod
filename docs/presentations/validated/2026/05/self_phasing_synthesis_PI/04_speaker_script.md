<!--
build_date: 2026-05-10
agent: pptx-build P5 speaker script
status: in_progress (awaiting C5 user ack)
parent: 03_slide_layout_script.md
target_duration_min: 20
total_slides: 25 (22 main + 3 backup)
budget_per_slide: 75-90 sec / 中 300-360 字
-->

# Self-Phasing PPT — Speaker Script (P5)

> 每 slide Tier 2 speaker note 75-90 sec ≈ 中 300-360 字。Tier 3 標 [ORAL-OPTIONAL] 依現場時間決定。

---

## Slide 01 — Cover [30 sec / 中 ~120 字]

「謝謝各位。今天 20 分鐘的報告主題是 self-phasing 整合觀察與 V5 Layer 1.5 設計缺陷揭露。這份內容整合 longphase-to-mod 5 commits 的修補成熟度，以及 5/9 paired cross-ref 新發現的 V5 Layer 1.5 設計問題。報告目的：協助 PI 決策 V5 是否作為 production tag baseline，以及是否啟動 F-paired-D3 follow-up cycle。所有素材來自 5/8 整合報告 1,211 行加 5/9 errata，以及 paired audit Step D。」

[ORAL-OPTIONAL] 5/8 主報告路徑 + 7 figures commit hash

→ next: slide 02 TL;DR

---

## Slide 02 — TL;DR [90 sec / 中 ~360 字]

「先看 thesis：今天有兩個焦點。**正面：修補主線確立。**baseline LongPhase-TO 全基因組 HP1:HP2 = 17.3:1，是個極端 systematic bias；V3F + V5 兩層修補在 read-level 對 chr19 752 個 priority bug victim 修正率 100%、全基因組 34,855 個 victim 一樣 100%，paired ground truth 對齊 +13.3 pp，20 指標 0 regression，caller F1 三版完全相同。**反面：V5 Layer 1.5 在 germline-absent 區仍未修對。**5/9 paired cross-ref 揭露：在 germline 缺席區域（占 chr19 events 約 5%），V5 與 baseline 的 4.19:1 偏 HP1 完全相同 — 也就是 V5 Layer 1.5 把 priority bug 從 bug 變成 designed feature，但偏移本質沒變；反而是 V3F 標 hp=33 保守不選邊更穩健。整體成熟度 12 維度 10 個 ✅、1 個 ⚠️、2 個待跑；V5 仍可作 production baseline，但 F-paired-D3 量化後可能要回 V3F default。」

[ORAL-OPTIONAL] V3F = 41ff147 / V5 = d0bcd8c + 938f0df / paired Step D = 766ec5f

→ next: slide 03 全基因組 17.3:1

---

## Slide 03 — 全基因組 17.3:1 偏移 [90 sec / 中 ~340 字]

「先解釋 17.3:1 為什麼是真實 systematic bias，不是樣本性質。baseline LongPhase-TO 在 HCC1395 5kHz 全基因組層觀察：HP1 reads 614,000 vs HP2 reads 35,500，HP1 占比 94.6%（vs 隨機預期 50%）。為什麼能斷定這是 engineering bias 不是 sample artifact？三條獨立論證：第一，生物學上 tumor 的 somatic ALT 屬某 sub-clone 的真實 haplotype，跨 23 條染色體不該系統性偏 HP1；第二，cnLOH artifact 只影響單一 chr，但這 94.6% 是跨 23 chr 一致；第三，paired tumor-normal 同樣 reads 走 paired pipeline HP1:HP2 ≈ 1:1 — 只有 TO 模式出現 17.3:1。這三條獨立佐證不可能同時都是 sample 效應，所以 17.3:1 是 LongPhase-TO 的 systematic engineering artifact。next slide 看具體位點層的鐵證。」

[ORAL-OPTIONAL] cnLOH 機制細節 / 23 chr 一致性表 / paired pipeline 程式碼差異

→ next: slide 04a SP1 IGV

---

## Slide 04a — SP1 chr19:17,565,944 IGV [60 sec / 中 ~250 字]

「全基因組 17.3:1 是平均值；用 IGV 6-BAM 並列在 chr19 找到三個近 100% 失衡位點，先看 SP1 chr19:17,565,944。baseline panel 113 條 reads 全部集中於 HP1，HP2 stack 為 0。並列觀察 V2b、V3F、V5 三個修補階段，最後 V5 翻轉至 HP2 主導，與 paired tumor 的 HP2 方向完全一致。為什麼能排除是噪音、caller、alignment 問題？因為 baseline 與 paired 不是衰減而是完全反向，V5 修正後與 paired ground truth 重合，這是 read assignment 強制集中的鐵證。next slide 看 SP2 / SP3 兩個同模式案例。」

[ORAL-OPTIONAL] 6-BAM 並列順序解釋 / V2b 中間階段意義

→ next: slide 04b SP2/SP3

---

## Slide 04b — SP2 + SP3 IGV 並列 [60 sec / 中 ~260 字]

「SP2 chr19:12,452,332 baseline 109:1，SP3 chr19:12,467,180 baseline 108:0，兩個位點與 SP1 同模式：baseline 全 HP1，V5 翻 HP2 對齊 paired ground truth，3/3 全對齊。三個 SP 位點都在 chr19:12-17M 區段，這個區段在 next slide 8 的 chr19 752 read-level victims 1Mb hotspot 散點圖中也是 enrichment 最高的區域，互相佐證 — read-level 個案與 IGV 截圖屬同一機制的不同層級觀察。三 SP 都 V5 修對引出四個自然問題：為什麼 baseline 把 read 都集中一邊？17.3:1 在 read 層級有多少個案、分佈在哪？baseline → V3F → V5 三版各修什麼？這些都是後面 §機制 / §量化 / §修補 / §驗證 章節要回答的。next 進入機制。」

[ORAL-OPTIONAL] paired_T 與 paired_N 對照細節

→ next: slide 05 phasing 層球員兼裁判

---

## Slide 05 — phasing 層球員兼裁判 [120 sec / 中 ~360 字]

「機制有兩層 bug。先講 phasing 層 — 球員兼裁判隱喻。正常 phasing graph 的物理基礎：germline het 在 HP1/HP2 兩 stream 50/50 隨機分佈，read 看到 (A,C) 就 assign HP1，看到 (G,T) 就 HP2，phasing graph 把 het 連成 haplotype，**設計上只 germline 該進 graph**。但 TO 模式沒 paired normal 可以區分 germline / somatic，只能用 PoN — Panel of Normals 這個多正常樣本建構的 reference set。**未在 PoN 內的位點預設當 germline 處理 → somatic mutation 也被當 germline 進 graph。**這裡致命：germline het 同 sub-clone 內 50/50 共現，但 somatic 屬某 sub-clone 的 100% 全帶 — 100% 共現比 50% 強得多 → 多個 somatic 之間 edge weight 暴漲 → graph 自我增強迴圈：somatic 越強 → 偏該 haplotype → 該 haplotype 越強。**結果：germline 真實訊號被 overrule。**「球員兼裁判」就是這意思：somatic 應該是被 phase 的對象（受 germline 仲裁），現在反過來主導 graph。修法叫 PON-only flag — commit 8b8c1fd — phasing 階段 graph 只放 PoN 內的 germline het，somatic 不進 graph，等 graph 拍板後再用結果反過來分類 somatic。但這只解 phasing 層，tag layer 的 priority bug 還是壞的，下張 slide 講。」

[ORAL-OPTIONAL] TO vs paired PoN 對照 / 自我增強迴圈視覺化 / Pass 1 vs Pass 2 預告

→ next: slide 06 tagging 層 priority bug

---

## Slide 06 — tagging 層 getVote priority bug [120 sec / 中 ~360 字]

「即使 phasing graph 修對了 PON-only flag，tagging 階段 getVote 函式還是會出錯。看 Figure F1 — baseline 程式碼用 vector 順序檢查：① somatic pair 在最前面 ② mixed pair ③ germline pair 在最後。for 迴圈遇到第一個非空 vector 就 break early，後面的 germline pair 永遠看不到。看右邊真實 read 範例：germline HP1=0、HP2=5（5 票主導）、somatic HP1_1=1（只 1 票）、HP2_1=0、HP3=0。baseline 的投票流程：第一輪檢查 (HP1_1, HP2_1) → HP1_1=1>0 → 決定 hp=11 → break。第二輪 (HP3, HP2_1) skipped；第三輪 (HP1, HP2) 也 skipped — germline 5 票完全被忽略。**正確答案應該是 hp=21**（germline HP2=5 主導，附 somatic HP1_1=1 標 21）。為什麼全偏 HP1？因為 tumor sub-clone 的 somatic 都同方向 — §3.2 解釋的 100% 共現 — baseline 看到的 somatic 票多偏 HP1_1 → priority bug 把所有受影響 reads 標 HP:i:11 系列 → 17.3:1 偏移在 tag layer 形成。注意：tag layer 與前一張 slide 講的 phasing layer 是不同層 bug，必須分別修補，不能合併單一 commit 解決。」

[ORAL-OPTIONAL] enum HAPLOTYPE1_1=2 vs HP tag int=11 比較 bug 細節 / 5-vote countMap 結構 / 真實 read_name

→ next: slide 07 兩層對應表

---

## Slide 07 — 兩層 bug 兩層修補對應 [90 sec / 中 ~310 字]

「總結兩層機制：phasing 層的球員兼裁判 — somatic 進 graph — 修補在 commit 8b8c1fd PON-only flag；tagging 層的 priority bug — vector 順序錯加 break early — 修補在 commit 41ff147 V3F two-layer getVote。但 V3F 還需要兩個 commit 才完成 V5：380e8d2 INDEL guard 補 OOB undefined behavior、d0bcd8c bundled Pass 2 ploidy fix 加 Layer 1.5 fallback、938f0df threshold 從 0.95 降到 0.9。為什麼不能只用 PON-only flag？因為解了 phasing 層循環依賴，但 tag layer 的 getVote 還是壞的，read 99.9% 仍然標 HP:i:11 系列 — 17.3:1 偏移仍在。為什麼 V3F 不夠還要 V5？因為 V3F 只有 Layer 1，germline 缺席區（cnLOH 或 amplicon hotspot）reads 全部 untagged 流失；V5 加 Layer 1.5 fallback。但 V5 Layer 1.5 在 germline-absent 區有自己的 caveat，後面 slide 16 詳述。」

[ORAL-OPTIONAL] 5 commit 順序 / V3F 命名來歷 / 跨層交互機制

→ next: slide 08 chr19 752 read-level

---

## Slide 08 — chr19 752 read-level 鐵證 [120 sec / 中 ~360 字]

「機制講完，看量化鐵證。對 baseline / V3F / V5 三版做 testing-only binary patch 加 --debug-vote-dump flag，dump 每條 read 經 getVote 後的 5-vote countMap 加 hpResult。chr19 子集（HCC1395 5kHz）三版各 549,206 rows × 3-way merged 1,069,832 events，篩 germline_majority ≠ somatic_majority 且 baseline hpResult 跟 somatic 方向 = **752 priority bug victims**。**V3F 修正比例 100%**（全改向 germline_majority），**V5 修正比例 100%**。全 752 條都是 baseline hp=11 → V3F=21 → V5=21 單向修正，無一條反向 — 這跟全基因組 17.3:1 偏移完全對應。4-path 驗證（T1.2 plan 設計嚴格驗證）：① 個案 trace 752 條 PASS；② 1Mb 區域聚集 chr19:30M (215) 加 27M (133) 共佔 chr19 victim 46% PARTIAL；③ Somatic density 共變 high vote ≥5 = 0 受害；low vote 才觸發 — 反向但有意義（高票 sub-clone 一致已對齊 graph，低票才被 priority bug 蓋過）；④ V3F/V5 修正後消失 PASS。**3 PASS + 1 PARTIAL = 機制因果確立。**Figure F4 顯示 chr19 1Mb hotspot 散點，30M 加 27M 區段對齊 SP1/2/3 的 chr19:12-17M。next 看全基因組擴展。」

[ORAL-OPTIONAL] 4-path detail / read_name 真實 case 5 條 / SP1/2/3 對應 chr19 hotspot 區段

→ next: slide 09 全基因組 34,855

---

## Slide 09 — 全基因組 34,855 victims (46×) [120 sec / 中 ~360 字]

「對同三版 binary 跑全基因組 vote dump（每版 ~40 min，總 18.9M tagged reads），dump 大小 744/687/687 MB（gzipped）。**全基因組 priority bug victims = 34,855**，是 chr19 752 的 46.4 倍。V3F / V5 修正率仍 100% / 100% 一致。但 per-chr 分佈推翻原 chr19 pilot 結論：**chr19 752 占全基因組僅 2.16%，rank 19**（24 條染色體中），不是主要 hotspot。主要 hotspot 是 chr7 (3,508 / rank 1) / chr2 (2,792 / rank 2) / chr1 (2,674 / rank 3) / chr16 / chr20。chr19 SP1/2/3 是「可重現案例」而非「主要分佈位置」。順帶說 chr8 — MEMORY 記錄 chr8 LOH+HPSig 是 7.4× FP enrichment 的 hotspot，但 chr8 priority bug enrichment 只 0.34× genome avg，rank 21 是冷區。**這兩個是不同 layer**：chr8 LOH+HPSig hotspot 是 ISM 下游 false-positive 富集，priority bug 是 longphase tagging 投票錯，兩者跨層獨立沒有直接因果關聯，所以 V3F/V5 修對 priority bug 後 chr8 hotspot 不會自動消失，要另案處理。Figure F2 雙 panel：左按 victim N 排序，右按 enrichment ‰ 排序，chr19 紅標 / chr8 藍標凸顯。」

[ORAL-OPTIONAL] 全 chr enrichment ‰ 表 / chrY 小 N 高 ‰ 解釋 / chr8 不同 layer 詳細分析

→ next: slide 10 5 commits timeline

---

## Slide 10 — 5 commits 時間軸 + 三色標 [120 sec / 中 ~340 字]

「修補設計演進。看 Figure F3 timeline — 5 commits 漸進完成。**8b8c1fd** 04-09 PON-only flag — phasing 層藍色；**41ff147** 04-10 two-layer getVote — tagging 層綠色，**這是修偏移的關鍵 commit ★**；**380e8d2** 04-25 INDEL guard — tagging 層綠色，補 countINDELHaplotype 在 somatic 位點的 OOB undefined behavior；**d0bcd8c** 04-30 跨兩層紫色 — Pass 2 ploidy fix 加 bundled Layer 1.5 加 countSNP guard，**這是唯一跨兩層的 commit**；**938f0df** 04-30 threshold 0.95 降 0.9 — phasing 層藍色 Pass 2 觸發。V3-Fixed = baseline + 41ff147 + 380e8d2，V5 = V3F + d0bcd8c + 938f0df。累計改動 ~155 行 tagging-layer 加 ~40 行 phasing-layer，HaplotagProcess.h:66-68 介面契約零變動。為什麼必須 stacking 不能合併？因為 layer 獨立，PON-only 解 phasing 但 tag 仍偏；V3F 解 tag 但 germline 缺席區 untagged；V5 補 Layer 1.5 fallback 加 Pass 2 觸發。」

[ORAL-OPTIONAL] 各 commit 對應 layer 細節 / 為何不合併 / cherry-pick 自 upstream zhenyu

→ next: slide 11 getVote 三版 code

---

## Slide 11 — getVote 三版程式碼對照 [120 sec / 中 ~360 字]

「三版程式碼 side-by-side。**baseline 紅底**：vector keys 順序 ① somatic FIRST ② mixed ③ germline LAST，for 迴圈遇到第一個非空就 break，priority bug 在這。**V3F 綠底**：拆 explicit Layer 1 germline only — `if (germlineHP1>0 || germlineHP2>0)` 決方向，gR=1 或 2；Layer 2 somatic annotation — `if (somaticTotal>0)` 標 hp=11/21/33。Layer 1 的 germline 永不被 somatic overrule。**V5 綠底加黃 highlight**：保留 Layer 1（同 V3F），新增 Layer 1.5 fallback — `else if (somaticHP1>0 || somaticHP2>0)` germline 缺席時用 somatic phased votes 決方向 — Layer 2 同 V3F。底部白話三段：baseline 順序錯 → V3F 拆 germline 主導 + somatic 標籤 → V5 補 germline 缺席 fallback。**但黃 highlight 是 caveat 預告**：Layer 1.5 在 germline-absent 區會繼承 priority bug 偏移，slide 16 詳述。」

[ORAL-OPTIONAL] enum HAPLOTYPE1_1=2 / HP tag int=11 比較 / V3F bonus 修：hpResult ≠ HAPLOTYPE1_1 比較失誤

→ next: slide 12 SP 修正 + zero-sum

---

## Slide 12 — SP1/2/3 修正後對齊 paired + 全基因組 17.3:1 → ~1:1 [90 sec / 中 ~340 字]

「驗證鐵證雙層：個案層 + 全基因組層。**個案層**：SP1/SP2/SP3 baseline 113:0 / 109:1 / 108:0，V5 修正後三位都翻轉至 HP2 主導，與 paired ground truth 對齊 3/3 ✅。**全基因組層**：HP1:HP2 從 17.3:1 → ~1:1 消除偏移；94.6% somatic→HP1 → ~50% balanced；15-site Problem PS（含 SP1/2/3）48.5% → 52.0% +3.5 pp 看似小但機制顯著（這幾個極端失衡位點本身是最難改善的）。Figure F5 4 象限視覺化 zero-sum 重分配：germline=0 區 +560,881 reads（V5 Layer 1.5 觸發 tagging）/ germline>0 區 -560,881 / 總和等於 0。**這個 zero-sum 不是 bug 是 Pass 2 reclassify 設計目標** — Pass 2 把 104,457 個 germline het 重新分類為 somatic 或 ./.，受影響 reads 從 Layer 1 路線 shift 到 Layer 1.5 路線。priority bug 修補的 100% verdict 不受影響。」

[ORAL-OPTIONAL] V2b / V3F 中間版本 / Layer 1.5 zero-sum 詳細機制 / Pass 2 reclassify 量化

→ next: slide 13 20 指標 no regression

---

## Slide 13 — 20 指標 no regression [120 sec / 中 ~360 字]

「對 V5 vs baseline 跑 5 大類別 20 個 pre-registered metrics — pre-registered 是說事先定義好不挑後選的避免 cherry-picking 嫌疑。① **ISM aggregate 3 項**：TP_rate +0.005、HP_Ratio median 0.788→0.574 — 注意 HP_Ratio 下降不是變差，是 tag bias 修正後 ratio 正常化 — Potential_LOH +3.5 pp。② **HP_Ratio AUC 2 項**：All -0.005 在隨機區間內、Inner +0.002。③ **Methylation 6 feature** 全 ±0.01 內持平。④ **Paired GT concordance 4 項 ⭐**：clean PS +8.3 pp / 15-site Aggregate +6.65 pp / 15-site Clean PS +13.3 pp / 15-site Problem PS +3.5 pp — 這 4 項都是 V5 顯著改善的 ground truth metric。⑤ **HP / LOH 結構 5 項 ⭐**：Phase block N50 +99.7%（4,061 → 8,109）、Phased rate +23.6 pp、執行時間 1.36× 快、LOH regions 完全相同（Jaccard=1.0）、LOH 總 bp 完全相同。**結論：20/0 no regression，6 項 ⭐ 顯著改善 +6.65 ~ +99.7%；其餘 ±0.01 內持平。**V5 全面 production-ready。注意這個 baseline 是 PI 報告 4-12 V5 BAM = Pass 1 only（ploidy bug 讓 Pass 2 從未觸發），主要功勞在 V3F 加 Layer 1.5。」

[ORAL-OPTIONAL] methylation 6 feat 列表 / HP_Ratio 0.788→0.574 詳解 / LOH bed Jaccard=1.0 機制

→ next: slide 14 caller F1 + cliffhanger

---

## Slide 14 — Caller F1 三版相同 + Cliffhanger [120 sec / 中 ~360 字]

「Caller F1 vs SEQC2 truth set v1.2.1。HCC1395 5kHz @ 0.93 purity：A1 baseline / A3 V3F / A5 V5 三版 TP=28,509 / FP=11,606 / FN=10,938 完全相同，**F1 = 0.7166 三版完全相同**。HCC1395 t30_n20 @ 0.6 purity：B1 / B3 / B5 三版 TP=24,190 / FP=13,487 / FN=15,257 完全相同，**F1 = 0.6273 三版完全相同**。為什麼相同？因果鏈：ClairS-TO snv.vcf 的 FILTER=PASS 集合由 caller 決定，longphase-to phase 改的是 GT、PS、GT2、GT3 欄位 — 也就是 read 怎麼分類到 HP1/HP2 — 不改 FILTER 欄位，所以 PASS set 不變 → TP/FP/FN 不變 → F1 不變。**V5 不改 caller**，不能用 F1 衡量 self-phasing 修補。正確 metric 是 read-level tag concordance（前面 slide 13 看到的 paired GT +13.3 pp）。ΔF1 (0.93→0.6) = -0.0893 是 ClairS-TO 在低 purity 本身偵測下降，與 V5 / baseline 無關。**注意這裡是 turning point — 但 5/9 paired cross-ref 揭露另一面，我們將進入 V5 Layer 1.5 設計缺陷。**」

[ORAL-OPTIONAL] PASS set / FILTER 機制細節 / purity 0.6 N50 微差 -14.5% 解釋 / V3F vs V5 0.6 對比

→ next: slide 15 paired mode 整體無偏移

---

## Slide 15 — paired mode 整體無 priority bug [120 sec / 中 ~360 字]

「paired mode 用不同的 binary — longphase-s 不是 longphase-to fork，**是獨立 codebase**，HP tag 用 HP:Z: 字串編碼（不是 HP:i: 整數）。對 paired tagged BAM chr19 primary reads 582K total / 354,919 tagged 統計分布：HP:Z:2 占 51.6%、HP:Z:1 占 40.5%、HP:Z:2-1 占 4.1%、HP:Z:1-1 占 3.5%、HP:Z:3 只 0.3%。核心 ratio：germline HP1:HP2 = **1:1.275** 接近隨機 — 對比 TO baseline 17.3:1 priority bug 完全沒偏移；somatic HP1-1:HP2-1 = **1:1.169** 也接近隨機。所以 **paired mode 整體沒 priority bug**。再看 chr19 1Mb window som_ratio 統計（57 windows）：mean 0.462 / median 0.494 / stdev 0.332 跨 windows 大變動 — stdev 0.332 跨 windows 大變動代表不同區域有不同 LOH 方向，是真實生物學。三個典型案例：chr19:3M 全 HP2-1 (755/0) LOH 方向特定；chr19:0M 全 HP1-1 (330/1) 反向區域；chr19:17M 對稱 0.500 (265/265) — 這就是 SP1 附近，**paired mode 認雙 sub-clone 共現**，對比 TO baseline 113:0 失衡，paired 看到的才是真實。next 進入 V5 Layer 1.5 設計缺陷的核心揭露。」

[ORAL-OPTIONAL] longphase-s codebase / SomaticHaplotagProcess.cpp:533 HP tag 字串 / paired 軸對齊 vs TO 軸

→ next: slide 16 V5 Layer 1.5 設計缺陷

---

## Slide 16 — V5 Layer 1.5 設計缺陷揭露 [150 sec / 中 ~420 字 ★ 最長 slide]

「這是整份報告最重要的一張 slide。對 paired chr19 read_name × T1.2 baseline / V3F / V5 vote dump 做 JOIN，篩 cnt_HP1+cnt_HP2=0 且 somatic>0 的 events — 也就是 germline-absent 區域 — 共 5,789 chr19 events。看 cross-tab：**baseline hp=11 (HP1 系列) 3,312 / hp=21 (HP2 系列) 791 → 4.19:1 偏 HP1**，這是 priority bug 的次峰偏移。**V3F 全 5,789 events 標 hp=33 somatic ambiguous — 保守不選邊** ✅。**V5 hp=11 = 3,313 / hp=21 = 790 → 4.19:1 與 baseline 完全相同！**機制詮釋：V5 Layer 1.5 設計是 germline=0 時用 somaticHP1 vs somaticHP2 票數決方向，但 self-phasing 機制下 sub-clone somatic 100% 共現 → graph 偏向同一 haplotype → 投票偏向同邊 → Layer 1.5 結果繼承 priority bug 偏移。**所以 V5 Layer 1.5 = priority bug 的 feature 化非修補，把 baseline 用 somatic vote 蓋過 germline 的 buggy 行為，改成 germline 缺席時才用 somatic vote 的 designed 行為，但該區域偏移本質沒變。**V3F 標 hp=33 是更穩健的選擇。caveat：這是 cross-binary axis alignment 問題，paired vs TO 各自獨立 phasing 軸；但 binary-internal 量化（baseline 自身 4.19:1）不受影響。F-paired-D3 follow-up 會量化 V5 改回 V3F 標 hp=33 對 ISM 下游的影響。為什麼 V5 設計者沒測到？因為 V5 設計時未對 paired ground truth 做 germline-absent cross-ref，5/9 paired audit Step D 才補上發現。」

[ORAL-OPTIONAL] cross-binary axis alignment 詳細 / phasing graph 機制連結 §3.2 / Layer 1.5 改回 V3F default 設計選擇 / F-paired-D3 工作量

→ next: slide 17 5 errata

---

## Slide 17 — 5 條 PI errata patched [90 sec / 中 ~330 字]

「PI 報告 4-29 5 條 errata 已 patch；主結論不撤回。**E1** §3.3.3 chr19 SP1/2/3 解讀降級：從「主要 hotspot」改為「可重現案例」（chr19 占全基因組 priority bug 僅 2.16%、rank 19）。**E2** §5.2 V5 working tree commit 狀態：原寫「未 commit」，現已 commit d0bcd8c + 938f0df 在 2026-04-30。**E3** §5.2 priority bug 證據強度升級：原依賴「commit message + 3 IGV 截圖」，補強為「+ 34,855 read-level victims 全基因組鐵證、V3F+V5 修正率 100%」三重佐證。**★ E4** §6.4/§6.5 V5 數值歸因精確化（最重要 errata）：原寫「V5 four-commit chain 整體效益」，精確化為「V5 BAM = Pass 1 only（ploidy bug 讓 purity=0），主要功勞 V3F + Layer 1.5；Pass 2 second round 二次效益尚未獨立量化」。**★ E5** 5/10 加：§5.2 V5 Layer 1.5 設計：原隱含「修補」，改為「germline-absent 區與 baseline 4.19:1 完全相同 — priority bug feature 化非修補；V3F 標 hp=33 反而穩健」。E5 來源是 5/9 paired audit Step D，5/10 amend 進主報告加 errata。E1-E3 是表述精確化、E4+E5 是核心 errata。errata commit chain：f17754f → 2553e96 → 71d21bd 全完成。」

[ORAL-OPTIONAL] 各 errata commit hash / errata patch 工作量 / 為何不撤回 vs 補 banner

→ next: slide 18 整體成熟度 + follow-up

---

## Slide 18 — 整體成熟度 + 5 follow-up [90 sec / 中 ~340 字]

「最後 take-away。整體成熟度 12 維度：10 個 ✅（機制因果、修補設計合理、chr19 SP 對齊、全基因組擴展、V5/V3F zero-sum 釐清、20 指標 0 regression、Caller F1 三版相同、purity 0.6 完整對照、三路徑算法不依賴 purity、Pass 2 量化 +3.51% N50、版本對齊 HEAD 938f0df、Paired Step A+C 整體 NEGATIVE）；1 個 ⚠️（V5 Layer 1.5 germline-absent 設計缺陷 E5）；2 個 ⏸ 待跑（Pass 2 second round 獨立貢獻量化 T1.3、跨樣本擴展 T3）。5 項 follow-up cycle 排序：**F-paired-D3 ★ 1-2 day 最重要** — Layer 1.5 改回 V3F 的 ISM 影響評估，決定 V5 vs V3F 哪個作 default；F-paired-D1 0.5 day — germline-absent 全基因組擴展驗證 4.19:1 是否跨 chr 一致；F-paired-D2 1 day — phase block 內 axis-aligned 分析去除 cross-binary axis caveat；T3 1-2 day — 7 樣本跨樣本 vote audit；T1.3 3 day — 4-cell ablation。**結論：V5 仍可作 production tag baseline；germline-absent 區占 chr19 events ~5% 不阻擋整體；F-paired-D3 ISM 影響量化後是否回歸 V3F default — 待 PI 決策觸發。**」

[ORAL-OPTIONAL] T1.3 4-cell ablation 設計 / 7 樣本 binary patch 工作量 / cnLOH 雙親同源待開放

→ Q&A backup（如 PI 提問則切換）

---

## Q&A Backup Slides

---

## Slide B1 — Pass 2 second round 機制 + N50 +3.51% [120 sec / 中 ~360 字]

「Q: Pass 2 為什麼只跑 2-point 不重跑 3-point？我從 source code 還原 longphase-to phase 真實流程。**不是「2-point vs 3-point 二選一」，而是兩個分別函式各自處理不同任務**：`somaticCalling` 用 3-point edges 加 patternMining first/second/third path，由 `!disableCalling` flag 控制與 purity 無關，Pass 1 內部呼叫一次；`edgeConnectResult` 用 2-point pairwise edges 永遠跑。**低 purity ≤ 0.9**：Pass 1 跑 2-point + 3-point 都跑、Pass 2 跳過。**高 purity > 0.9**：Pass 1 跑 2-point + 3-point 都跑、Pass 2 只重跑 2-point — 不重跑 3-point 因 Pass 1 已產出穩定 origin 分類。Pass 2 incremental 量化（HCC1395 5kHz 0.93）：phased var 1,848,538 → 1,756,339，差 -92,199（-2.90 pp）；blocks 1,808 → 1,631，差 -177（-9.79%）；N50 11,388,114 → 11,788,053，差 +399,939（**+3.51%**）。Pass 2 是 polish/merge：blocks 變少 -10% 但每塊變更長 +3.51%。失去 92K phased var 是 Pass 2 reclassify 為 somatic/./. 的設計目標非 regression。**常見誤解澄清：「低 purity 用 3-point」倒過來；高 purity 才多做事 = Pass 2 多 2-point。**」

[ORAL-OPTIONAL] patternMining first/second/third path 詳細 / Pass 2 不重跑 somaticCalling 原因

---

## Slide B2 — purity 0.6 完整對照 [90 sec / 中 ~320 字]

「Q: 低純度樣本是否變差？purity 0.6 樣本 baseline vs V5 完整對照（HCC1395 t30_n20）。**6 caller F1 全完全相同**：TP / FP / FN / Precision 0.6420 / Recall 0.6132 / F1 0.6273 三版一致（V5 不改 caller）。**9 結構指標**：phased% +4.01 pp（61.82 → 65.83）、n_blocks +18.1%（9,748 → 11,514）、HP:i:33 +20（V5 conservative 標 somatic ambiguous）、AMB% +3.12 pp（合理 conservative）、HP1/HP2 ratio 0.48 → 0.38（修正 tag bias）、purity 計算 0.607 → 0.634（更接近真實 0.6 +0.027）、Pass 2 兩者都不觸發（< 0.9）、LOH regions 完全相同。**唯一 N50 -14.5%（798,903 → 683,296）但仍 ≥ 600 K — P1 設計閾值；blocks 變多但平均長度仍夠**。結論：0 critical regression；V5 conservative tagging 是好事；ploidy bug 在低 purity 自我治癒（polynomial 在 low-input 不會崩到 0；baseline 0.607 接近真實 0.6 不需 d0bcd8c 修正）。」

[ORAL-OPTIONAL] V3F vs V5 0.6 對比 / N50 -14.5% 為何接受 / Pass 2 不針對低 purity

---

## Slide B3 — 7 樣本擴展 + cnLOH 雙親同源 [60 sec / 中 ~250 字]

「Q: 跨樣本一致性？cnLOH 區呢？T3 跨樣本擴展：HCC1395 5kHz 已驗證為本報告主案例；6 樣本待跑 — HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829 — 每樣本 ~3 hr，共 1-2 day（含 binary patch + dump 6 × 3 版本）。cnLOH 雙親同源情境：parent 1 vs parent 2 染色體在 LOH 區因物理上同源難 phase；Layer 1.5 設計選擇對 cnLOH 區影響待量化，與 V5 Layer 1.5 設計缺陷 E5 連動。**Follow-up 排序**：F-paired-D3（1-2 day）改 binary 較直接 actionable，決定 V5 vs V3F default；T3（1-2 day）跨樣本擴展是 generalizability check；T1.3（3 day）4-cell ablation 量化 Pass 2 second round 獨立貢獻。」

[ORAL-OPTIONAL] 7 樣本各自特性 / cnLOH 機制與 LOH 區別 / F-paired-D3 詳細

---

## P5 整體 self-check

| 項目 | 結果 |
|------|------|
| 25 slide × 75-90 sec speaker note | ✅（1+ slide cover 30 sec / 1+ Q&A 60-150 sec 浮動）|
| 字數 ≈ 中 300-360 字 / slide | ✅（Tier 2 budget 滿足）|
| Tier 3 oral-optional 標 | ✅ 每 slide 末已標 |
| Transition (→ next) | ✅ 22 main slide 全標連結 |
| Cliffhanger (slide 14 → slide 15/16) | ✅ |
| 數字一致 vs 03_slide_layout_script.md | ✅（Agent-N grep verified）|

## 時長預估（25 slide × 75-90 sec）

| 段 | slide 數 | 時長 |
|----|---------|------|
| S0 cover + TL;DR | 2 | 2 min |
| S1 觀察起點 (含 04 拆 4a/4b) | 4 | 4.5 min |
| S2 機制 | 3 | 5.5 min |
| S3 量化鐵證 | 2 | 4 min |
| S4 修補設計 | 2 | 4 min |
| S5 驗證 + cliffhanger | 3 | 5.5 min |
| S6 V5 缺陷 | 2 | 4.5 min |
| S7 errata + follow-up | 2 | 3 min |
| **Main 總計** | **22** | **~33 min** |
| Q&A backup | 3 | ~5 min（依問才用）|

⚠ **總時長預估 33 min 超出 target 20 min 約 65%**

## C5 必停討論

22 main slide 33 min 超時。3 個處理選項：

1. **嚴守 20 min**：每 slide 從 ~90 sec 降到 ~55 sec / 中 220 字 → 大幅刪 Tier 2 細節（PI 可能聽不懂機制）
2. **彈性 30 min**：接受 33 min；補 PI 可中斷追問空間
3. **拆兩場**：20 min Part 1 (S0-S5) 修補主線 + 20 min Part 2 (S6-S7+Q&A) V5 缺陷與 follow-up

→ 待 C5 用戶決定後產 build_pptx.py。

