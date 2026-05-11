# Speaker Script — Self-phasing V5 provenance audit (18 slides, 25 min target)

## Timing assumption
- 中文 400 字/min（PLOS Rule 1: ~ 1 min/slide）
- Tier 2 必講：300-360 字/張 = 75-90 sec
- Tier 3 [ORAL-OPTIONAL]：依現場時間決定，可削減

---

## Slide 01: Cover — Main thesis 鎖定
- 字數: 530 中字
- 預估時長: 80 sec (1.3 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (530 > 360)，部分需標 [ORAL-OPTIONAL]

```
開場 30 秒抓 main thesis：本週最關鍵的發現是 5/05 provenance audit 揭露4/29 PI 報告的全部 V5 數值都源自 4/12 產的 BAM，那份 BAM 因 ploidy bug 讓 purity=0、Pass 2 從未觸發。所以「V5 為當前最佳 phasing baseline」的主結論需暫停為「Pass 1 only 條件下觀察」，Pass 2 完整重驗是本週最高優先級 (P0)。同時 5/01 三條路 audit + force_path2only ablation 實證「路 3 second round 抵消路 2 反轉」，新 baseline self-phasing 反而更甚。本次彙報聚焦兩件 critical findings 與下一步決策。

[ORAL-OPTIONAL] 若教授問 deliverable 數量，可補：本週共 9 份 deliverable，包含 4/28 V5 audit、4/29 技術報告、4/30 V3F ablation × 2、5/01 三條路 audit、5/01 force_path2only ablation、5/05 provenance audit。
```

## Slide 02: Background — Thread D 切換 + 本週主線
- 字數: 567 中字
- 預估時長: 85 sec (1.4 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (567 > 360)，部分需標 [ORAL-OPTIONAL]

```
本週上下文先架構好。上週週報明確設定三件事：Thread D 主軸切換、Thread B 撤回、V5 audit chain 推進。本週前 5 天大致按計畫完成，9 份 deliverable 涵蓋 V5 audit / 技術報告 / V3F ablation / 三條路 audit / force_path2only ablation。但 5/05 做provenance audit 時發現一個 critical issue：4/29 PI 報告引用的 V5 數據實際來自4/12 那份 BAM，那份 BAM 是 Pass 1 only 結果，所以「V5 為當前最佳」這個結論的現實基礎需要重新確認。

敘事框架是混合主線：主線是 problem (V5 caveat)，但 sub-thread 是 progress (9 份 deliverable)，所以結尾仍能 actionable。教授視角優先序我建議先聽完 critical finding 與 Pass 2 重驗計畫，再回看其他進度。

[ORAL-OPTIONAL] 若教授追問為何 4/29 報告未及時發現 Pass 1 only，可解釋：當時log 只看「flag 是否觸發」沒看「purity 是否 > 0.9」，這是 5/05 audit 才補上的盲點。
```

## Slide 03: V5 修補鏈 — 5 commits 含 4/30 兩 commit
- 字數: 565 中字
- 預估時長: 85 sec (1.4 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (565 > 360)，部分需標 [ORAL-OPTIONAL]

```
V5 不是 4 commits 是 5 commits，這是本週第一個校正點。原本記錄的 V5 是前三個 commit (PON-only flag、tag layer fix、somatic fallback)，4/30 那天補了兩個關鍵 commit：d0bcd8c 修 ploidy bug，938f0df 從上游 cherry-pick threshold 0.95→0.9 調整。這兩個都是紅色節點，因為 d0bcd8c 解開 Pass 2 觸發機制，938f0df 引入路 3 抵消效應。

之前 audit 留有一個 caveat R1：working tree 有未 commit 改動。5/05 確認 git diff --stat HEAD 為空，全部邏輯都已 commit、可追溯。這個 caveat 隨 938f0df 一起被解決，是本週副效益。

[ORAL-OPTIONAL] threshold cherry-pick 細節：上游 zhenyu 把 highPurity 觸發閾值從 0.95 降到 0.9，讓更多 case 進入 Pass 2 second round。這個改動在我們 case 是雙刃劍 — 後面 slide 6 會看到路 3 second round 反而抵消路 2 反轉效果。
```

## Slide 04: ★ 5/05 Critical — PI 報告 = Pass 1 only
- 字數: 744 中字
- 預估時長: 112 sec (1.9 min)
- Tier 2 必講範圍判定: ❌ 大幅超出 (744 > 700)，必須拆 Tier 3

```
這張是本週最重要的 slide。5/05 做 provenance audit 時發現一個被忽視的盲點：4/29 PI 報告所有 V5 數值都引用 4/12 產的那份 BAM，但那份 BAM 因為 ploidy bug 讓 purity 算出 0，整個 highPurity 條件 false，Pass 2 second round 從未觸發。所以那些好看的數字 — clean PS +13.3pp、AMB% 從 17.5% 降到 8%、HP:i:33 完全消失 — 全部是 「PON-only Pass 1 + tag layer fix」的功勞，不是 V5 完整版的功勞。

Evidence 兩條：A 是 BAM provenance，log 顯示 purity:0、無 second round 字串；B 是機制因果，purity=0 → highPurity=false → Pass 2 跳過。兩條交叉確認。

為什麼這很 critical？因為「V5 為當前最佳 phasing baseline」這個主結論在現有 audit chain（4/29 技術報告、4/28 audit、Memory project_v5_somatic_fallback_verification）裡到處被引用，但實際上這個結論的現實基礎是 Pass 1 only 條件。我們需要 Pass 2 真實觸發版重驗才能維持。

[ORAL-OPTIONAL] highPurity 三次多項式公式在 PhasingProcess.cpp:142-220，purity 從 read coverage / depth 算出來，ploidy bug 讓計算 chain 在 ploidy=0 時短路。
```

## Slide 05: 機制 — ploidy bug 因果鏈
- 字數: 582 中字
- 預估時長: 87 sec (1.5 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (582 > 360)，部分需標 [ORAL-OPTIONAL]

```
把 critical finding 視覺化成 5 階段因果鏈。從左到右：4/12 BAM → ploidy=0 → purity=0 → highPurity=false → Pass 1 only result。中間三步 (紅色) 是 ploidy bug 的影響鏈，右邊的橘色節點是最終結果，也就是 4/29 PI 報告引用的數值。

4/30 d0bcd8c 解了這個 bug，purity 計算重新正常，現在 Pass 2 second round 會真的觸發。但問題是 4/29 PI 報告寫的時候用的還是 4/12 那份 BAM，所以數字維持原樣。要更新就必須在 4/30 修補後 binary 下整套重跑：HCC1395 5kHz sanity、concordance、ISM benchmark，估計 ~25 hr 含 7 樣本平行。

Decisive next step P1 就是這個重跑。如果不做，整條 V5 audit chain 都卡在 Pass 1 only 條件。

[ORAL-OPTIONAL] ploidy 計算分支在 PhasingProcess.cpp:142-220 那個區段，d0bcd8c 改了 5 行 — 詳細 diff 在 audit_01.md §2.3。可以追加討論為何上游沒一起 cherry-pick。
```

## Slide 06: ★ 5/01 三條路 audit — 路 3 抵消路 2
- 字數: 683 中字
- 預估時長: 102 sec (1.7 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (683 > 360)，部分需標 [ORAL-OPTIONAL]

```
本週第二個 critical finding：5/01 三條路 audit。我們先把 PhasingProcess 拆成三條路徑來理解。路 1 是原始 baseline，沒有任何 flag，single pass，self-phasing 17.3:1 殘留。路 2 是 V5 flag PON-only，4 步流程含 somaticCalling 重跑，這是反轉 self-phasing 的核心 — 舊 V5 在 ploidy bug 條件下走純路 2，HP1:HP2 = 0.735 ✅ 反轉。路 3 是新加的 highPurity > 0.9 second round，只重 phase 不重 call，self-phasing 偏移殘留。

問題是新 binary 下 V5 flag 會走「路 2 + 路 3」，路 3 的 second round 把路 2 反轉的東西重新偏掉，最終 HP1:HP2 = 1.400，比新 baseline (路 1+3) 1.400 還高、比舊 baseline (路 1) 1.328 更糟。

教授視角的 implication：V5 flag 在新 binary 下「機制成立」(路 2 仍能反轉) 但「最終效果消失」(路 3 抵消)。這直接影響「V5 flag 還要不要留」的決策，後面 slide 11 會回到這個問題。

[ORAL-OPTIONAL] HaplotagProcess.cpp:484-563 是 somaticCalling 重跑的程式碼，如果想看具體哪幾行讓路 2 反轉，可以打開該段。
```

## Slide 07: 5 版本對比表 — NEW noPath3 ≈ OLD V5
- 字數: 690 中字
- 預估時長: 104 sec (1.7 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (690 > 360)，部分需標 [ORAL-OPTIONAL]

```
這張表是本週最重要的數值對照。5 個版本：OLD baseline / OLD V5 / NEW baseline / NEW V5 / NEW noPath3。前兩個是 4/12 BAM 結果（ploidy bug 條件），後三個是 4/30 修補後 binary 結果。

綠色行：成功反轉 (HP1:HP2 < 1)。OLD V5 (0.735) + NEW noPath3 (1.127) 都成功。紅色行：抵消或更糟。NEW baseline (1.400) + NEW V5 (1.400) 都失敗。

關鍵數字 — HP_33。OLD V5 與 NEW noPath3 的 HP_33 都是 14,524，完全相同。這個等價性證明 NEW noPath3 (強制 highPurity=false 跳過路 3) 復現了舊 V5 的反轉行為。這個假說 PASS — 路 2 仍然成立，問題只在路 3 抵消。

教授視角的 implication：如果我們本地維護 longphase-to-noPath3 變體 binary，等於回到舊 V5 的反轉效果但無 ploidy bug。這是 slide 11 三選項中的 (b)。

[ORAL-OPTIONAL] HP1:HP2 = 1.127 vs 0.735 的差距：兩者都反轉但程度不同。在 noPath3 下因為 baseline 已從 1.328 升到 1.400，反轉後到 1.127；舊 V5 是從 1.328 反轉到 0.735。差距是 baseline 偏移量造成的，不是 V5 flag 本身效果差異。
```

## Slide 08: force_path2only ablation — 假說 PASS
- 字數: 652 中字
- 預估時長: 98 sec (1.6 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (652 > 360)，部分需標 [ORAL-OPTIONAL]

```
這張是負控制實驗 (negative control)。我們想驗證「路 3 抵消路 2 反轉」這個假說。做法：直接在 PhasingProcess.cpp patch 一行，bool highPurity = false，強制讓 second round 條件失敗、跳過路 3，不論實際 purity 計算結果。

如果路 3 真的是抵消來源，跳過路 3 應該復現舊 V5 反轉。結果：HP1:HP2 = 1.127 反轉成功，且 HP_33 = 14,524 與舊 V5 完全相同。等價性 + 反轉成功兩條都符合預期，假說 PASS。

副產物：phase 時間從 2881s 降到 813s，節省 71%。這個沒有實際生產意義 (因為 noPath3 binary 並非 default)，但顯示路 3 是計算密集區段。

教授視角的 implication：這個實驗證明「機制完全清楚」— V5 flag 本身沒壞，問題只在路 3。如果決定維護本地 longphase-to-noPath3 變體 binary，22.55 MB 可直接 deploy。

[ORAL-OPTIONAL] 為什麼選 force_path2only 這種 hack 而非改 threshold？因為 threshold 改動(0.95 vs 0.9) 是上游 zhenyu 那邊維護的，本地改動會 diverge；強制 highPurity=false 是一行 patch，可作為 ablation 工具，不影響上游。
```

## Slide 09: Caller F1 全 6 版本相同
- 字數: 619 中字
- 預估時長: 93 sec (1.5 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (619 > 360)，部分需標 [ORAL-OPTIONAL]

```
這張是 supporting evidence — 確認 caller F1 完全不受 phasing 改動影響。0.93 pass quality 6 個版本都是 0.7166，0.6 pass quality 3 個版本都是 0.6273。0.6273 這個數字是本週首次發布 (4/30 V3F purity ablation)，之前只有 0.93 點數據。

為什麼 caller F1 對 phasing 改動完全不敏感？因為 longphase-to 只動 phasing graph (decide HP tag 怎麼分配)，不動 FILTER 欄位 (PASS/LowQual/...)。Caller F1 看的是 FILTER，所以 phasing 改動 invisible。

教授視角 implication：4/29 PI 報告寫「caller F1 不變」是對的，但要附註 — 這個不變不代表「phasing 沒改動」，只代表「caller 沒受傷」。phasing 品質要用 ISM SuggestFilter 下游 F1 衡量(下游 F1 才會看到 HP tag 分配差異)。

[ORAL-OPTIONAL] 0.6273 這個數字之前未公布的原因：4/30 V3F ablation 第一次跑出 0.6 pass quality 點，作為 cross-purity 確認。詳見該 ablation 報告。
```

## Slide 10: ⚠ V5 結論 caveat banner
- 字數: 615 中字
- 預估時長: 92 sec (1.5 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (615 > 360)，部分需標 [ORAL-OPTIONAL]

```
這張是把所有需要更新的 artifacts 列出來，作為 action item 的視覺整合。不是新發現，是把 critical finding 對 audit chain 的影響具體化。

Banner 內文：「V5 為當前最佳 phasing baseline」這個結論基於 Pass 1 only 條件，目前需暫停為「Pass 1 only 條件下觀察」。等 Pass 2 重驗 ~25 hr 完成才能恢復為「驗證後結論」。

受影響 artifacts 4 個：4/29 longphase TO vs V5 技術報告、4/28 V5 audit、Memory project_v5_somatic_fallback_verification (status 改 needs_rerun)、Thread D 主軸切換文件 (因為它引用 V5 為 baseline)。每個都要補 caveat 或更新狀態。

教授視角的 implication：要決定是「立即補 caveat」還是「等 Pass 2 重驗完一次更新」。兩個選擇 trade-off 不一樣 — 立即補 caveat 透明但會讓上游讀者困惑；等重驗完一次更新清晰但需 1 週透明度落差。這是 Top Asks #1。

[ORAL-OPTIONAL] 也可以考慮折衷：先在 Memory 改 needs_rerun (內部記錄)，等重驗完才更新外部技術報告。
```

## Slide 11: Decision tree — V5 留存三選項
- 字數: 531 中字
- 預估時長: 80 sec (1.3 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (531 > 360)，部分需標 [ORAL-OPTIONAL]

```
這是 Top Asks 的第二個必判斷決策點。新 binary 下 V5 flag 不再有效反轉 self-phasing，我們有三個選項。

(a) 接受新 default：什麼都不做。優點是與上游同步無維護成本，缺點是下游 ISM 必須容忍 self-phasing 1.40 殘留。0 hr。

(b) 維護 longphase-to-noPath3 變體 binary：22.55 MB 已建好。可復現舊 V5 反轉效果。缺點是上游每次更新需週期性同步，~1 hr/week。

(c) 上游 PR：提案 zhenyu 修路 3 邏輯讓 second round 重 call。最徹底但時程不可控，~2 週。

我的建議：短期 (b) 作為 ISM 重驗備援，中期推 (c)，(a) 作 fallback。但這要教授判斷 — 因為 (c) 涉及上游溝通的政治成本，(b) 涉及 binary divergence 風險。

[ORAL-OPTIONAL] 22.55 MB binary 大小其實很合理，longphase 主執行檔本來就 ~22 MB。已 staged 在本地 build/，未 commit；commit 與否視 (b)(c) 決策。
```

## Slide 12: 跨樣本 938f0df 影響 [U]
- 字數: 506 中字
- 預估時長: 76 sec (1.3 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (506 > 360)，部分需標 [ORAL-OPTIONAL]

```
這張是 risk surface — 把 cross-sample 不確定性視覺化。目前 HCC1395 5kHz 一個樣本的 baseline HP1:HP2 從 1.328 (OLD) 升到 1.400 (NEW)，+5.4%。其他 6 個樣本還沒測。

為什麼 cross-sample 重要？因為「新 baseline 比舊 baseline 更偏 HP1」可能是 systematic (每個樣本都受 938f0df 影響) 或樣本特異 (只是 HCC1395 5kHz 的 noise)。若 systematic，slide 11 的決策 (b) 或 (c) 都需要走；若樣本特異，可暫接受新 default。

Priority 2 預估 4 hr (跨 6 樣本平行)，是 Pass 2 重驗 (P1, 25 hr) 的子集。可以 P1+P2 一起做。

[ORAL-OPTIONAL] HCC1395 DORADO 是 Dorado basecaller 重 call 過的版本，與 5kHz 是同一樣本不同 basecaller。這兩個算 1.5 個樣本 — 若兩者一致，仍需要其他樣本確認。
```

## Slide 13: HPFineNGroups 在新 baseline 重驗
- 字數: 608 中字
- 預估時長: 91 sec (1.5 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (608 > 360)，部分需標 [ORAL-OPTIONAL]

```
這張是 audit chain 的另一條延伸 — HPFineNGroups subclone marker。4/23 週報把它從 ⭐4 降到 ⭐3，重詮釋為 phasing signature。但所有 audit 都是 4/12 BAM (Pass 1 only)。新 baseline 下需要重驗。

需要驗的三件事：(1) NG≥3 = 0 是否仍成立 (若是 → phasing signature 假說 ✓)；(2) flag=on (NEW V5 路 2+3) 是否變化 (預期類似因為仍走路 3)；(3) noPath3 binary 是否復現舊行為 (若是 → 完全對稱 paired check)。

Priority 3 估時 0.5 hr，與 R1 (Pass 2 重驗) 共用 binary，所以可一起做。

教授視角 implication：4/28 audit 結論「marker tier ⭐3」目前暫定，等 R3 驗完才能確認。若 R3 結果與 4/28 不一致，可能要再降一級或重新評估機制。

[ORAL-OPTIONAL] HPFineNGroups 是 src/include/HPFineNGroups.hpp 定義的特徵，計算 HP tag 內部精細分群的群數。詳細定義已在 Memory feedback_feature_name_vs_definition_rule.
```

## Slide 14: Thread B 撤回背景
- 字數: 701 中字
- 預估時長: 105 sec (1.8 min)
- Tier 2 必講範圍判定: ❌ 大幅超出 (701 > 700)，必須拆 Tier 3

```
這張是把 Thread B 的演進與目前狀態講清楚，因為它與本週主軸 (Thread D 切換) 直接相關。

Thread B 原本是 4/19 提的假說 — NG=2 LOH-AF interaction 產生 methylation 訊號。4/23-26 重新審視時發現 self-phasing 會先污染 HP_Ratio LOH，所以「methylation 訊號」其實是phasing artifact。5/01 機制 pivot 為 phasing signature。

TO 層：撤回 methylation 機制，重詮釋為 phasing signature，等 Pass 2 重驗。
Paired 層：保留 POSITIVE 但加 caveat。Inter AF→NGroups +0.705 ~ +0.787 (7/7) 是 paired 層的觀察，paired 不受 self-phasing 影響所以保留，但需要獨立的 phasing-vs-methylation 驗證才能確認 mechanism。

教授視角 implication：Thread B 不是完全失敗 — paired 層仍有訊號，TO 層需重驗。這也說明 self-phasing artifact 對下游分析的污染範圍很廣。

[ORAL-OPTIONAL] Memory 兩條：project_loh_subclone_af_methylation_positive (paired 保留) 與 project_hpfinengroups_subclone_marker (TO 降級)。詳情可追問。
```

## Slide 15: Future Priorities P1-P5
- 字數: 467 中字
- 預估時長: 70 sec (1.2 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (467 > 360)，部分需標 [ORAL-OPTIONAL]

```
把 Future 5 件事視覺化成優先序卡片。P1 Pass 2 重驗是決定性 25 hr — 不做的話「V5 為當前最佳」結論失去現實基礎，所有引用 V5 的 audit chain 都卡住。P2 跨樣本 4 hr 與 P1 共用 binary 並行。P3 PI 報告 caveat 2 hr 看 P1 完成度決定立即補還是等重驗完一起。P4 V5 留存決策文件 1 hr 等教授判斷三選項後再寫。P5 Memory 校正 0.5 hr 收尾。

本週執行計畫：P1+P2 並行 ~25 hr wall time，週日到週四完成。P3 視 P1 結果決定。

預期交付：Pass 2 重驗報告、cross-sample 938f0df 評估、caveat 補丁。下週週報就能基於 Pass 2 真實數據更新 V5 結論。

[ORAL-OPTIONAL] 25 hr 不是純計算時間，是 wall time。實際運算 7 樣本平行 ~21 hr，分析整理 ~4 hr。如果 cluster busy 可能延到 30 hr，但仍週內。
```

## Slide 16: Take-home 3 件事
- 字數: 594 中字
- 預估時長: 89 sec (1.5 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (594 > 360)，部分需標 [ORAL-OPTIONAL]

```
Take-home 3 件事，對應 main thesis 的三個面向。第一件事：commit 鏈完整化。V5 不是 4 commits 是 5 commits，4/30 補了兩個關鍵 commit (d0bcd8c ploidy fix + 938f0df threshold cherry-pick)。working tree caveat R1 同時解決。這部分是 ✅ 已完成。

第二件事：Pass 1 only caveat。4/29 PI 報告所有 V5 數值都基於 Pass 1 only 條件，因為 4/12 BAM 的 ploidy bug 讓 Pass 2 從未觸發。這部分是 ⚠ 主結論暫停為「Pass 1 only 條件下觀察」。

第三件事：Pass 2 重驗 P0。在 4/30 修補後 binary 下重跑 sanity/concordance/ISM benchmark，~25 hr 含 7 樣本平行。本週執行，下週週報就能恢復為「驗證後結論」。

三件事按 main thesis 順序：commit chain (5 commits) → Pass 1 only → Pass 2 重驗 P0。

[ORAL-OPTIONAL] 強調本週的 critical finding 不是 V5 失敗，是「V5 結論基礎需校正」— 更謹慎的態度而非結論翻轉。
```

## Slide 17: Q&A 預備 7 個
- 字數: 340 中字
- 預估時長: 51 sec (0.8 min)
- Tier 2 必講範圍判定: ✅ 在 75-90 sec target 內

```
Q&A 預備頁是 backup — 不一定講，看時間與教授提問方向。必問 3 個是直接挑戰主結論的：PI 結論還能用嗎、新 baseline 偏移 systematic 嗎、V5 flag 還要留嗎。三個都已在前面 slide 處理過，這裡只是 quick-reference。

可能問 4 個是技術細節：caller F1 為何重要、4→5 commits 怎麼回事、HPFineNGroups 還有效嗎、下週時間夠嗎。每個都有 1 句 evidence-backed 回答。

[ORAL-OPTIONAL] 如果教授問其他不在這 7 個的問題，預設回應：「我把 evidence 整理到 audit_01.md 之後 follow up」— 不要硬答超出已驗證範圍的東西。
```

## Slide 18: Acknowledgments + References
- 字數: 553 中字
- 預估時長: 83 sec (1.4 min)
- Tier 2 必講範圍判定: ⚠ 略超出 (553 > 360)，部分需標 [ORAL-OPTIONAL]

```
結尾頁面 — Acknowledgments + References。本週 9 份 deliverable 完成，包括 4/26 Thread D 主軸切換、4/28 V5 audit、4/29 兩份 (技術報告 + supplement)、4/30 兩份 V3F ablation、5/01 兩份 (三條路 + force_path2only)、5/05 critical finding provenance audit。Memory 同步校正一條。

上游 zhenyu collaboration 兩個 commit：938f0df threshold cherry-pick 與 d0bcd8c ploidy fix。兩個都是 4/30 同一天落地，讓本週 audit 正好卡在 caveat 邊界。

Reference 10 項涵蓋本週與上週相關文件，完整 master_draft.md 在 InterSubMod/docs/reports/validated/...

[ORAL-OPTIONAL] 上週 4/23 週報還有 3 份 (Thread B 撤回、HPFineNGroups 降級、LOH-AF 保留)，與本週 9 份合計 ~12 份。Memory 與證據鏈索引可以追問。
```

---

## 總計時長估算
- 總字數: 10547 中字
- 估算時長: 26.4 min (1582 sec)
- Target: 25 min
- Delta: +1.4 min

**⚠ 超時警告**：建議把 Tier 3 [ORAL-OPTIONAL] 段落主動削減；以下是建議拆遷清單
