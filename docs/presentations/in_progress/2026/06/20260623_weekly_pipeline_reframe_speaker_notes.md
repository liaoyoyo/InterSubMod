<!--
title: 週報口頭講稿（reframe 新版）— 位點分類方法 pipeline（6/15–6/23）
date: 2026-06-23
status: 定稿（reframe 新版，三輪 CP-D 雙審 0 BLOCKING）
deck: 20260623_weekly_pipeline_reframe_lab_deck.html
build_branch: research/subclonal-reconstruction-202606
verified_by: wc1v407en（5 主線查驗）+ w9456xpig/ws9yqa9pn/wqo60wnyi（三輪 CP-D 雙審）
-->

# 週報口頭講稿 · 位點分類方法 pipeline（reframe 新版）

> **總時長**：約 14 min + 問答。**敘事**：流程型（方法學說明 → 驗證 gate → 分類 → 判別機制 → chr2 解釋 → 誠實邊界）。
> **全程不鬆口的紅線**：① 方法是 **SOUND_WITH_CAVEATS（非 ad-hoc），判別力（分 subclone）尚未測**，不是 fully validated。② 對齊的位點**多為 cis-ASM**，subclone 是**候選需確認**。③ 單樣本 HCC1395 ⭐2-3 PARTIAL。
> **念法總則**：reclassify 念「Python pilot 46.6%、report_claim」不講 committed verified；對齊念「~18%、多為 cis-ASM」不當核心發現；chr2 絕不講「5 個 subclone」。

---

## S1 · 封面（⏱ 0.5 min）

🎤 「教授好。這禮拜我把研究的呈現方式換了個角度——不是講『方法揭示了什麼侷限』，而是**『我建立了一套位點分類的方法 pipeline，能把每個位點清楚分類、並用來解釋 read 跟可能的 clone 狀況』**。」

🔴 念法：標題下那行一定要念出 caveat——「方法是 **SOUND_WITH_CAVEATS**，**判別力還要等 normal cis-control**」。

➡ 過場：「我先用一頁講這套方法做了哪四件事。」

---

## S2 · 本週一句話 BLUF（⏱ 1.5 min）

🎤 「四步 pipeline：① 建甲基距離樹 ② 驗證方法學健全 ③ 分類每個位點 ④ 交叉特徵解釋。」

🎤（指底部誠實定調，**必念**）「先把邊界畫清楚：**方法是 SOUND、非 ad-hoc，但判別力（能不能分 subclone）還沒測，要 normal cis-control**。對齊的位點**多為 cis-ASM**；subclone 是**候選、需外部確認**。**方法的價值在分類、解釋、把判別問題變 well-posed，不在『找到 subclone』。**」

🔴 紅線：先講邊界，後面講分類就不會被當 overclaim。

➡ 過場：「先講方法怎麼運作。」

---

## S3 · 方法學說明 · 流程圖（⏱ 2 min）

🎤（指流程圖左到右）「四步：**① 取材**——每個 sSNV 取周圍 ±5kb 的 reads；**② 建甲基距離樹**——read 跟 read 算 BERNOULLI 距離、UPGMA 成樹；**③ 切群**——用顯著性閘，嚴閘不易誤切、寬閘多撈候選；**④ 分類**——每個位點走決策流落進 5 個桶之一。」

🎤（指中間綠色 gate）「**最重要的是這個關卡：四步做完之後，方法必須先通過驗證，下面的分類才可信**——這就是下一頁。」

➡ 過場：「先講驗證。」

---

## S4 · 驗證 gate（⏱ 1.5 min）· 你指定的順序核心

🎤 「在我給分類數字之前，先證明這個方法 SOUND。三項：」
- 「**驗證 1**：read×read 距離→分群→PERMANOVA 是成熟框架（Legendre/Anderson），不是我亂湊的。」
- 「**驗證 2**：換三閘、5 桶、純幾何切，結論都一樣——但我要誠實講，這是**對切法不敏感（robustness），不是獨立證據**，因為它們用的是同一份距離矩陣；真正跨樣本的獨立性是那個全 7 樣本的 Jaccard。」
- 「**驗證 3**：舊判定漏報了顯著結構，reclassify 修了。**這裡數字要小心**——救回 46.6% 是 **Python pilot、report_claim、還沒落 C++**，不是已驗證的 committed 結果；C++ 全基因組版 53.9% 還在另一個 branch 待落檔。」

🔴 紅線：「這三項過了，下一頁的分類**（characterization）才站得住**——但判別力（分 subclone）還是沒測，那要 normal cis-control。」

➡ 過場：「好，方法 SOUND 了，來看分類結果。」

---

## S5 · 分類結果 5桶（⏱ 2 min）· 核心數字頁

🎤 「方法 SOUND 之後，把全基因組 34,736 個位點清楚分成 5 類，每個位點落唯一一桶。」

🔢 念法（**順序很重要，先講多數**）：
- 「**最大的一塊是 S6 乾淨單群、65%——三分之二的位點其實沒有結構。**」
- 「**有結構且對齊的 S1 只有約 18%**（5,538 個）。」
- 「有結構但不對齊 S2 6.24%；不穩定、次閾值各約 4-6%。」

🔴 紅線（**對齊數字一定要配 cis-ASM**）：「但這個 **18% 對齊，多數是 cis-ASM、不是 subclone**——證據是下一頁會講的：FP（germline 假陽性位點）的對齊率 35.65% **還比** TP 的 24.66% 高。所以對齊是 cis-ASM 指紋，不是 subclone。」

🔴 補充：「S4 是 chance-level 的候選、S5 殘留是 0 所以跳號、這批 5 桶數字在另一個 worktree（我標了來源）。」

➡ 過場：「那『對齊』是怎麼判的？」

---

## S6 · 判別機制（⏱ 1.5 min）

🎤（指分岔）「切群只給『有沒有結構』。要判它是什麼，得跟**獨立的標籤**交叉——HP、等位基因這種**事先就知道、跟甲基無關**的標籤。」

🎤 「對齊 germline → 就是 **cis-ASM**，等位效應、多數情況、不是 subclone；不對齊、而且其他特徵（AF/CN/LOH）支持 → 才是 **subclone 候選**。」

🔴 紅線（**兩個關鍵**）：「① 用獨立標籤判 = **非循環**（不是 double-dip）；② 配特徵 AF/CN/LOH **不是用來『找』subclone，是用來排除 cis-ASM、給候選一個定義**。而且**候選不等於已確認，要 cis-control**。」

➡ 過場：「拿一個真實位點走完整流程看看。」

---

## S7 · chr2:18,086,020 例子（⏱ 2 min）

🎤（指熱圖）「拿 chr2:18M 走完整 pipeline，看甲基訊號**能不能被解釋**。熱圖每一格是一個 CpG 在兩套 basecaller（HKU/DORADO）× α/β 兩條 haplotype 的甲基值，紅色高、藍色低。」

🎤（指綠框）「重點：**只有 {3.1, 3.2} 這兩個 CpG 是跨 basecaller 最強複製的乾淨候選**；而 **{3.3, 3.4, 3.5, 4.1} 被既存的 germline-ASM confound——這些不能當 subclone 的新生甲基**。」

🔴 紅線（**絕不鬆口**）：「誠實講——這個位點的 **locus-subclone 還沒證實**；它落在 SEQC2 的 high-confidence 空隙、**無法用外部真值評估**；整體最多支持 **約 3 條 regional states、不是 5 個 subclone**；證據天花板是 **L2**。⛔ 我**絕不會講『找到 5 個 subclone』**。」

➡ 過場：「所以這套方法的誠實邊界在哪。」

---

## S8 · 誠實邊界（⏱ 1.5 min）

🎤 「核心一句：**分類不等於確認**。」

🎤 「為什麼對齊多是 cis-ASM——**FP（germline 假陽性位點）的 confident 結構率 35.65%，比 TP 的 24.66% 還高**。FP 是 germline，它結構更多，正好說明這些結構是等位效應、不是 tumor subclone。**這是 tumor-only 下分不開的根本原因。**」

🔴 紅線：「所以 **normal cis-control 是『判別力到底成不成立』的決定性測試、目前還沒做**。這套分類是 **characterization 觀察層，不是 filter、不是自動 subclone caller**；單樣本 ⭐2-3、沒有 FDR 校正。」

➡ 過場：「總結。」

---

## S9 · 收束（⏱ 1.5 min）

🎤 「總結這週的正面貢獻：**建立了一套方法學健全的位點分類 pipeline、一個交叉標籤+特徵的解釋架構，並且把『需要 normal 才能判別』這件事定位成一個 well-posed 的明確測試。**」

🎤（指下一步）「下一步按優先序：**normal cis-control 最優先**——這不只是下一步，是判別力成不成立的決定性測試、現在還沒做；然後 chr2 的多 sSNV CCF 做非循環確認；single-cell 當最終裁判。COLO829 補甲基 normal 是**解鎖跨樣本之後才談 ⭐4，現在還不是**。」

🎤 「標題我會界定成 **regional LOH-constrained partition、不是完整的 tree**；骨幹是 somatic haplotag、甲基是 corroborate、不是 driver。」

➡ 「以上是這週進度，開放討論。」

---

## 教授問答預備（機動）

**心法：任何『你是不是找到 subclone 了』的問題，一律回到「我們建立了能 characterize、能分類的方法，subclone 是候選、確認需要 normal cis-control / single-cell」。**

**Q · 方法既然 SOUND，是不是就能找 subclone 了？**
> SOUND 只證明方法非 ad-hoc、流程合理。能不能判別 subclone 是另一回事——那要 normal cis-control 把 cis-ASM 扣掉，現在還沒做。

**Q · 對齊 ~18% 是不是就是 subclone？**
> 不是。對齊多數是 cis-ASM；FP（germline）的對齊率還比 TP 高，證明對齊是等位效應指紋、不是 subclone 檢定。

**Q · chr2 那個位點到底有沒有 subclone？**
> 有 somatic 事件、甲基能分出約 3 條 regional states，但 locus-subclone 未證實、落 SEQC2 空隙無法外部驗證、天花板 L2。是 proof-of-concept，不是確認。

**Q · 三種方法都同意，不是很強的證據嗎？**
> 它們共享同一份距離矩陣，是 robustness（對切法不敏感）、不是統計獨立。真正獨立的是全 7 樣本的 Jaccard。

---

> **驗證足跡**：deck 數字經 wc1v407en（5 主線查驗）+ 三輪 CP-D 雙審（w9456xpig/ws9yqa9pn/wqo60wnyi，末輪 0 BLOCKING、三紅線 consistent）。reclassify 46.6% = Python pilot report_claim；5 桶數字 @ism-review-infra worktree（待落檔）。
