<!--
title: 週報口頭講稿 — Subclonal Reconstruction 判別方法論建立（6/15–6/22）
date: 2026-06-22
status: in_progress（CP-D 修訂版，待 CP-E 定稿）
deck: 20260622_weekly_0615_0622_subclone_methodology_lab_deck.html
build_branch: research/subclonal-reconstruction-202606
verified_by: wc1v407en（5 主線查驗）+ w9456xpig（CP-D 雙審）
data_sources: 同 deck frontmatter
-->

# 週報口頭講稿 · Subclonal Reconstruction 判別方法論建立

> **總時長**：本體約 14–15 min + 問答。**敘事弧**：進展型 + 辨別力雙核心（M2 a-priori 判別軸 + M3 切群≠subclone）。
> **三條全程不鬆口的紅線**：① proof-of-concept ⭐3（單樣本 PARTIAL）② 甲基 = corroborate 非 driver ③ 確認 subclone 需外部 ground truth。
> **念法總則**：報「初步觀察 / 已驗證」分清楚；reclassify 念 46.6%（committed），53.9% 標「待落檔」；82% 念「CN-獨立 8.1–17.6%」原文口徑。

---

## S1 · 封面（⏱ 0.5 min）

🎤 「教授好。這禮拜的主軸是 **subclonal reconstruction**。一句話總結：我這週把『**怎麼辨別真假 subclone**』的方法論建成了一條完整的鏈條。先把定位講清楚——**骨幹是 somatic haplotag、甲基化只做佐證（corroborate），現在是 proof-of-concept、單樣本為主**。」

➡ 過場：「我先用一頁講這週到底建了哪四個能力。」

---

## S2 · 本週一句話 BLUF（⏱ 1.5 min）

🎤 「這週建了四個能力：① 修正 ISM 的假陰性；② **確立正確的判別軸**；③ 切群三閘重設計；④ 文獻定位加論文基石。其中**第二跟第三是今天的雙核心**，其他三個是支撐。」

🎤（指底部誠實邊界，**一定要主動念**）「但我先把邊界畫清楚：**甲基化是 corroborate、不是 driver；單一 bulk 樣本只能 characterize，真正確認 subclone 需要外部 ground truth，像 single-cell 或 multi-region。**」

🔴 紅線：先講邊界，後面就不會被質疑 overclaim。

➡ 過場：「為什麼這週要做判別方法？要從上週的接點講起。」

---

## S3 · 背景與接點（⏱ 1.5 min）

🎤 「上週我證明了『**切得出結構**』——無監督 UPGMA 分群、PERMANOVA 配 PERMDISP、做了狀況分類。但這禮拜的問題升級了：**有結構，不等於是 subclone**。」

🔢 念法：「全基因組大約 **91.7% 的 region 落在 Weak/Noise**」——**講『落在 Weak/Noise』，不要講『91.7% 假陰性率』**（那是定義過嚴、不是 bug）。

🎤 「所以這週要回答兩個問題：哪些切出來的群是真 subclone、哪些是 cis-ASM 假象？以及——怎麼用**非循環**的方法去判別。」

➡ 過場：「先講第一個核心：判別軸。」

---

## S4 · 核心 1 · a-priori 判別軸（⏱ 2.5 min）★重點

🎤（指左圖紅色）「先講為什麼**無監督切群不能拿來判別**。左邊這個流程：我拿 read×read 的距離矩陣，先用 silhouette 挑一個『最分得開』的切法，然後**再用同一份距離**去做 PERMANOVA 檢定它顯不顯著。問題在這個迴圈符號——**資料被用了兩次**：先選、再測。這叫 double-dip，**結果一定顯著、但那是假的**。所以無監督的 Noise 在這裡會衝到 97–100%。」

🎤（指右圖綠色）「正解在右邊：**群的定義改用獨立的外部標籤**——haplotag、等位基因、扣掉 normal。群不是從距離挑出來的，PERMANOVA 才合法。」

🔴 紅線（**最重要的一句**）：「但我要很小心地講——a-priori 的優勢是『**合法的 null、免 collider、可解釋**』，**不是『偵測率更高』**。那個『0% noise』是定義使然的 tautology，**不是判別力**；真正非循環的判別率大概只有 **7.6%**。而且 within-HP 看到的甲基分離是 **cis-ASM 候選**，不是已確認的 subclone。」

➡ 過場：「這帶到第二個核心——切群切出來的，到底是不是 subclone。」

---

## S5 · 核心 2 · 切群 ≠ subclone（⏱ 2.5 min）★重點 · 數字最硬

🎤（指中間分岔）「同一個『甲基可以分成兩群』的觀察，其實有**兩種解釋**，判讀規則就是**看 normal 有沒有**。」

🎤（指左橘）「左邊：如果 **normal 也有**，那就是 **cis-ASM**——甲基綁在 germline 的 haplotype 上，是正常的等位效應。**這是多數情況。**」

🎤（指右灰）「右邊：只有 **normal 沒有**，才**可能**是 subclone，甲基綁在 tumor 的演化譜系上。注意我用灰色、寫『候選不等於已確認』。」

🔢 念法（數字最硬，可以有底氣）：
- 「**全 7 個樣本的 Jaccard 都很低，0.091 到 0.161**——意思是無監督切出的群，跟標籤對齊的群幾乎不重疊。**這是這週唯一跨樣本一致的結果。**」
- 「對齊 germline 的群，TP 富集 **3.29 倍**——這是 cis-ASM 的 characterization 訊號。」
- 「而真正的『subclone 候選』、不對齊 germline 的，TP 反而**貧化到 0.65 倍、FP 還更多**——代表這些結構**不是 somatic 特異的**。」
- 「而且這個結論**對方法不敏感**——換三閘、5 桶 null95/90、還是純幾何切，FP（germline）的結構都不少於 TP（5 桶版本 FP 35.65% 還大於 TP 24.66%）。但要誠實講：這些是**同一份距離矩陣的不同切法、不是獨立證據**；真正跨樣本的獨立性是那個全 7 樣本的 Jaccard。」

🎤 收尾：「所以結論很清楚：**tumor-only 設定下，單一 bulk 樣本沒辦法把 cis-ASM 跟 subclone 分開，需要 normal cis-control。**」

➡ 過場：「接下來三頁是支撐——方法本身的健全性、文獻定位、跟論文基石。」

---

## S6 · 支撐 · ISM 假陰性修正（⏱ 1.5 min）

🎤 「ISM 方法本身是 **SOUND** 的——這點外部文獻也支持。問題出在舊的 verdict：它把 91.7% 判成 Weak/Noise，但**其中 Noise 有 74.3%、Weak 有 98% 其實帶著顯著的 PERMANOVA，只是最終判定從來沒去引用它**。我用 reclassify-v2 修了這個。」

🔢 念法：「救回率我引 **46.6% 這個 committed 的版本**（3810/8179）。另外有一個 C++ 全基因組加 level 軸的 53.9%，但**那個還在等落檔，先不當 headline**。兩版的 FP 都是 0——原本的 Strong 一個都沒被降級。」

🔴 紅線（指右下）：「但要講清楚——『**救回**』是把被誤判的『有結構』region 標回 valid，是**結構存在性**，**不是 TP/FP 的判別力提升，這不是一個 filter**。」

➡ 過場：「方法健全之後，下一個問題是——我們在整個領域裡站哪？」

---

## S7 · 支撐 · 方法論定位 + 天花板（⏱ 2 min）· 最易被打臉

🎤（先給前提，指頂部橫幅）「文獻盤點了 74 個來源、0 個真衝突。這頁我講一個本來想做、但發現行不通的事：**我們本來想用『CN 分群和甲基分群對齊』來當 subclone 證據**。」

🎤（指 DAG）「**但這個對齊是假的。** 因為有一個共同上游——『哪一條 allele 存活下來』、也就是 CN/LOH 事件——它同時驅動了 CN 軸和甲基 β。所以兩邊看起來相關，其實是被共同上游製造出來的偽相關。**講白一點：甲基只是 CN 的一面鏡子。**」

🔢 念法（**原文口徑，不要念 82%**）：「外部最強的證據是 Martin-Trujillo——他發現印記 DMR 裡，**真正 CN-獨立的表觀變化只有 8.1 到 17.6%，也就是逾八成是被 CN/LOH 解釋的**。但我一定要補一句限定：**這是限定在 37 個印記 DMR、array 資料，不能外推成『所有 ASM 都是 CN 假象』**——那樣反而會傷到我們自己的 ASM 存在性。」

🎤（指綠盒）「那有沒有乾淨的 anchor？有——**SV 軸**。一條 read 有沒有物理跨過斷點，跟它的甲基值沒有因果關係，所以它是唯一非循環的。」

🔴 紅線：「還有，**停止講『first read-level methylation lineage』這種首創措辭**——2026 年的 PCDH 已經很接近了。」

🎤（指右側天花板）「所以確認的黃金標準還是 single-cell / multi-region / longitudinal，單 bulk 只能 characterize。」

➡ 過場：「最後一個支撐——為什麼這整套需要甲基。」

---

## S8 · 支撐 · 論文基石 motivation（⏱ 1.5 min）

🎤（指圖）「這是論文的 motivation。somatic SNV 在基因組上**太稀疏**了——密度只有每百萬鹼基 10.6 個，相鄰兩個的中位間距 51.4 kb。可是一條 read 中位才 4.9 kb，**大概只有間距的十分之一**（指藍 read 條）。」

🔢 念法：「結果就是——**一條 read 跨到兩個以上 sSNV 的機率只有 0.74%，94.2% 的 read 連一個都跨不到**。」

🎤 「所以這就是為什麼我們**需要甲基化——因為它在每條 read 上處處都有**，可以補上 sSNV 連不起來的訊號。」

🔴 紅線：「但這是**單樣本 HCC1395 的 PARTIAL 結果**；而且我只論『稀疏所以需要補位』，**不會滑到『甲基驅動偵測』或當 filter、建 tree**。」

➡ 過場：「總結這週。」

---

## S9 · 收束 · 誠實定位 + 下一步（⏱ 1.5 min）

🎤 「總結：這週把判別方法論建成完整鏈條，定位是 **proof-of-concept ⭐3，骨幹是 somatic haplotag、甲基是 corroborate**。**正面講——我們建立了一個方法框架：在 bulk ONT 上把 cis-ASM 跟 subclone 候選分類開來，並把『需要 normal 才能判別』定位成一個 well-posed 的明確測試。**」

🎤（指左，**主動畫線**）「我先主動講我們**不宣稱**什麼：不宣稱甲基驅動 genome-wide 重建、不宣稱完整的 phylogeny tree、不宣稱 tumor-only 無監督能偵測 subclone、也不講『首創』。」

🎤（指右下一步）「下一步按 ROI：**最優先是 normal cis-control**——這不只是下一步，是『**判別力到底成不成立**』的決定性測試。因為現在的『不判別』是 tumor-only 下的結論；把 normal 的 cis-ASM 扣掉之後，殘餘的 tumor-specific 結構才是 subclone 的真正候選。然後是 chr2:18M 的多 sSNV CCF 做非循環確認；single-cell 當最終裁判；COLO829 是解鎖跨樣本之後才談 ⭐4，現在還不是。」

🎤 「標題的 reconstruction 我會界定成 **regional LOH-constrained partition，不是 full tree**。」

➡ 過場：「以上是這週的進度，接下來開放討論。」

---

## S10 · 教授問答預備（機動）

**心法：任何『你是不是找到 subclone 了』的問題，一律回到「我們 characterize 了候選，確認需要外部 ground truth」，不鬆口。**

**Q1 · 論文還站得住嗎？**
> 站得住，但要嚴格 frame 成 proof-of-concept。不宣稱完整 tree；骨幹是 somatic haplotag，甲基是 corroborate。賣點是『**在 bulk ONT 上，用 normal-baseline 的 cis-test + 無監督結構檢定，把 cis-ASM 和 subclone 候選乾淨地分類**』這個方法整合。

**Q2 · 真 subclone 訊號還在嗎？**
> 在 a-priori 條件化軸上有候選——within-HP 的 cis-ASM 候選，上界大約 18%。但**這是候選不是確認**；單 bulk 在原理上沒辦法自我確認，要靠外部 ground truth。

**Q3 · 下一步怎麼確認？**
> 三步：normal cis-control 先把 cis-ASM 分離掉；chr2:18M 的多 sSNV CCF 做非循環確認；single-cell 當 terminal arbiter。

**可能的追問 · 為什麼無監督不行？**
> double-dip——同一份距離矩陣先用來選分群、又用來測顯著，循環論證，一定顯著。要用獨立外部標籤定義群才合法。

**可能的追問 · 那個 82% / CN confound？**
> 原文口徑是 CN-獨立真表觀只有 8.1–17.6%，限 37 個印記 DMR、array 資料。我引它是當『CN confound 已有經典證據』的背書，所以我們在非印記位點也先固定 CN 再談甲基，**不外推成所有 ASM 都是 CN 假象**。

---

> **驗證足跡**：deck 全部數字經 wc1v407en（5 主線查驗，M3 9/9 + M5 11/11 機械重算 MATCH）+ w9456xpig（CP-D 雙審，overclaim 三紅線全 consistent）。reclassify 53.9% 為 report_claim 待落檔；其餘量化單樣本 HCC1395 PARTIAL（除 M3 Jaccard 全 7 樣本）。
