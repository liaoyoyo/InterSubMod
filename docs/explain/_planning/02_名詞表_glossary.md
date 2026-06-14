---
title: 統一名詞表（Unified Glossary）— ISM Subclonal Reconstruction 論文解釋頁複用
doc_type: glossary
status: planning
scope: 名詞表為單一真相來源（SoT），供所有解釋頁 / HTML / slide 複用；避免各頁定義漂移
created: 2026-06-12
source: 10 條研究線盤點 JSON 各自 glossary 去重合併（ism-core / ism-evolution / ism-vs-external / main-axis / methyl-phasing-assist / asm-zar1l-brca2 / asm-cis-confound / negatives / loh-phasing / background-glossary）
note: 名詞表只定義「概念」。專案特定數字（如某位點的 Δβ、某樣本 AUC、各線統計結果）一律不放此表，數字請引各研究線報告。
---

# 統一名詞表（Unified Glossary）

> **用法**：本表是所有解釋頁 / HTML / slide 的名詞**單一真相來源**。任何頁面引用某名詞時，定義以本表為準，不在各頁重寫定義（避免漂移）。
> **分組**：① haplotype 與 phasing ② 變異與 CN ③ 甲基化 ④ 統計與評估 ⑤ 工具與 pipeline。
> **標記**：🧩 = 高認知負荷、建議配 SVG 示意圖（讀者直覺易誤，需圖輔助）。
> **每條欄位**：嚴謹定義（資訊工程師科學敘述）／ 直覺類比 ／ 具體例子 ／ 常見誤解（若有）／ 關聯名詞。
> **數字政策**：本表不含專案特定數值（位點 Δβ、樣本 AUC、各線統計結論等），這些屬各研究線報告；本表只定義概念與計算口徑。

---

## ① Haplotype 與 Phasing

### Haplotype（單倍體型）
- **嚴謹定義**：同一個體在同一條染色體上，從單一親代遺傳到的一組等位基因組合。二倍體個體的每對同源染色體各帶一個 haplotype（記為 HP1 / HP2），可用 germline heterozygous SNP 作為錨點（anchor）推斷其歸屬（phasing）。
- **直覺類比**：你左腳和右腳穿的鞋來自兩個親代，每隻鞋上有不同記號（SNP）；看記號就知道這塊 DNA 屬於哪一條親代染色體。
- **具體例子**：一個 read 同時覆蓋三個 germline het SNP，這三個 SNP 的等位組合一致指向 HP1，則該 read 被歸為 HP1。
- **常見誤解**：把 HP1/HP2 當成「第幾號染色體」；實際是「同一對同源染色體的父本 vs 母本版本」。
- **關聯名詞**：HP-tag、Phasing、HP1-1（somatic sub-haplotype）

### HP-tag（HP:i: tag / HP:Z: tag）🧩
- **嚴謹定義**：BAM 格式中標記每條 read 所屬 haplotype 的整數 tag。在 longphase-S（paired）下：HP=1 為 germline HP1、HP=2 為 germline HP2；在 somatic haplotagging（longphase-S/-TO 6-state）下另有 HP=11（HP1-1）、HP=21（HP2-1）、HP3、unTag/HP0 等狀態。
- **直覺類比**：每條 DNA 分子讀完後在旁邊貼標籤：「來自父本（HP1）」「來自母本（HP2）」「來自帶突變的腫瘤亞株（HP1-1）」「身份不明（unphase/HP0）」。
- **具體例子**：一條帶 somatic ALT 的 read 落在 HP1 染色體鏈上，被標為 HP=11（即 HP1-1）。
- **常見誤解**：把 HP=11 當成「第三條染色體」或「copy number 標記」；它是 haplotype × somatic-event 的交集（見 HP1-1）。
- **關聯名詞**：Haplotype、HP1-1、longphase-S/-TO、6-state tag

### HP1-1（somatic sub-haplotype）🧩
- **嚴謹定義**：longphase-S/-TO 分配的 HP tag 值 11，代表「位於 HP1 germline 染色體鏈上、且攜帶 somatic ALT allele 的 reads 子群」。它是 HP1 鏈上出現體細胞突變後的子集，並非第三條染色體，也不是 copy number 標記。對稱項為 HP2-1。
- **直覺類比**：把父親那條染色體（HP1）上後天突變的那些 reads 單獨標出來；HP1-1 = 父親那條 + 突變，HP1 = 父親那條 + 沒突變。
- **具體例子**：比較 HP1（germline-HP1，無 somatic 突變）vs HP1-1（germline-HP1 + somatic ALT）的甲基化差異，此時染色體身份完全 hold constant，copy number 不是混淆因素。
- **常見誤解**：(1) 誤為第三單倍型；(2) 誤為 copy loss / copy number tag；(3) 誤為一定能被乾淨偵測（依賴 longphase-S somatic tag，LOH 區反而最難 phase）。
- **關聯名詞**：HP-tag、somatic haplotagging、HP-axis ASM、Phasing

### somatic haplotagging
- **嚴謹定義**：用 longphase-S 在 ONT reads 上同時標記體細胞突變來源（somatic / germline）與所屬單倍型（HP1/HP2/HP1-1/HP2-1），產生多態 haplotype tag 的流程。
- **直覺類比**：給每條 DNA read 蓋上雙重印章：「它來自哪條染色體副本」+「帶不帶腫瘤突變」。
- **具體例子**：longphase-S 對一個 paired tumor BAM 跑完後，每條 read 帶 HP tag，下游 ISM 才能按 HP 分群。
- **常見誤解**：以為 haplotagging 只分 HP1/HP2；somatic haplotagging 還能拆出 somatic sub-haplotype（HP1-1）。
- **關聯名詞**：HP1-1、Phasing、longphase-S/-TO

### Phasing（定相）
- **嚴謹定義**：依 germline heterozygous SNP 的共現模式，判斷同一 read 或染色體片段上的等位基因屬於 HP1 或 HP2 的過程。長讀因單條 read 夠長、能同時覆蓋多個 SNP，phasing 精度遠優於短讀。
- **直覺類比**：考古學家根據「哪些陶片出現在同一層土」判斷它們屬於同一個罐子；phasing 是根據「哪些 SNP 出現在同一條 DNA 分子」判斷它們屬於同一親代版本。
- **具體例子**：一條 read 跨越 5 個 het SNP，其等位組合與已知 phase-block 一致，即被定相到該 haplotype。
- **常見誤解**：以為 phasing 在所有區域都可行；在 LOH 區 germline het SNP 消失，phasing 失去錨點。
- **關聯名詞**：PS phase-block、LOH、LOH-constrained phasing、self-phasing circular dependency

### PS phase-block（phase-set block）
- **嚴謹定義**：germline VCF 中同一 phase-set（PS）的 SNP 群組，由 longphase / WhatsHap 在正常細胞中定相，組內 SNP 的相位互相一致。可作 held-out ground truth：遮住 read 自身的 germline 證據後，用其他訊號預測 HP，再對照 PS block 真值評估。
- **直覺類比**：考試時把某些題目答案遮起來，只用對其他題目學到的知識來猜，再對答案——比看著答案重做同樣題目更誠實。
- **具體例子**：在 methyl-assisted phasing 評估中，PS block 提供 read 真實 HP 歸屬，用以驗證「只用甲基化能否救回 unphase read」。
- **常見誤解**：以為 PS block 是腫瘤特異結構；它是 germline phasing 的產物。
- **關聯名詞**：Phasing、unphase read、held-out 驗證

### unphase read
- **嚴謹定義**：longphase-S 在 somatic haplotagging 後未能分配 HP tag 的 read，通常因局部缺乏足夠 germline het SNP 覆蓋而無法判斷歸屬。
- **直覺類比**：地圖上一條沒有地址的路，既沒被標「北邊」也沒被標「南邊」，因為附近沒有足夠路標。
- **具體例子**：落在 LOH 或 imprinting 區（germline SNP 稀疏）的 read 常被標為 unphase，這正是 methyl-assisted phasing 想救援的對象。
- **常見誤解**：以為「救援 unphase」在所有 unphase read 上都可行；多數 unphase read 因無本地 germline 錨點而無法建甲基參考——機制恰在最需要的 LOH/imprinting 區失效（chicken-egg）。
- **關聯名詞**：HP3、Phasing、PS phase-block、methyl-assisted phasing

### HP3（HP3 tag）
- **嚴謹定義**：longphase-S 6-state tag 之一，代表 read 帶有 somatic 突變但其所在 germline haplotype 未知（無足夠 germline het 來定相）。
- **直覺類比**：知道這孩子帶腫瘤身份（有 somatic 突變），但查不到它是父系還是母系族譜（germline haplotype 未知）。
- **具體例子**：在 H3→H1-1/H2-1 歸屬研究（T2）中，多數 H3 read 屬 no_germline，本身無 ground truth，故歸屬宣稱需謹慎。
- **常見誤解**：以為 H3 能直接外推到 H1-1/H2-1；H3 定義上無 germline 真值，直接外推是邏輯跳躍。
- **關聯名詞**：unphase read、HP1-1、somatic haplotagging

### 6-state tag（longphase-S/-TO 六態 tag）
- **嚴謹定義**：longphase somatic haplotagging 的 read 標記狀態集合：unTag(HP0) / H1 / H2（germline 雙倍型）/ H1-1 / H2-1（somatic sub-haplotype）/ HP3（somatic 突變但 germline 未知）。
- **直覺類比**：六個抽屜，每條 read 依「來自哪條染色體 + 帶不帶突變 + 身份明不明」放進對應抽屜。
- **具體例子**：methyl-phasing-assist 線在做任何 T1/T2/T3 分析前，先讀 longphase-S 原始碼確認此 6-state cascade。
- **常見誤解**：以為只有 HP1/HP2 兩態。
- **關聯名詞**：HP-tag、HP1-1、HP3、longphase-S/-TO

### LOH-constrained phasing
- **嚴謹定義**：在雜合性缺失（LOH）區域，腫瘤物理上只保留單一 haplotype，somatic SNV 因此形成「同-haplotype 子家族（Inner，HP1+HP1-1）」vs「跨-haplotype 比較（Outer，HP1+HP2-1）」的不對稱結構；以此方向性區分 somatic vs germline-like 變異。
- **直覺類比**：LOH 把其中一個染色體副本刪掉，突變只能落在剩下那條，讀出來的子族自然偏向同一邊。
- **具體例子**：LOH Inner 區的 NG=2 reads 以 same-haplotype 組成為主；非 LOH（Outer）區的 NG=2 以 cross-het 組成為主。
- **常見誤解**：把此 phasing signature 當「甲基化雙峰」訊號（HPFineNGroups 曾被誤解為甲基化，實為 phasing × allele 的 occupancy count）。
- **關聯名詞**：LOH、HPFineNGroups、Same-haplotype NG=2、by-construction circularity

### Same-haplotype NG=2（同單倍型雙群）
- **嚴謹定義**：一個 region 中 HPFineNGroups=2，且兩個被填滿的 bucket 屬同一單倍型（HP1+HP1-1 或 HP2+HP2-1）。在 LOH Inner 區此組成佔絕大多數，因 LOH 物理上只保留單一 haplotype，somatic SNV 必然產生 REF/ALT 子族均在同一條 haplotype 上。
- **直覺類比**：同一家人的兩個房間（父房 HP1 + 子房 HP1-1），而不是分屬兩家（父房 HP1 + 母房 HP2）。
- **具體例子**：對照 cross-het NG=2（HP1+HP2-1）= germline het 偽裝；same-hap NG=2 = LOH-driven somatic 訊號。
- **常見誤解**：把 same-hap rate 高直接當「真生物發現」，忽略 by-construction circularity（見下）。
- **關聯名詞**：HPFineNGroups、LOH-constrained phasing、by-construction circularity

### HPFineNGroups（HP 細分群數）🧩
- **嚴謹定義**：ISM 特徵（由 LabelTest.cpp 的 hp_to_fine_labels 計算），統計一個 region 中 {HP1, HP1-1, HP2, HP2-1} 四個 sub-haplotype bucket 有多少個被 reads 填滿（occupancy count）。純 phasing × allele presence，與甲基化矩陣無直接計算關係。NG=1 僅單一 sub-family 出現，NG=4 四個 bucket 皆有 reads。
- **直覺類比**：把 reads 依 HP tag 分成四個抽屜（HP1/HP1-1/HP2/HP2-1），HPFineNGroups 就是「有幾個抽屜有 reads」的計數——不是甲基化有幾個峰。
- **具體例子**：LOH Inner 的 NG=2 通常只填 HP1+HP1-1 兩格；Outer 的 NG=2 填 HP1+HP2-1 兩格。
- **常見誤解**：因名字像「HP 細分 N 個甲基化 group」而誤為甲基化雙峰程度；C++ 原始碼確認為 occupancy count（已更正）。
- **關聯名詞**：Same-haplotype NG=2、LOH-constrained phasing、Feature name 陷阱

### by-construction circularity（R-SELFREF，自參照偏差）🧩
- **嚴謹定義**：在 LOH-phasing 分析中，Inner same-HP1 bucket（HP1+HP1-1）的 HP1-1 成員本身由 somatic-attributed reads 定義，使「TP 比 FP 更多同單倍型分裂」的觀察不完全獨立於 pipeline 的標記決策，形成定義層面循環。決定性負控需用 --germline-hp-only flag 關閉 somatic HP tagging 後重觀察 bucket 是否塌陷。
- **直覺類比**：用一個人自己的護照照片來判斷他的身份，再聲稱「他長這樣」——在沒有獨立身份來源時是循環。
- **具體例子**：要打破循環，需跑 R-SELFREF（flag-on 全基因組 C++ 重跑）作 null 對照；此負控未完成前，phasing 脊柱證據封頂 Grade B+ 而非 Grade A。
- **常見誤解**：把循環依賴下的高 same-hap rate 當純生物學發現。
- **關聯名詞**：Same-haplotype NG=2、LOH-constrained phasing、Grade B+、self-phasing circular dependency

### self-phasing circular dependency（自相位循環依賴）🧩
- **嚴謹定義**：longphase 在 tumor-only（TO）模式下的結構性缺陷——somatic variant 的 ALT reads 同時被當作 phasing anchor（決定 HP 標記）與被評估對象（判斷 LOH），導致 ALT reads 系統性偏向同一 HP group，製造統計上的虛假 LOH 訊號。
- **直覺類比**：讓嫌疑人自己當法官決定自己有沒有罪，判決當然有利於自己；突變既是「定位標誌」又是「被定位對象」。
- **具體例子**：移除 self-phasing 後，大量 TO-mode TP LOH 消失；這是 TO 模式 QS 失效、LOH over-calling、甲基方向反轉等下游問題的共同根源。
- **常見誤解**：把 TO 模式的 LOH 訊號當真實生物結構。
- **關聯名詞**：by-construction circularity、Tumor-only (TO) pipeline、longphase-S/-TO、LOH

---

## ② 變異與 CN（Copy Number）

### Somatic mutation vs Germline mutation
- **嚴謹定義**：Germline mutation 是受精卵/生殖細胞即存在的序列變異，所有體細胞皆帶、會遺傳；Somatic mutation 是體細胞分裂中後天獲得、不遺傳，在腫瘤中大量累積，是 ISM 分析的目標突變。paired tumor-normal caller（ClairS）可區分兩者。
- **直覺類比**：germline 是出廠設定（所有細胞都有），somatic 是使用中某個細胞獲得的小傷痕（只傳給它的子代細胞）。
- **具體例子**：ClairS-TO 在無 normal 對照時 call somatic SNV，germline-somatic 分離不完全；ClairS paired full 在有 normal BAM 下精確分離。
- **常見誤解**：以為 caller 報出的都是 somatic；tumor-only 模式仍可能混入 germline。
- **關聯名詞**：caller_af、TP / FP、ClairS-TO、germline-het negative control

### TP / FP（True Positive / False Positive，真陽性 / 假陽性）
- **嚴謹定義**：相對 truth set（如 SEQC2、NYGC），caller 報出的 somatic SNV 中與真值相符者為 TP、不符者為 FP；truth 中漏掉者為 FN（False Negative）。ISM 多項研究的核心問題是「甲基化/ASM 能否判別 TP vs FP」。
- **直覺類比**：truth set 是標準答案，caller 是學生作答；對的是 TP、錯的是 FP、漏的是 FN。
- **具體例子**：用某特徵的 AUC 評估它能否分開 TP 與 FP；ISM 多輪結論為「甲基化判別力近隨機」。
- **常見誤解**：把「TP 上某現象更常見」直接當「可用過濾器」；需同時看靈敏度與跨樣本穩定性（見 anti-discriminative）。
- **關聯名詞**：caller_af、AUC、anti-discriminative、germline-het negative control

### caller_af（variant allele frequency，等位基因頻率）
- **嚴謹定義**：caller（如 ClairS-TO）報告每個候選位點的 ALT 等位基因頻率（ALT reads / total reads）。低 AF 位點常為次克隆突變或假陽性，是 TP/FP 最強的單一判別特徵之一。
- **直覺類比**：班級投票中「有多少比例同學舉手」；只有 5% 舉手（低 AF）很可能是噪音，80% 舉手（高 AF）才較可信。
- **具體例子**：在甲基→filter 研究中，甲基化訊號被證實主要是 caller_af 的代理（proxy），非獨立訊號。
- **常見誤解**：把甲基化與 caller_af 視為完全獨立的特徵；實際在低 AF 區兩者高度共線。
- **關聯名詞**：TP / FP、proxy signal（代理訊號）、vestigial covariate

### LOH（Loss of Heterozygosity，雜合性缺失）
- **嚴謹定義**：腫瘤在某染色體區段喪失原本的等位多樣性，僅保留單一單倍型；可因物理缺失（copy loss，CN=1）或 cnLOH（CN 不變但兩份 copy 同源）。LOH 區 germline het SNP 消失，使 haplotype phasing 失去錨點。ISM 以 LOH.bed 標記，用於區分 Inner（LOH 內）vs Outer（非 LOH）。
- **直覺類比**：兩本書（HP1/HP2）各有不同版本的某頁，LOH 像其中一本被撕掉只剩一版，或兩本都換成同一版。
- **具體例子**：LOH 區只能走 ALLELE-axis 比較（HP-axis 在一條 haplotype 丟失時不適用），故該區 ALLELE-axis ASM 不可宣稱 somatic 特異。
- **常見誤解**：把 LOH 一律當成 copy loss；cnLOH 的 copy number 仍是 2。
- **關聯名詞**：cnLOH、LOH-constrained phasing、ALLELE-axis ASM、Inner vs Outer

### cnLOH（Copy-Neutral LOH，拷貝數中性 LOH）
- **嚴謹定義**：LOH 子類型——基因座物理 copy number 維持 2（diploid），但兩份 copy 皆為同一等位（uniparental disomy 或 mitotic recombination 所致）。在純 copy number 分析中不顯示為缺失，需 allele-specific 方法才能偵測。
- **直覺類比**：本來兩本書各一本，cnLOH 是把 HP2 那本撕掉、再印一本 HP1 複本填補——架上仍兩本，但內容一樣了。
- **具體例子**：copy-neutral cnLOH 的 ASM 不是 copy dosage artifact（其 |Δβ| 與 dosage 方向可相反），是區分 ASM 真實性 vs copy 假象的關鍵案例。
- **常見誤解**：以為 LOH 一定伴隨 copy 數下降；cnLOH 反例打破此假設。
- **關聯名詞**：LOH、Copy Number (CN) confound、cis vs trans

### Copy Number (CN) confound / dosage confound 🧩
- **嚴謹定義**：當某 haplotype 的 copy 數增多時，其 reads 比例上升，會在表觀上造成「等位甲基差異」，但實際只是 copy 數量比例問題，而非真正的 allele-specific 調控。文獻報告 apparent ASM 有相當比例可被 copy 數解釋，故任何 ASM 宣稱須去 CN confound。
- **直覺類比**：投票看起來「A 派壓倒性多」，但其實只是 A 派人數本來就多（copy 多），不代表 A 派立場更強（真調控）。
- **具體例子**：HP-axis（HP1 vs HP1-1）設計上 held-constant CN/ploidy/alignment（兩組 reads 來自同一條 germline 染色體），是 de-confounded estimator；ALLELE-axis 則無此保護。
- **常見誤解**：把高 |Δβ| 直接當生物訊號，忽略它可能是 copy/coverage artifact。
- **關聯名詞**：HP-axis ASM、ALLELE-axis ASM、normal-anchored cis-test、copy-partition

### Inner vs Outer（LOH 內 vs 非 LOH 區）
- **嚴謹定義**：以 LOH.bed 劃分的區域類別。Inner = 落在 LOH 區段內的位點（只剩單一 haplotype，somatic ALT 只能同單倍型分裂）；Outer = 非 LOH 區（兩條染色體都在，germline het 與 somatic SNV 共存，可產生 cross-het 分裂）。
- **直覺類比**：Inner 像只剩一條路的村莊（車流只能走那條）；Outer 像雙線道（車流可分兩邊）。
- **具體例子**：Inner same-HP1 rate > Outer cross-het rate 的方向性，是 LOH-constrained phasing 的核心 signature。
- **常見誤解**：以為 Inner/Outer 只是位置標籤；它承載「somatic 分裂模式必然不同」的物理意義。
- **關聯名詞**：LOH、LOH-constrained phasing、Same-haplotype NG=2

### Feature name 陷阱（feature 名稱直覺解讀風險）
- **嚴謹定義**：特徵命名的字面語意可能與其 C++ 操作定義不符（如 HPFineNGroups 字面像「甲基化群數」實為 phasing occupancy count）。方法紀律要求分析新特徵前先讀 src/include 的定義，不依字面推論生物語意。
- **直覺類比**：菜單上「夫妻肺片」沒有夫妻也沒有肺；看名字點菜會誤解。
- **具體例子**：HPFineNGroups 被字面誤解為甲基化雙峰，C++ 原始碼確認為 {HP1,HP1-1,HP2,HP2-1} bucket 計數後才更正論文主軸。
- **常見誤解**：直接用特徵名推論其生物意義。
- **關聯名詞**：HPFineNGroups、方法紀律

---

## ③ 甲基化（Methylation）

### CpG（Cytosine-phosphate-Guanine）
- **嚴謹定義**：基因組中胞嘧啶（C）與鳥嘌呤（G）相鄰的二核苷酸（5'-CpG-3'）。哺乳動物多數 CpG 的 C 可被甲基化；CpG island（CGI）是 CpG 密度高的區段，常與啟動子相關。ISM 分析的甲基化位點均為 CpG（從 BAM MM/ML tag 定位）。
- **直覺類比**：CpG 是 DNA 上可被螢光筆畫線的特定字元；甲基化就是對這個字元做標記。
- **具體例子**：ISM 的 read×CpG 矩陣每一欄是一個 CpG 位點、每一列是一條 read。
- **常見誤解**：以為甲基化發生在任意 C；主要在 CpG context。
- **關聯名詞**：5mC vs 5hmC、MM/ML tag、β（甲基率）

### 5mC vs 5hmC（5-甲基胞嘧啶 vs 5-羥甲基胞嘧啶）
- **嚴謹定義**：5mC = 胞嘧啶 C5 位甲基化（最主要的 DNA 甲基化形式，與基因沉默相關）；5hmC = 5mC 的氧化形式（TET 酶催化，可能為主動去甲基化中間態）。ONT 可在 MM/ML tag 分別偵測（C+m? 為 5mC、C+h? 為 5hmC），是目前唯一能同時區分兩者的主流測序平台。
- **直覺類比**：5mC 是「標記」，5hmC 是「標記被改成另一符號」；都是 CpG 的化學修飾但生物意義不同。
- **具體例子**：ISM 的 C++ MethylationParser 目前只解析 C+m?（5mC-only），跳過 C+h?；Python 後處理層曾用 max-collapse 合成 any-modification。
- **常見誤解**：把 5mC 與 5hmC 當同一訊號；C++ 層 5mC-only 與 Python 層 max-collapse 是兩種不同的 5hmC 處理口徑，混用會出錯。
- **關聯名詞**：MM/ML tag、max-collapse、β（甲基率）、Δβ

### MM/ML tag（modification tag）
- **嚴謹定義**：BAM 中標記鹼基修飾的 SAM tag。MM 記錄修飾在序列上的位置（C+m? / C+h? 等，delta-encoded）；ML 記錄各修飾的機率（值 /255 = 機率）。ISM 用 MethylationParser 解析 MM/ML，配合 CIGAR 把 CpG 定位回參考座標。
- **直覺類比**：MM 是「在第幾個 C 上做了記號」的清單，ML 是「這個記號有多確定」的信心分數。
- **具體例子**：modkit 等工具也讀 MM/ML 提取 per-site 甲基率；ISM 則保留 per-read pattern。
- **常見誤解**：以為 ML 是 0/1；它是 0–255 的連續機率，硬二值化會丟棄信心資訊。
- **關聯名詞**：5mC vs 5hmC、CpG、β（甲基率）、soft 機率距離

### β（beta，甲基率）
- **嚴謹定義**：某 CpG 位點上的甲基化比例 = 甲基化 reads 數 / 該 CpG 被覆蓋的總 reads 數，範圍 0–1。是「位點層」甲基化的標準量。
- **直覺類比**：某題「答對的人數比例」。
- **具體例子**：β = 0.8 表示該 CpG 上 80% 的 reads 呈甲基化。
- **常見誤解**：把 β（位點層比例）與 read-level pattern 混為一談；β 已對 reads 取平均，丟失 between-read 結構。
- **關聯名詞**：Δβ、CpG、per-position 聚合、max-collapse

### Δβ（delta-beta，甲基率差）🧩
- **嚴謹定義**：兩個比較組在每個 CpG 上甲基率 β 的配對差值取平均，有方向、單位是甲基率（0–1，差值範圍 [−1,+1]）。負值代表後組偏低甲基（hypo）。⚠ 與 Δ（距離差）是完全不同的量。
- **直覺類比**：傳統意義的「甲基化差異」——同一 CpG，A 組 reads 有多少比例甲基化、B 組多少，相減。
- **具體例子**：HP-axis Δβ 比 HP1 vs HP1-1；ALLELE-axis Δβ 比 ALT vs REF。負 Δβ = 後組更 hypo。
- **常見誤解**：把 Δβ（有方向、率空間）與 Δ（無方向、NHD 距離空間）混淆；兩者算的是不同空間。
- **關聯名詞**：Δ（距離差）、β（甲基率）、HP-axis ASM、ALLELE-axis ASM

### Δ（距離差，LabelTest）🧩
- **嚴謹定義**：between_mean（異群 read 對的平均 NHD）減 within_mean（同群 read 對的平均 NHD），無方向、單位是 NHD（0–1）。只在 Δ>0（異群比同群更分散）才有生物意義並跑 permutation。
- **直覺類比**：問「不同群的 reads 甲基模式是否更分散」；Δ>0 才代表分群有意義。
- **具體例子**：HP-axis Δ 大表示「HP1-1 群 reads 彼此相似、與 HP1 群距離大」（結構分離），與同位點 Δβ（率差方向）是兩回事。
- **常見誤解**：把 Δ 與 Δβ 都當「差值」混為一談；Δ 講結構有無，Δβ 講率差方向。
- **關聯名詞**：Δβ（甲基率差）、NHD、PERMANOVA、軸 C

### ASM（Allele-Specific Methylation，等位特異性甲基化）
- **嚴謹定義**：同一基因組位置的兩條等位/單倍型上，甲基化程度有統計顯著差異的現象。本專案的 somatic ASM 特指因 somatic 突變誘發的 HP1 vs HP1-1 甲基差異，需與 germline het ASM（imprinting / random monoallelic）明確區分。
- **直覺類比**：DNA 像兩條平行鐵軌（父、母染色體），ASM 是其中一條鐵軌的「防鏽漆（甲基）」明顯比另一條多或少。
- **具體例子**：研究結論——ASM 真實存在、位點私有（somatic 因果）、跨癌種現象復現，但不是可用的 TP/FP 判別器（「存在」≠「可判別」）。
- **常見誤解**：把「ASM 存在」直接推論「ASM 可當 variant filter」；兩者是分開的 claim。
- **關聯名詞**：HP-axis ASM、ALLELE-axis ASM、strong-ASM、credible ASM、anti-discriminative

### HP-axis ASM（somatic-controlled 軸）🧩
- **嚴謹定義**：比較同一 somatic 位點的 HP1（germline-HP1 reads）vs HP1-1（同一 germline 染色體上攜 somatic ALT 的 sub-haplotype reads）的甲基差異。兩組來自同一條 germline 染色體，germline allele 效應與 CN/ploidy/alignment 被 held-constant，是 somatic-controlled estimator。
- **直覺類比**：只比「同一條鐵軌（HP1）上，突變區段（HP1-1）vs 非突變區段（HP1）的防鏽漆」，排除兩條鐵軌天生不同的影響。
- **具體例子**：HP-axis 是去 confound 的正確軸；任何 somatic ASM 宣稱應走此軸並附 germline-het negative control。
- **常見誤解**：以為 LOH 位點也能走 HP-axis；LOH 區一條 haplotype 丟失，HP-axis 不適用，只能走（受 confound 的）ALLELE-axis。
- **關聯名詞**：HP1-1、ALLELE-axis ASM、Copy Number confound、normal-anchored cis-test

### ALLELE-axis ASM（confounded 軸）🧩
- **嚴謹定義**：比較 somatic 位點的 ALT-allele reads vs REF-allele reads 的甲基差異。在 germline het 位點，ALT/REF reads 分屬不同 germline haplotype（HP1 vs HP2），此差異混入「兩條 germline 染色體本就存在的 baseline allelic methylation（imprinting / random monoallelic）」，無法單純歸因 somatic 事件。
- **直覺類比**：比較「兩條染色體」，但兩條本來就有甲基差異（如 imprinting），所以知道「有差」卻分不清是突變造成還是天生。
- **具體例子**：LOH 位點只能走 ALLELE-axis，故其 ASM 不可宣稱 somatic 特異；需與 germline-het null 對照（若 TP 率不高於 null 即被否定）。
- **常見誤解**：把 ALLELE-axis 顯著當 somatic-specific ASM。
- **關聯名詞**：HP-axis ASM、germline-het negative control、Copy Number confound、ASM

### strong-ASM（強 ASM 位點）
- **嚴謹定義**：同時滿足統計顯著（如 Bonferroni 校正後顯著）與效果量門檻（如 |Δβ|≥0.1）的 ASM 位點。全基因組占比極低，方向（hypo/hyper）通常無偏好。
- **直覺類比**：在所有甲基差異中，挑「差異夠大且非偶然」的那些位點。
- **具體例子**：strong-ASM 在 FP 中反而富集（anti-discriminative）；其膨脹根因是 FP 多落在低覆蓋/近飽和 baseline/高 LOH 的 regime。
- **常見誤解**：把 strong-ASM 當「強 TP 訊號」；它反而 anti-discriminative。
- **關聯名詞**：credible ASM、anti-discriminative、regression-to-extreme、ASM

### credible ASM（可信 ASM）
- **嚴謹定義**：通過 regime filter（如 n_cpg 足夠 + germline baseline 非極端 0.5±0.3 + 優先 nonLOH）與 blind-ARI gate 後的 ASM 位點，排除低覆蓋 / 極端 baseline / 高 LOH 的 artifact regime。
- **直覺類比**：在 strong-ASM 中再篩掉「看起來差異大但其實只是 reads 太少、baseline 近飽和、或 LOH 假象」的位點。
- **具體例子**：credible ASM 多落在 promoter±2kb 附近，是 characterization 用、非 filter 用。
- **常見誤解**：把所有 strong-ASM 都當可信；需 regime gate 過濾。
- **關聯名詞**：strong-ASM、LOH-coverage-baseline 三角、blind-ARI、ASM

### regression-to-extreme（極端回歸，FP 甲基富集機制）
- **嚴謹定義**：高 |Δβ| 位點在 FP 中過度富集的機制解釋——這些位點的大甲基差是 read 分群結構缺失 / LOH imbalance / 低有效覆蓋的 artifact，而非 somatic 驅動的真實差異（測量不穩定性放大假象）。
- **直覺類比**：最「顯眼」的甲基差異恰好集中在雜訊最多的地方（FP），是不穩定放大的假象，不是生物訊號。
- **具體例子**：FP 的 strong-ASM 多伴隨無 clustering / 高 LOH / 小 subhaplotype，三維互相關聯。
- **常見誤解**：把大 |Δβ| 當強訊號。
- **關聯名詞**：strong-ASM、anti-discriminative、LOH-coverage-baseline 三角

### LOH-coverage-baseline 三角（ASM artifact regime）🧩
- **嚴謹定義**：FP 的 strong-ASM 多落在三個互相關聯的維度——高 LOH 比例、低有效覆蓋（n_cpg 少）、極端 germline baseline（近 0 或近 1）。機制鏈：LOH → 只剩單一單倍型 reads → 有效覆蓋降低 → baseline 可能飽和 → 微小噪音造成看似大的 |Δβ|。
- **直覺類比**：三個看似獨立的發現，其實是同一機制鏈的三個切面。
- **具體例子**：可信 ASM 篩選準則正是反向避開此三角（足夠覆蓋 + baseline 非極端 + 優先 nonLOH）。
- **常見誤解**：把三個維度當三個獨立 confound 分別處理。
- **關聯名詞**：credible ASM、regression-to-extreme、Copy Number confound、LOH

### normal-anchored cis-test 🧩
- **嚴謹定義**：以匹配正常樣本（normal）的 HP 甲基為 germline baseline，判斷腫瘤的 allelic 甲基差是否為 somatic-cis 驅動。三方/嚴格設計：A=normal-HP1、B=tumor-HP1、C=tumor-HP1-1，計算 d_somatic=C−B、d_cis=C−A、d_drift=B−A；以「d_cis 大且相對 d_drift 夠大」為 cis 判準。設計 held-constant copy/ploidy/alignment。
- **直覺類比**：對照組實驗——normal 是基線，tumor-HP1 是「同單倍型但正常組」，tumor-HP1-1 是「突變那群」；只有突變群與 normal 的差遠大於正常群與 normal 的漂移，才算真 cis。
- **具體例子**：是 ISM 宣稱的差異化貢獻之一（已發表工具中罕見同時具備）；但只在有 matched normal BAM 的樣本才完整。
- **常見誤解**：(1) 把 ISM 的 HP_Residual（廣義 HP-family 殘差）等同窄義三方 cis-test；(2) 以為高甲基差就是 cis（可能是 drift 或 copy）。
- **關聯名詞**：HP_Residual、cis vs trans、d_cis / d_drift / d_somatic、copy-partition、HP-axis ASM

### cis vs trans（驅動方向）
- **嚴謹定義**：cis-driven 指甲基差由同一染色體上局部序列/突變直接驅動，不受對側等位影響；trans-driven 指由遠端調控、基因表現量、或對側等位事件間接驅動。判斷 cis 困難，因 copy number 變化可造成表觀等位差異但實為 dosage。
- **直覺類比**：cis 像「同一條街施工影響這條街的交通」；trans 像「別條街施工透過改道影響這條街」。
- **具體例子**：normal-anchored cis-test 用以在腫瘤中分離 cis-driven ASM 與 germline allelic / copy 假象。
- **常見誤解**：把任何 allelic 甲基差都當 cis。
- **關聯名詞**：normal-anchored cis-test、Copy Number confound、cnLOH

### HP_Residual（ISM 欄位）
- **嚴謹定義**：ISM significance_summary 的 HP_Residual_Delta 欄位 = tumor 樣本 HP delta（HP1 vs HP2 甲基差）減 normal 樣本 HP delta，量的是 germline-HP-family 層級的 tumor-specific 殘差，作 cis 代理。注意：屬「廣義」殘差，非窄義三方 cis-test（tumor HP1-1 vs tumor HP1 vs normal HP1）。
- **直覺類比**：「tumor 裡兩條染色體的差」減掉「normal 裡兩條染色體的差」= 癌症特有的差異；但仍是兩條染色體相比，不是 somatic subclone vs 正常 clone 的精準三方比較。
- **具體例子**：用 HP_Residual 篩 cis 候選時，FP 的 cis 候選率可能最高——即 cis 殘差不判別 TP/FP。
- **常見誤解**：把 HP_Residual 當窄義 cis-test。
- **關聯名詞**：normal-anchored cis-test、cis vs trans、ALLELE-axis ASM

### copy-partition（copy 分解）
- **嚴謹定義**：把觀察到的 HP 甲基差拆解為「copy 數差驅動（d_copy）」與「within-copy cis 殘差（d_within）」的分析，用以判斷甲基差是否被 copy 主導。塌陷判準如 |d_within| < 0.5|d_HP| 表示 copy 主導、非真 cis。
- **直覺類比**：把「兩堆球數量不同造成的平均差」與「同一堆內球的真實差」拆開算。
- **具體例子**：某啟動子位點 d_copy 大但 d_within 邊際 → 判定 copy/subclone 主導、非乾淨 cis。
- **常見誤解**：把通過 cis-test 候選的位點都當乾淨 cis，未做 copy-partition。
- **關聯名詞**：normal-anchored cis-test、Copy Number confound、d_cis / d_drift / d_somatic

### max-collapse（5mC/5hmC 最大值合成）
- **嚴謹定義**：對同一 read/CpG 的 5mC 機率與 5hmC 機率取最大值，得 any-modification 代理值。是修正 MSA Level1「對雙修飾 BAM 每 (read,CpG) 發 5mC+5hmC 兩列分別計數導致 β 砍半」artifact 的 Python 層做法。
- **直覺類比**：有兩個甲基訊號（紅燈與黃燈），取最亮的代表「有燈亮」；好處是合成 any-mod，壞處是不知亮的是哪盞、且 5hmC 弱訊號被 5mC 蓋掉。
- **具體例子**：C++ 層為 5mC-only（跳過 5hmC），Python 層為 max-collapse；兩層 5hmC 處理口徑不同，混用會誤判。
- **常見誤解**：以為 C++ 也 max-collapse；C++ 是 5mC-only。
- **關聯名詞**：5mC vs 5hmC、β（甲基率）、MSA、Dirichlet-Multinomial 分軌

### epipolymorphism（甲基表觀多態性 / WSH disorder）
- **嚴謹定義**：在小窗口（如 4-CpG）內計算所有 epiallele 組合的機率分佈（如 Simpson 多樣性 1−Σpₖ²），度量同一群 reads 甲基 pattern 的隨機失序程度（stochastic disorder），而非系統性 allele-specific 差異。文獻（如 Landau 2014 PDR）明文用以排除 ASM 來源。
- **直覺類比**：問「同一條染色體上的 reads 甲基模式是否亂成一鍋粥」；高 = 像亂碼、低 = 像整齊兩類。
- **具體例子**：ISM 把 NME/epipolymorphism 當「對照/disorder 軸」而非主引擎——ISM 量的是「結構（structure）」不是「失序（disorder）」。
- **常見誤解**：把高 epipolymorphism 當 ASM 訊號；它量的是 disorder，目的恰是排除 ASM。
- **關聯名詞**：軸 E（disorder/WSH scalar）、ASM、PERMANOVA、結構 vs 失序

### germline-het negative control
- **嚴謹定義**：宣稱 somatic ASM 時，同時對 germline heterozygous SNP 位點（純 germline，非 somatic）跑同一 pipeline，比較其 strong-ASM 率（null）與 somatic 位點的 strong-ASM 率；若 somatic 率不高於 null，則 somatic ASM 宣稱被否定。是方法紀律的必要對照。
- **直覺類比**：問「隨機找一個普通的 germline 雜合 SNP，同樣分析會有多少假陽性？」somatic 必須比這更高才算真。
- **具體例子**：ALLELE-axis 的 somatic 率可能低於 germline-het null（confounded）；HP-axis 的 somatic 率高於 null（modest 但成立）。
- **常見誤解**：宣稱 somatic ASM 而不跑此對照。
- **關聯名詞**：ALLELE-axis ASM、HP-axis ASM、null（虛無分佈）、ASM

### per-position 聚合 vs per-read 表示
- **嚴謹定義**：per-position 聚合（傳統，modkit/pycoMeth/DSS 等軸 A）對每個 CpG 在所有 reads 上取平均甲基率，丟失 read 之間的共甲基化結構；per-read 表示（ISM 軸 C）建 read×CpG 矩陣，保留每條 read 的完整甲基 pattern，使 between-molecule 結構可被分析。
- **直覺類比**：per-position 是「看每題平均答對率」；per-read 是「看每份問卷」——平均 70 分可能是全 70，也可能一半 0 一半 100。
- **具體例子**：若一半 reads 全甲基、一半全未甲基，per-position 均值看不到差異，但 read×CpG 矩陣與距離法能看出兩個分離的雲。
- **常見誤解**：以為「HP1 vs HP2 取平均再相減」就夠；那丟失了共甲基化結構（亞群訊號）。
- **關聯名詞**：read×CpG 矩陣、軸 C、PERMANOVA、β（甲基率）

### methyl-assisted phasing（甲基救援 phasing）
- **嚴謹定義**：用每條 read 的甲基化 pattern 輔助判斷其 haplotype 歸屬，以救援無法被 germline SNP 定相的 read（unphase）或糾正打錯的 tag。機制天花板：甲基是 germline-haplotype 層級訊號——分「不同 haplotype 強、within-haplotype 弱」。
- **直覺類比**：甲基像「家族口音」（germline haplotype 層級）——能分北方 vs 南方口音（不同 haplotype），但分不出同是南方人的「老家族 vs 新移入」（同 haplotype 內亞群）。
- **具體例子**：T1（unphase→H1/H2）有利、T3（拆同一 haplotype 內亞群）不利；且真實 unphase 多在 LOH/imprinting 區（無本地錨點）無法救援——機制在最需要處失效（chicken-egg）。
- **常見誤解**：把「甲基攜帶 haplotype 資訊」當一刀切特性；它只在 between-haplotype 強、within-haplotype 弱。
- **關聯名詞**：germline-haplotype 層級甲基、unphase read、anchor AUC、PS phase-block

### germline-haplotype 層級甲基
- **嚴謹定義**：DNA 甲基化在同一個體兩條同源 haplotype 間的系統性差異，由 germline cis 調控（cis-mQTL、imprinting）決定，在正常二倍體組織與腫瘤中皆穩定存在。可用 matched normal（copy-clean 二倍體）驗證其 HP 甲基分離度，排除 copy artifact。
- **直覺類比**：父本與母本染色體在某些基因上天生就有不同的開關設定，細胞複製時被傳承，非腫瘤才突然出現。
- **具體例子**：matched normal 的 HP 甲基分離度可 ≥ tumor，且 LOH 區 tumor 分離度反而最低——方向與「copy artifact 假說」相反，證實是 germline 層級訊號。
- **常見誤解**：把 HP 甲基分離當 somatic-specific；它主要是 germline 層級。
- **關聯名詞**：methyl-assisted phasing、ASM、HP-axis ASM、Copy Number confound

---

## ④ 統計與評估

### NHD（Normalized Hamming Distance，標準化漢明距離）
- **嚴謹定義**：對兩條 read 的二元甲基 pattern（每 CpG 甲基=1 / 未甲基=0）算 Hamming 距離後除以共同覆蓋的 CpG 數，範圍 0–1（0=完全一樣、1=完全相反）。需雙方共同覆蓋 ≥ 一個最小 CpG 數（如 C_min=5），不足者記 MAX_DIST=1.0（懲罰式）。
- **直覺類比**：兩條讀段各自的甲基「指紋」有多少 CpG 位置不一樣，除以可比較的 CpG 數。
- **具體例子**：是 ISM read-read 距離矩陣的預設 metric，下游接 UPGMA 聚類與 PERMANOVA。
- **常見誤解**：把低重疊對的 MAX_DIST=1.0 fallback 當真實大距離；它可能造假雙峰結構。
- **關聯名詞**：read×CpG 矩陣、Δ（距離差）、PERMANOVA、軸 C、soft 機率距離

### read×CpG 矩陣 → 距離矩陣
- **嚴謹定義**：ISM 先建 reads（列）× CpG（欄）的二元甲基矩陣（二值化：高>0.8=甲基、低<0.2=未甲基、中間=ambiguous），再對每對 read 算距離（NHD 等）得 N×N 顯式距離矩陣，作為 UPGMA 聚類與 PERMANOVA 的輸入。
- **直覺類比**：把每條讀段想成一份答卷（每個 CpG 填 0/1），再算「任兩份答卷有多不像」。
- **具體例子**：這是 ISM 軸 C（between-read 距離）的核心資料結構，保留 between-molecule 結構。
- **常見誤解**：以為 ISM 直接比 per-CpG 平均；它保留 read-level pattern。
- **關聯名詞**：NHD、per-read 表示、UPGMA、PERMANOVA、軸 C

### UPGMA（階層聚類）
- **嚴謹定義**：以距離矩陣做的凝聚式階層聚類（Unweighted Pair Group Method with Arithmetic mean）；ISM 用 TreeCutter 配 silhouette score 在 k=2..10 自動選最佳群數。
- **直覺類比**：把最像的兩份答卷先併在一起，反覆合併成一棵樹，再從樹上剪出最自然的群數。
- **具體例子**：UPGMA 樹 + silhouette 給出 read 的甲基亞群分群，再疊 HP label 評估與 HP 對齊程度。
- **常見誤解**：把 UPGMA 群數當固定 2；ISM 自動選 k。
- **關聯名詞**：read×CpG 矩陣、距離矩陣、silhouette、Bernoulli-mixture EM

### PERMANOVA（Permutational MANOVA）🧩
- **嚴謹定義**：在 read-read 距離矩陣上計算 pseudo-F =（組間距離 SS/(k−1)）/（組內距離 SS/(N−k)），再用標籤置換（如 999 次）建 null 分佈求 p。不假設多變量常態，問的是「按 HP 標籤分組後，組間甲基距離是否顯著大於組內」。對稀疏列聯表比 Cramér's V 穩健。
- **直覺類比**：把 reads 的 HP 標籤隨機洗牌 N 次，看真實分群是否比 99% 的隨機洗牌都更整齊；整齊才說明 HP 真的在解釋甲基差異。
- **具體例子**：若一半 reads 全甲基、一半全未甲基，per-CpG 均值看不到差，但 PERMANOVA pseudo-F 大（有結構）；epipolymorphism 高的混雜情況 pseudo-F 小（無結構）。
- **常見誤解**：把它當「per-CpG 率差檢定」；它是多變量「結構有無」檢定（結構 ≠ 率差）。
- **關聯名詞**：NHD、Δ（距離差）、軸 C、Cramér's V reliability gate、結構 vs 失序

### Fisher exact test（per-CpG）+ over-dispersion 問題 🧩
- **嚴謹定義**：ISM per-CpG ASM 用 2×2 列聯表 Fisher exact test，假設每條 read 對同一 CpG 的 call 為獨立 Bernoulli（dispersion φ=0），within-region 做 BH-FDR。長讀中同 haplotype/clone 的 reads 共享來源、非獨立（φ>0），Fisher 低估 variance → p 偏小 → 顯著 CpG 偏多（ASM 存在率偏樂觀）。
- **直覺類比**：擲 100 枚硬幣若獨立會分散在平均附近；但若硬幣受同一批次影響（彼此相關），實際方差更大——同 clone 的長讀正是這種。
- **具體例子**：建議改 beta-binomial（顯式吸收 over-dispersion），或保留 Fisher 但明標 anti-conservative + 加 effect-size 門檻。
- **常見誤解**：以為 Fisher 在長讀也保守；其獨立性假設使它在 ASM 存在率上偏樂觀。
- **關聯名詞**：over-dispersion、beta-binomial、BH-FDR、多重檢定校正

### over-dispersion（超分散）
- **嚴謹定義**：二項模型中實際觀測 variance 超過理論 Binomial np(1−p)，常因觀測間正相關（ρ>0）引起。甲基化中同 haplotype/clone 的 reads 甲基狀態相關，使 read 間變異超過獨立 Bernoulli 預測（φ>0）。
- **直覺類比**：每次擲的是同一枚有偏硬幣（reads 來自同一 clone），正反面結果相關；Fisher 以為你擲了更多次獨立硬幣，低估了不確定性。
- **具體例子**：是 Fisher 假陽性偏多、ASM 存在率偏樂觀的根因之一。
- **關聯名詞**：Fisher exact test、beta-binomial、ASM

### beta-binomial（貝塔二項分布）+ dispersion shrinkage
- **嚴謹定義**：Binomial(n,π) 中令 π 服從 Beta(α,β)，邊緣即 beta-binomial，variance = np(1−p)(1+(n−1)φ)，φ 為 over-dispersion 參數（φ→0 退化 Binomial）。DSS 風格用 lognormal prior 跨 CpG 做 empirical-Bayes shrinkage，使低覆蓋 CpG 的 φ 向全基因組均值收縮，再用 Wald 檢定甲基差。
- **直覺類比**：把「硬幣偏差率」本身當隨機變數，允許同組 reads 各帶不同甲基傾向，比 Fisher 更貼近長讀生物學；shrinkage 像「你家 3 個樣本不夠，借隔壁一千家平均水電費來穩定估計」。
- **具體例子**：是業界黃金標準 DSS 的核心，建議用以取代或補強 ISM 的 Fisher。
- **常見誤解**：以為 beta-binomial 只是換個檢定；它顯式量化 over-dispersion 並穩定單樣本估計。
- **關聯名詞**：over-dispersion、Fisher exact test、DSS、ASM

### BH-FDR（Benjamini-Hochberg False Discovery Rate）
- **嚴謹定義**：多重檢定中對排序後 p 值做 q_i = p_i × m / rank_i 調整並單調修正，控制預期假陽性比例（FDR）在閾值（如 5%）以下。ISM 在 within-region per-CpG 層有 BH-FDR，但跨 region / genome-wide 層缺多重檢定校正。
- **直覺類比**：同時做 1000 個檢定若都用 p<0.05，平均 50 個假陽性；BH-FDR 依排名動態調小閾值，使假陽性「比例」（非個數）控制在 5%。
- **具體例子**：掃全基因組時若 global_p 取多 p 最小值又不校正，ASM 存在率假陽性會膨脹——建議跨 region BH-FDR + global_p 用 Sidak/Bonferroni。
- **常見誤解**：以為 within-region BH-FDR 就涵蓋全基因組；跨 region 仍需校正。
- **關聯名詞**：多重檢定校正、Fisher exact test、PERMANOVA

### 多重檢定校正（跨 region / genome-wide）
- **嚴謹定義**：在掃描數萬位點時，對 region 級 p 跨位點做 FDR，並對組合 p（如 global_p = min(p1,p2,p3)）用 Sidak（1−(1−min_p)^k）或 Bonferroni 校正，避免「取多個 p 最小值」造成的假陽性膨脹。
- **直覺類比**：抽很多次獎，總會抽到「中獎」；校正是把「抽很多次」這件事算進去，免得把運氣當實力。
- **具體例子**：被列為 ISM 偏樂觀來源中「除 Fisher 外最重要」的統計 gap。
- **常見誤解**：以為 min-of-k p 值不需校正。
- **關聯名詞**：BH-FDR、Fisher exact test、global_p

### Cramér's V reliability gate（Cochran gate）🧩
- **嚴謹定義**：對 HP-cluster 列聯表計算 Cramér's V（關聯強度），並用 Cochran 準則（所有格期望值≥5）判斷卡方近似是否有效；不滿足時 gated 值在 summary 輸出 reliable?v:0（歸零），但 significance.json 仍保留原始未閘控值（雙口徑）。稀疏表下 gated 與 raw 可天差地遠。
- **直覺類比**：Cramér's V 是用卡方推導的，需每格樣本夠多才近似成立；Cochran gate 像把關門說「你的表太稀疏，這個強度數字不算數」。
- **具體例子**：一個位點 gated CramérV=0 但 PERMANOVA 顯著，不是矛盾——是稀疏表使卡方不可靠（見 latent 真結構）。
- **常見誤解**：(1) 把 gated 與 raw 當同一欄混用；(2) 把 CramérV=0 當「無分群」（可能只是 reliability gate 觸發）。
- **關聯名詞**：latent 真結構、PERMANOVA、Cochran 準則、雙口徑

### latent 真結構（latent real structure）
- **嚴謹定義**：gated Cramér's V=0（被 Cochran 閘控歸零）但 PERMANOVA（距離法，對稀疏穩健）顯著的 ASM 位點——reads 確實按 HP 分群，只是卡方統計因稀疏列聯表判定不可靠而被過濾。屬「被保守統計過濾掉的真訊號」。
- **直覺類比**：reads 的分群（距離角度）確實存在且與 HP 對齊，但因樣本少讓卡方說「我沒把握」——真信號被保守方法濾掉。
- **具體例子**：這類位點適合以 PERMANOVA gate 做 characterization，但不可寫成 TP/FP filter。
- **常見誤解**：把 latent 位點當 noise 丟棄（只看 Cramér's V）。
- **關聯名詞**：Cramér's V reliability gate、PERMANOVA、credible ASM

### AUC（Area Under the ROC Curve）
- **嚴謹定義**：ROC 曲線下面積，衡量二元分類器在所有閾值下的區分力，範圍 0–1；0.5=隨機、1.0=完美。ISM 用以評估特徵能否分 TP vs FP；多輪研究以 AUC ≤ 約 0.58 確立 filter NEGATIVE。
- **直覺類比**：把一顆 TP 球與一顆 FP 球混抽，AUC 是「TP 球分數較高」的機率。
- **具體例子**：甲基化 |Δβ| 預測 is_TP 的 AUC 近 0.5（隨機）→ 不判別。
- **常見誤解**：AUC 高 ≠ 可用閾值（見「AUC 高但 FP removal=0%」）；且絕對 AUC 在某些設計會系統性膨脹（見 anchor AUC）。
- **關聯名詞**：TP / FP、anchor AUC、null（虛無分佈）、安全約束下 FP removal=0%

### anchor AUC（leave-one-out centroid AUC）🧩
- **嚴謹定義**：用有 HP1/HP2 tag 的 reads 建兩群甲基質心，遮住目標 read 標籤後用「距最近質心」預測 HP，與真實 HP 對照算 AUC。注意：在有充足 anchor read 的區域會系統性膨脹（germline het null median 遠高於 0.5），不應直接當 effect size，須引相對 null（delta vs shuffle）或 held-out rescue rate。
- **直覺類比**：用班上其他人的成績分出 A/B 組平均（質心）再預測新人——準確但因組內差異小而顯得「過於準確」。
- **具體例子**：methyl-assisted phasing 的絕對 anchor AUC 膨脹是方法樂觀；誠實量是 delta-over-null 或 held-out 救援率。
- **常見誤解**：把絕對 anchor AUC 當真實 effect size。
- **關聯名詞**：AUC、null（虛無分佈）、methyl-assisted phasing、germline-het negative control

### null（虛無分佈 / permutation null）
- **嚴謹定義**：透過隨機化（label shuffle、HP-shuffle permutation 等）建立「無真實效應」下的統計量分佈，作為觀測值的對照基準。許多絕對指標（如 anchor AUC、ASM rate）的 null 本身就遠離 0.5/0，故須用「觀測 − null」（excess/delta）而非絕對值比較。
- **直覺類比**：賭博裡「比莊家多贏的比例」——不是看總共贏多少，而是扣掉純靠運氣本來就能贏的部分。
- **具體例子**：excess-over-null（見下）正是把覆蓋/樣本大小造成的 null 漂移扣掉後的誠實效應量。
- **常見誤解**：假設 null 一定接近 0.5/0；許多設計下 null 很高。
- **關聯名詞**：excess-over-null、anchor AUC、permutation test、germline-het negative control

### excess-over-null（超出隨機基準）
- **嚴謹定義**：某位點集合的觀測指標（如 somatic-controlled HP-axis ASM rate）減去 permutation null 的同指標，量化「比隨機多出多少」，排除覆蓋與樣本大小對 raw rate 的影響。不同樣本須用 excess 而非 raw rate 比較。
- **直覺類比**：高覆蓋與低覆蓋樣本的 null 不同；只看 raw rate 會把低覆蓋樣本本來就高的 null 誤當訊號。
- **具體例子**：跨樣本 ASM 復現須看 excess>0 是否一致，而非 raw |Δβ| 大小。
- **常見誤解**：直接比各樣本 raw rate。
- **關聯名詞**：null（虛無分佈）、ASM、跨樣本復現

### LOSO（Leave-One-Sample-Out，留一樣本法）🧩
- **嚴謹定義**：跨樣本留一驗證——N 個樣本中用 N−1 個訓練、對留下 1 個預測，重複 N 次取平均，評估模型跨樣本的真實泛化。對照 within-sample k-fold（只在同一樣本內分割，無泛化挑戰）。
- **直覺類比**：用 A/B/C/D 四校的題訓練、去考 E 校的題；而不是用 E 校 80% 題訓練再考 E 校剩 20%（後者作弊成分高）。
- **具體例子**：甲基 filter 在 in-distribution k-fold 看似有提升，但 LOSO held-out 近零、跨樣本 transfer 甚至大負——揭露效益全來自 sample-level circularity。
- **常見誤解**：把 within-sample k-fold 當跨樣本驗證。
- **關聯名詞**：sample-level circularity、proxy signal、vestigial covariate、TP / FP

### sample-level circularity（樣本層循環偏差）
- **嚴謹定義**：當訓練集與驗證集來自同一生物樣本（即使做 k-fold 分割），模型學到的是該樣本特有分佈而非可泛化規律，效能虛高。LOSO 與 within-sample k-fold 的巨大落差即為其量度。
- **直覺類比**：用甲同學 100 題訓練、再用甲同學另 20 題考試——只是在背答案，跟去考乙同學完全不同難度。
- **具體例子**：甲基→FP filter 方向的 in-distribution 提升被證實 100% 來自此循環。
- **常見誤解**：以為交叉驗證自動排除循環；同樣本內分割不行。
- **關聯名詞**：LOSO、proxy signal、vestigial covariate

### proxy signal（代理訊號）
- **嚴謹定義**：某特徵的預測力其實來自它與另一主導特徵的共線性，而非自身獨立資訊。甲基化在 filter 框架中被證實主要是 caller_af 的 proxy（低 AF subclone TP 恰有特定甲基模式），模型實際在重學 AF 分佈。
- **直覺類比**：甲基化是依附在 caller_af 上的「搭便車乘客」，不是司機。
- **具體例子**：消融甲基化特徵後模型效能幾乎不變（甚至略升）→ 確認其為 proxy / vestigial。
- **常見誤解**：把共線特徵的表觀預測力當獨立貢獻。
- **關聯名詞**：caller_af、vestigial covariate、sample-level circularity

### vestigial covariate（殘餘共變數）
- **嚴謹定義**：多變數迴歸中與主導特徵高度共線、獨立預測力極低（消融後效能幾乎不變甚至略升）的特徵。在 L2 ridge 下其小係數被壓縮但殘餘噪訊仍可輕微干擾決策邊界，故移除反而可能略好。
- **直覺類比**：蛇的退化後肢——結構仍在，對行動毫無貢獻。
- **具體例子**：甲基化在 LR filter 中排第 5 順位，主導者為 caller_af → LOH_inner → Coverage → NG。
- **常見誤解**：以為「加特徵不會更差」；共線殘餘噪訊可使效能略降。
- **關聯名詞**：proxy signal、caller_af、LOSO、消融（ablation）

### 消融（ablation，消融實驗）
- **嚴謹定義**：系統性移除模型的一組特徵後比較效能變化，量化該組特徵的獨立貢獻。是判定 vestigial covariate / proxy signal 的標準方法。
- **直覺類比**：拆掉零件看機器還能不能跑，藉此判斷零件是否必要。
- **具體例子**：移除全部甲基化特徵後 ΔF1 只掉極小量（占原提升的個位數百分比）→ 甲基化 vestigial。
- **常見誤解**：以為消融只用於模型壓縮；它也是因果歸因工具。
- **關聯名詞**：vestigial covariate、proxy signal

### Cohen ribbon（效果量門檻帶）
- **嚴謹定義**：依 Cohen's d 或類似標準預先設定「具實用意義」的最小效益閾值帶（如 ΔF1 ≥ +0.005）。低於此門檻即使統計顯著仍視為 marginal。
- **直覺類比**：不問「有沒有差」，問「差到值不值得用」；血壓降 0.1 mmHg 統計顯著但醫師不會換藥。
- **具體例子**：甲基 filter pilot 的 ΔF1 落在 Cohen ribbon 之下，判為 marginal/PARTIAL。
- **常見誤解**：把統計顯著等同實用顯著。
- **關聯名詞**：effect size、AUC、Wilcoxon signed-rank

### Wilcoxon signed-rank（符號秩檢定）+ 方向一致性 ≠ effect 大小 🧩
- **嚴謹定義**：對配對樣本的差值符號秩做的非參數檢定。當 n 個樣本差值全同號時，檢定統計量達理論極值、p 達理論最小值（如 n=7 時 p=1/2⁷）；此 p 衡量「方向一致性」而非「效果大小」（各樣本 magnitude 可差很大）。
- **直覺類比**：7 個樣本每個都丟硬幣全正面的機率是 1/128 = p；它說「方向極一致」，不說「每個 gap 多大」。
- **具體例子**：跨樣本 phasing gap 全為正 → 極小 p，但 magnitude 跨樣本可變數倍——須分開報「方向一致性」與「effect size」。
- **常見誤解**：把極小 p 誤讀為強效果量。
- **關聯名詞**：Cohen ribbon、null（虛無分佈）、跨樣本復現、effect size

### 安全約束下 FP removal=0%（AUC 高但無可用閾值）🧩
- **嚴謹定義**：分類器 AUC 可顯示「能區分」，但在嚴格安全約束（如 TP loss ≤2%）下任何閾值都會同時移除過多 TP（因 TP/FP 分數分佈重疊），導致實際可移除的 FP 為 0%。即「可區分」與「可用過濾」是兩回事。
- **直覺類比**：要保留 98% TP，閾值必須極寬鬆，幾乎所有 FP 也一起留下。
- **具體例子**：read-level germline FP 模型 AUC 約 0.72，但 TP loss ≤2% 約束下 FP removal=0%。
- **常見誤解**：把 AUC>0.7 當「可用 filter」。
- **關聯名詞**：AUC、TP / FP、anti-discriminative

### anti-discriminative（反判別性）
- **嚴謹定義**：某特徵在欲判別類別（FP）中比對照類別（TP）富集，即用它「判別 TP」不但無效、反指錯方向（如 strong-ASM 在 FP 富集，OR<1）。
- **直覺類比**：用「血壓高」來篩健康人，反而篩到病人。
- **具體例子**：strong-ASM 在 FP 富集 → 越強的 ASM 越可能是 FP；但連續 |Δβ| 整體 AUC 仍近隨機（區分「子集 FP 富集」與「整體判別力」兩層）。
- **常見誤解**：把「strong-ASM 在 FP 富集」與「整體甲基訊號的判別力」混為一談。
- **關聯名詞**：strong-ASM、regression-to-extreme、AUC、ASM

### blind-ARI（盲化 Adjusted Rand Index）
- **嚴謹定義**：以 ARI（衡量兩種分群一致性的指標，校正隨機期望）評估「甲基化驅動的 read 分群」與「HP / somatic 結構」的對齊程度，並以 placebo（隨機標籤）為對照（blind）。高 blind-ARI 表示甲基分群確與目標結構對齊（真結構），非隨機。
- **直覺類比**：把甲基分群與 HP 分群兩張座位表比對，ARI 高代表兩表座位安排很像；blind 對照確保不是隨便都很像。
- **具體例子**：blind-ARI ≫ placebo → 確認甲基結構真實存在（即使它不能判別 TP/FP）。
- **常見誤解**：把 ARI 與 AUC 混用；ARI 量分群一致性、AUC 量二元判別。
- **關聯名詞**：null（虛無分佈）、PERMANOVA、credible ASM、germline-haplotype 層級甲基

### partial correlation（偏相關）
- **嚴謹定義**：在控制第三變數後，量化兩變數的線性關聯。本專案用以檢查 HP-axis |Δβ| 與 CN 在控制其他因素後的偏相關，若接近 0 或反向，支持「ASM 非 CN 驅動」。
- **直覺類比**：扣掉「天氣」這個共同因素後，問「冰淇淋銷量」與「溺水數」是否還相關。
- **具體例子**：HP-axis |Δβ| 與 CN 的偏相關反向且極小 → de-confounded 證據。
- **常見誤解**：把原始相關當因果，未控混淆。
- **關聯名詞**：Copy Number confound、HP-axis ASM、collider bias

### collider bias（對撞偏差 / collider）
- **嚴謹定義**：當對「共同果（collider）」做 residualize / 條件化時，會在其兩個因之間引入虛假關聯。本專案警示：residualize on AF 等可能產生虛假訊號，需 AF-bin 分層交叉驗證；且控制某變數後關聯反而增強可能是 collider 反向。
- **直覺類比**：只看「進了同一所名校的學生」，會誤以為「會考試」與「家世好」負相關（其實是被名校門檻 collider 篩出來的假象）。
- **具體例子**：AF→NGroups 控甲基後幾乎不動且方向反向，被判為 collider 反向（非真甲基效應）。
- **常見誤解**：把 residualization 一律當去 confound；對 collider 做條件化反而引入偏差。
- **關聯名詞**：partial correlation、pooled OLS 陷阱、AF-bin 交叉驗證

### within-group OLS vs pooled OLS（殘差化陷阱）
- **嚴謹定義**：pooled OLS residualization 會在殘差中保留分組資訊（產生虛假組間訊號），須改用 within-group OLS（組內回歸）並以 AF-bin 交叉驗證。是本專案反覆強調的方法紀律。
- **直覺類比**：把全班一起算殘差，會把「班別」資訊偷渡進殘差；應該各班分開算。
- **具體例子**：空間 auto-correlation（chr+pos 聚合特徵）的 AUC 須在 mid-TP-rate window 內驗證，否則純 artifact。
- **常見誤解**：以為任何 residualization 都能去 confound。
- **關聯名詞**：collider bias、partial correlation、spatial autocorrelation confound

### 結構（structure）vs 失序（disorder）🧩
- **嚴謹定義**：ISM 的核心區別——「結構」指 reads 在甲基模式空間中形成系統性分離的亞群（supervised PERMANOVA 可檢出、與 HP/somatic 標籤對齊）；「失序（disorder/entropy）」指 reads 甲基模式的隨機異質性（epipolymorphism/PDR/NME 量化）。ISM 量的是結構、把 disorder 當對照。
- **直覺類比**：「結構」像整齊分成兩疊的考卷（有規律），「失序」像每張都不一樣的亂碼（無規律）。
- **具體例子**：兩種情況的 per-CpG 均值可能相同（都 50%），但 PERMANOVA pseudo-F 一大一小——區分結構 vs 失序。
- **常見誤解**：把高 disorder 當有結構訊號。
- **關聯名詞**：PERMANOVA、epipolymorphism、軸 C、軸 E

### Grade B+（證據強度等級）
- **嚴謹定義**：本專案框架下的證據強度——某軸達 A-strength（如跨樣本 7/7 方向一致、Wilcoxon p 達理論下限）但決定性負控（如 R-SELFREF 全基因組 flag-on 重跑）未完成，故整體封頂 B+ 而非 A。Grade A 要求所有 self-ref 循環的負控均 PASS。
- **直覺類比**：申請大學——成績單 A（cross-sample 一致）已交，但推薦信（self-ref 負控）還差一封；錄取把握高（B+）但不到確定（A）。
- **具體例子**：LOH-phasing 脊柱目前為 Grade B+，是否投入負控跑完是論文關鍵待決議點。
- **常見誤解**：把方向一致性（A-strength 的一軸）等同整體 Grade A。
- **關聯名詞**：by-construction circularity、self-phasing circular dependency、tier ⭐

### tier ⭐（證據成熟度分級）
- **嚴謹定義**：研究結論成熟度的 ⭐1–5 標記。封頂規則之一：single-pipeline（共用同一 caller/phasing pipeline）無法排除共用偏差，故 tier 封頂 ⭐3；升 ⭐4 需獨立 pipeline / 額外樣本（如 COLO829 對照）。
- **直覺類比**：星級越高代表越能在不同條件下站得住；只跑一條 pipeline 像只在一家餐廳吃過，不能說「全城最好」。
- **具體例子**：跨 6 樣本 ASM 復現雖跨 3 癌種，但共用 ClairS-TO/longphase → 封頂 ⭐3。
- **常見誤解**：把樣本數多等同 tier 高；pipeline 獨立性才是封頂關鍵。
- **關聯名詞**：Grade B+、single-pipeline caveat、跨樣本復現

---

## ⑤ 工具與 Pipeline

### ISM（InterSubMod，Inter-Subclonal Methylation）
- **嚴謹定義**：本專案的 C++17 工具，輸入 haplotagged tumor BAM（含 MM/ML 甲基 tag + HP tag）、somatic SNV VCF、LOH BED、可選 normal BAM、參考基因組；以每個 somatic SNV 為中心開視窗，五階段處理（① read 擷取 → ② read×CpG 甲基矩陣 → ③ read-read 距離矩陣 → ④ UPGMA 聚類 → ⑤ 顯著性檢定套件）輸出 significance_summary.csv。核心定位為 read-level epigenetic characterization，非 variant filter。
- **直覺類比**：在每個突變附近把所有 read 攤成一張大表，算「哪些 read 甲基模式最像/最不像」，再統計驗證是否與腫瘤亞克隆結構有關。
- **具體例子**：ISM 坐落於上游（Dorado base/methyl calling + longphase haplotagging）與下游 Python 後處理（Δβ、normal-anchored cis-test、copy-partition）之間。
- **常見誤解**：把 ISM 當 variant filter；其誠實定位是 characterization + 去 confound。
- **關聯名詞**：軸 C、read×CpG 矩陣、PERMANOVA、significance_summary.csv、MSA

### 軸 A–F（業界甲基化方法 6 軸 MECE 分類）
- **嚴謹定義**：業界甲基化方法的 MECE 分類——A=per-position 率差（modkit dmr/DSS/methylKit，業界主流）；B=within-read CpG-CpG LD（分子內連動，MHB/MHL/Metheor）；C=between-read 距離+clustering（read 對 read，ISM 主場）；D=cohort 跨樣本網路；E=disorder/WSH scalar（epipolymorphism/PDR/NME）；F=長讀 phasing/ASM。ISM 主站軸 C、含軸 E 作對照、反向消費軸 F 的 HP tag。
- **直覺類比**：軸 B 是「同一份問卷裡不同題目的連動」（水平）；軸 C 是「不同問卷之間的整體差異」（垂直）；兩者都用 reads 但幾何軸正交。
- **具體例子**：沒有任何外部工具同時具備 ISM 的「顯式距離矩陣 + UPGMA/silhouette + PERMANOVA + normal-anchored cis-test + LOH/CN 耦合」五項組合（但組合無人佔 ≠ 更好/有用）。
- **常見誤解**：把軸 B（分子內 LD）與軸 C（分子間距離）混為「都是 read 層分析」。
- **關聯名詞**：軸 C、per-position 聚合、ISM、PERMANOVA、epipolymorphism

### 軸 C（between-read read-read 距離+clustering）
- **嚴謹定義**：以 read×CpG 矩陣為輸入，算每對 read 的距離（如 NHD）形成 N×N 顯式距離矩陣，再用階層聚類（UPGMA）或 EM 找亞群，並用 PERMANOVA 評估群間結構顯著性。核心問題：「reads 是否在高維甲基模式空間中形成結構上分離的雲」。是 ISM 的主場。
- **直覺類比**：把每條讀段當一份答卷，問「這些答卷分成幾個明顯不同群體」，而非「第 5 題平均答對率多少」（軸 A）。
- **具體例子**：軸 C 能在 per-CpG 均值看不到差異時揭示亞群結構（bimodal vs disorder）。
- **常見誤解**：以為軸 C 也是算率差；它是多變量結構檢定。
- **關聯名詞**：軸 A–F、read×CpG 矩陣、PERMANOVA、NHD、per-read 表示

### longphase-S / longphase-TO
- **嚴謹定義**：ONT 長讀 phasing/haplotagging 工具。longphase-S（paired mode）需 tumor+normal BAM，用 germline SNP 為錨點做精確 somatic haplotag（含 HP1-1/HP2-1）；longphase-TO（tumor-only mode）無 normal BAM、self-phasing，reads 系統性偏向 HP1（已知 bias）。
- **直覺類比**：longphase-S 是「有地圖的導航」（normal BAM 當地圖），longphase-TO 是「無地圖的導航」（只靠腫瘤自身 SNV 猜路，會系統性走錯）。
- **具體例子**：ISM 需要 longphase 提供 HP tag（尤其 somatic sub-haplotype HP1-1）才能做 HP-axis 分析；canonical benchmark 用 paired（longphase-S）。
- **常見誤解**：以為 TO 與 paired 的 HP tag 等價可比；TO self-phasing 有循環偏差。
- **關聯名詞**：somatic haplotagging、HP1-1、self-phasing circular dependency、Tumor-only (TO) pipeline

### ClairS / ClairS-TO（somatic variant caller）
- **嚴謹定義**：ONT 長讀 somatic 變異 caller。ClairS（paired）用 tumor+normal 精確分離 germline/somatic；ClairS-TO（tumor-only）依賴 Panel of Normals（PON）過濾 germline，germline-somatic 分離不完全，LOH annotation 殘缺風險較高。
- **直覺類比**：ClairS-TO 像沒有健康對照組的疾病篩查——能偵測異常，但不確定哪些是先天（germline）哪些後天（somatic）。
- **具體例子**：ISM 的 somatic SNV VCF 來自 ClairS/ClairS-TO；跨 6 樣本 ASM 研究共用 ClairS-TO → single-pipeline 封頂 ⭐3。
- **常見誤解**：把 ClairS-TO 輸出全當乾淨 somatic。
- **關聯名詞**：Somatic mutation vs Germline mutation、Panel of Normals、single-pipeline caveat、TP / FP

### Panel of Normals (PON)
- **嚴謹定義**：由一組正常樣本建立的變異資料庫，用以在 tumor-only caller（如 ClairS-TO）中過濾常見 germline / 技術假陽性變異。
- **直覺類比**：一份「正常人常見變異」黑名單，看到名單上的就先剔除。
- **具體例子**：ClairS-TO 在無 matched normal 時靠 PON 近似 germline 過濾。
- **常見誤解**：以為 PON 等同 matched normal；PON 是族群層近似，無法完全取代個體 normal。
- **關聯名詞**：ClairS-TO、Somatic mutation vs Germline mutation、Tumor-only (TO) pipeline

### modkit（Oxford Nanopore modkit）
- **嚴謹定義**：ONT 官方工具，從 MM/ML tag 的 BAM 提取 base modification（5mC/5hmC），輸出 per-site 甲基率（bedMethyl）。操作在 bulk pileup 層（軸 A 率差）。其 dmr 內部無 --haplotype flag，需先 pileup --phased 拆 hp1/hp2 再餵 dmr pair（等同把 haplotype 當靜態樣本，失去分子連結）。
- **直覺類比**：modkit 是「計算每個位置平均被螢光筆畫幾次」，ISM 是「比較兩條 DNA 分子的螢光筆標記模式有多不像」——層次不同。
- **具體例子**：可用 modkit 與 ISM 的甲基提取做一致性外部驗證（高相關），確認 ISM 甲基讀取正確。
- **常見誤解**：以為 modkit dmr 能做 somatic-controlled HP1 vs HP1-1 軸；它失去分子連結，無此能力。
- **關聯名詞**：per-position 聚合、軸 A–F、MM/ML tag、DSS

### DSS（Differential methylation via beta-binomial）
- **嚴謹定義**：偵測差異甲基位點/區域（DMS/DMR）的 R 套件，用 beta-binomial 建模 CpG 甲基的 over-dispersion + lognormal dispersion shrinkage + Wald，比 Fisher 更能控假陽性。操作在 pileup 層（軸 A），適群體層差異甲基，非 read-level 亞克隆分析。是業界統計黃金標準。
- **直覺類比**：DSS 是「統計學家找兩組人平均心率差異的方法」，ISM 是「醫師看每個病人心電圖找異質模式」。
- **具體例子**：建議 ISM 的 per-CpG Fisher 學 DSS 的 beta-binomial + shrinkage。
- **常見誤解**：把 DSS 的群體層 DMR 能力當 read-level 亞群分析。
- **關聯名詞**：beta-binomial、Fisher exact test、軸 A–F、over-dispersion

### MSA（MethylSomaticAnalysis，外部抽取工具）
- **嚴謹定義**：外部甲基抽取工具，≠ ISM binary。其抽取合理但 Level-3 統計不可用（比較 VCF 來源非 HP group），改用 Level-1 + Python。已知 dual-row bug：對雙修飾 BAM 每 (read,CpG) 發 5mC+5hmC 兩列分別計數，使 β 砍半（max-collapse 修正後與 5mC-only 收斂）。
- **直覺類比**：MSA 像一個「把甲基原料抽出來」的前處理器，但它附帶的統計報表（Level-3）不能直接信，要自己用 Python 重算。
- **具體例子**：早期某位點 Δβ 偏小是 MSA Level1 雙列砍半 artifact，max-collapse 修正後翻倍收斂。
- **常見誤解**：把 MSA 當 ISM binary；或直接信 MSA Level-3 統計。
- **關聯名詞**：max-collapse、5mC vs 5hmC、ISM、Δβ

### cvlr / ASMS / qFDRP（軸 C 近鄰工具）
- **嚴謹定義**：與 ISM 最相近的軸 C 工具。cvlr = Bernoulli-mixture EM（soft posterior，無顯式距離矩陣、無顯著性檢定、固定 k）；ASMS（Raineri 2024 preprint，⚠UNVERIFIED）= ONT no-phasing EM 2-component，最像但未 peer-review；qFDRP = 與 ISM 用相同 normalized-Hamming kernel 但塌成 per-CpG scalar，丟掉幾何。
- **直覺類比**：cvlr/ASMS 用「機率模型分群」（soft），ISM 用「顯式距離矩陣 + 顯著性檢定」；qFDRP 把矩陣壓成一個數字。
- **具體例子**：沒有任何單一工具同時具備 ISM 的五項組合（距離矩陣 + UPGMA/silhouette + PERMANOVA + normal-anchored cis + LOH/CN 耦合）。
- **常見誤解**：把「最像」當「等價」；各工具缺 ISM 的某項關鍵能力。
- **關聯名詞**：軸 C、Bernoulli-mixture EM、PERMANOVA、normal-anchored cis-test

### DAMEfinder（per-CpG ASM 參考工具）
- **嚴謹定義**：per-CpG / CpG-pair tuple 的 ASM 偵測工具。⚠ 引用訂正：ISM 源碼曾將作者誤標「De Waele 2020」，正確為 Orjuela et al. 2020（Epigenetics & Chromatin 13:25, DOI 10.1186/s13072-020-00346-8, PMC7268773）。
- **直覺類比**：一個專門找「某 CpG 或某對 CpG 在兩等位上甲基不同」的偵測器。
- **具體例子**：建議 ISM 借鑑 DAMEfinder 的 CpG-pair tuple log-odds（within-read 連動）。
- **常見誤解**：沿用舊引用「De Waele 2020」（捏造作者，已訂正為 Orjuela 2020）。
- **關聯名詞**：軸 B、ASM、epipolymorphism、引用訂正

### Pipeline 三條（paired_full / paired_pileup / tumor_only）🧩
- **嚴謹定義**：ISM 的三條評估 pipeline。paired_full = canonical benchmark（需 normal BAM，longphase-S phasing，ClairS paired），結果正向；paired_pileup = 副次評估；tumor_only (TO) = 無 normal BAM、longphase-TO self-phasing、ClairS-TO，因 self-phasing bias 結果 NEGATIVE（已凍結）。⚠ 三條的 F1 口徑不同，不可直接互比。
- **直覺類比**：三套不同規則的考試，分數不能直接比高低，要看用哪套規則。
- **具體例子**：canonical benchmark 選 paired_full 而非 TO，因 TO 有 self-phasing 循環偏差。
- **常見誤解**：把三條 pipeline 的 ΔF1 並列比較。
- **關聯名詞**：Tumor-only (TO) pipeline、self-phasing circular dependency、longphase-S/-TO、口徑差（provenance）

### Tumor-only (TO) pipeline
- **嚴謹定義**：無 matched normal BAM 的 pipeline——longphase-TO self-phasing + ClairS-TO（PON 過濾）。因 self-phasing circular dependency 產生虛假 LOH、haplotag 系統性偏向 HP1，結果不可靠（ΔF1 NEGATIVE，Phase 1A 凍結）。
- **直覺類比**：沒有「正常對照組」的整套分析，導航與判讀都會系統性走偏。
- **具體例子**：TO 模式所有下游問題（QS AUC 近隨機、LOH over-calling、甲基方向反轉）的根源是 self-phasing。
- **常見誤解**：以為 TO 只是「少一個樣本」；它的循環偏差是結構性缺陷。
- **關聯名詞**：self-phasing circular dependency、longphase-S/-TO、ClairS-TO、Pipeline 三條

### significance_summary.csv（ISM 輸出）
- **嚴謹定義**：ISM 對全位點的滙總輸出（每 region 一列，含 Significant gate、HPMergedDelta、CramérV、Fisher_Frac_Sig、PERMANOVA p、HP_Residual_Delta 等多欄）。欄位數依版本（曾記 59 / 117 欄，需以當前版本為準）。是下游 Python 後處理（F1、AUC、cis-test）的輸入。
- **直覺類比**：每個突變位點一行的「體檢報告表」，匯總各項統計指標。
- **具體例子**：display 頁 / aggregate 由此 CSV 衍生；CramérV 欄已 Cochran 閘控（與 significance.json raw 值雙口徑）。
- **常見誤解**：把 summary 的 gated CramérV 與 significance.json 的 raw 值混用。
- **關聯名詞**：ISM、Cramér's V reliability gate、HP_Residual、PERMANOVA

### single-pipeline caveat（共用 pipeline 偏差）
- **嚴謹定義**：當多樣本/多分析共用同一條 pipeline（同一 caller + 同一 phasing 工具）時，pipeline 特有偏差無法被排除，故結論 tier 封頂 ⭐3；升 ⭐4 需獨立 pipeline 或正交驗證。
- **直覺類比**：所有菜都用同一個壞掉的秤量過——量出來一致，但可能一致地偏。
- **具體例子**：跨 6 樣本 ASM 雖跨 3 癌種復現，但共用 ClairS-TO/longphase → 封頂 ⭐3。
- **常見誤解**：把跨樣本一致當已排除 pipeline 偏差。
- **關聯名詞**：tier ⭐、Grade B+、跨樣本復現、Pipeline 三條

### 口徑差（provenance / data tiering）🧩
- **嚴謹定義**：同一現象在不同 build / 計算口徑（如 Tmode vs paired、5mC-only vs max-collapse、read-instance vs unique、gated vs raw CramérV）下的數值不可直接並列；每筆數據須隨時標 tier（原始可信 P / 方法存疑 P-caveat / 二次紀錄 S / 已捏造 F），引用前查級。
- **直覺類比**：用公斤與磅量同一袋米，數字不同不代表米變了；先確認單位（口徑）再比較。
- **具體例子**：MSA Level1 5mC+5hmC 雙列 vs max-collapse、summary gated CramérV vs significance.json raw、三條 pipeline 的 F1——皆為不可並列的口徑差。
- **常見誤解**：把不同口徑的數字直接比較或並列。
- **關聯名詞**：max-collapse、Cramér's V reliability gate、Pipeline 三條、tier ⭐

---

## 附錄：建議配 SVG 的高認知負荷名詞速查（🧩）

下列名詞讀者直覺易誤、建議在解釋頁配 SVG 示意圖：

| 名詞 | 為何高認知負荷 | SVG 建議 |
|---|---|---|
| HP-tag | 多態 tag 易誤為染色體編號 | integer→語意對照表（1/2/11/21/33/0） |
| HP1-1（somatic sub-haplotype） | 易誤為第三染色體 / copy 標記 | 染色體鏈 + somatic ALT 分裂出 HP1-1 |
| by-construction circularity (R-SELFREF) | 循環隱蔽難察 | HP1-1 定義 → Inner bucket 篩選的閉環 + flag-on 打破循環 |
| self-phasing circular dependency | anchor 與被評估對象同一 | TO 模式 ALT reads 既是錨點又被評估的環形圖 |
| Δβ vs Δ | 兩個都叫「差」卻是不同空間 | 左 N×N 距離熱圖（Δ）/ 右 per-CpG β 條形（Δβ） |
| HP-axis vs ALLELE-axis | 切法不同 → confound 程度不同 | germline HP1/HP2 + somatic HP1-1 三層，標哪軸 held-constant |
| normal-anchored cis-test | 三 delta 相互依存 | 三盒（normal-HP1/tumor-HP1/tumor-HP1-1）+ d_drift/d_somatic/d_cis 箭頭 |
| Copy Number confound | dosage 假象難直觀 | copy 多 → reads 多 → 表觀等位差 vs 真調控 |
| LOH-coverage-baseline 三角 | 三維互相關難拆 | LOH→單倍型→低覆蓋→極端 baseline 因果鏈 |
| Cramér's V reliability gate / latent 真結構 | gated vs raw 雙口徑、CV=0 但 PERMANOVA 顯著 | 2×2（CV reliable/gated × PERMANOVA 顯著/不顯著）標 latent 格 |
| PERMANOVA（結構 vs 失序） | 「結構」非「率差」 | bimodal（有結構）vs disorder（無結構）兩矩陣，均值相同 pseudo-F 不同 |
| Fisher over-dispersion | 跨「統計假設」+「長讀生物學」兩域 | 同 clone reads 非獨立 + 2×2 表 vs beta-binomial φ |
| LOSO / sample-level circularity | k-fold 易誤為跨樣本驗證 | within-sample 圓內切 vs LOSO 跨圓箭頭 |
| anchor AUC 膨脹 | 絕對 AUC 系統性偏高 | real AUC 分佈 vs shuffle null 分佈，標 delta 才是 effect |
| AUC 高但 FP removal=0% | 可區分 ≠ 可用閾值 | ROC + TPR≥0.98 橫線對應 FPR≈1.0 |
| Wilcoxon 方向一致性 ≠ effect | 極小 p 易誤讀為強效果 | 7 樣本 gap 全正但高低不一 + null 直方圖標極值位置 |
| methyl-assisted phasing 機制天花板 | between vs within haplotype 分層 | 樹狀圖標甲基在哪層可分（HP 間強、亞群弱）+ unphase 救援漏斗 |
| per-position vs per-read | 均值丟失共甲基化結構 | 4 read×5 CpG 矩陣，均值相同但距離矩陣分群不同 |
| Pipeline 三條 / 口徑差 | F1 不可並列 | 三 pipeline 並排 + 「口徑不同不可直接比」警示 |
