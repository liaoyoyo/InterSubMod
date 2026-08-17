# InterSubMod Wiki

**Local mutation-state candidate reconstruction from single-molecule sSNV co-occurrence。**

本專案重建的是**局部、允許遞迴、最小的突變狀態候選 arborescence／topology**；它不是已確認的細胞 clone、subclone 或腫瘤譜系。

這個 Wiki 是專案的完整說明：**有哪些部件、各自吃什麼吐什麼、現在哪些真的跑得動、以及怎麼動手跑。**
給指導教授、實驗室同學，以及任何第一次接觸這套系統的人。

---

## 先看這張圖

![系統架構](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/develop/docs/images/architecture-overview.png)

由下而上五層。方框顏色代表現在的可用狀態：綠＝可跑、黃＝有限制、紅＝被鎖住。

---

## 這個專案在解決什麼問題

一顆腫瘤不是單一細胞群，而是好幾群帶著不同突變組合的細胞混在一起。
要理解抗藥性與病程發展，就需要知道**哪個突變先發生**、**哪些突變住在同一群細胞裡**。

若輸入只有各位點的邊際變異等位頻率，沒有 linkage 或額外模型假設，可能有多個聯合結構
產生相同邊際值，因此無法只從這些邊際值識別唯一聯合結構。

ONT 長讀長改變了可觀測資料：一條分子可以同時跨過好幾個體細胞突變，
因此可以**直接觀測同一條物理分子上的共現**。兩個突變是否屬於同一細胞、哪個先發生或是否形成細胞譜系，仍是受模型與資料限制的推論，不能由一條 read 直接確認。

---

## 三分鐘掌握重點

| # | 重點 | 一句話 |
|---|---|---|
| 1 | **骨幹是什麼** | 同一條物理分子上的體細胞突變共現。不依賴任何待推論的標籤，因此非循環。 |
| 2 | **甲基化的角色** | **bounded-auxiliary（有界輔助）**。候選突變狀態拓撲定好之後才計算，只做 association-only 註記，**動不了任何一條邊**。 |
| 3 | **為什麼不能用甲基化重建** | 以目前 single-bulk measurement set、且沒有 orthogonal data 或額外假設時，四種成因（germline ASM / LOH 解遮蔽 / 拷貝數劑量 / 真譜系差異）不可識別；拿同一訊號自我確認會形成循環論證。 |
| 4 | **實際產出** | 7 個資料集、chr1–22：71,955 個可排序且 family-complete 單元中，63,506 個（88.26%）只有一種 rooted-unlabeled **數學拓撲**；不是細胞譜系普及率。 |
| 5 | **能力天花板** | 單一 bulk 只能 characterize 不能 confirm。**確認的細胞亞群 = 0**。這是資訊論界限，不是實作缺口。 |
| 6 | **哪條線產生 exact-PS funnel** | 獨立的 research `exact_ps_topology_af` C++ solver ＋ Python runners。`inter_sub_mod` 本身產生 per-region 甲基化／統計輸出；**不是** LongLineage 主線（見下）。 |

---

## 各頁導覽

| 頁面 | 你會得到什麼 |
|---|---|
| **[System Overview 系統全景](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview)** | 五層架構、兩 repo 關係、**誠實狀態表**、全 7 樣本 funnel、能與不能回答什麼 |
| **[InterSubMod Engine](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine)** | 3 個必填輸入、8 個內部階段、17 種輸出檔（附真實 header）、3 條實跑指令、9 條陷阱 |
| **[LongLineage Engine](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine)** | 4 個子命令狀態、M1→M2→topology 鏈、artefact 契約、**為什麼輸出 0 棵樹** |
| **[Upstream & Data](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data)** | Dorado / ClairS / LongPhase-S / SAVANA、sidecar 格式、7 樣本、三個算錯數字的陷阱 |
| **[Analysis & Presentation](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation)** | 該用哪些 Python 腳本、缺少已宣告必填欄位時的 fail-closed 設計、這層的四個坑 |
| **[How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)** | 六個步驟，每步附驗收條件與真實輸出；常見狀況排除表 |

---

## 最快的上手路徑

```bash
git clone https://github.com/liaoyoyo/InterSubMod.git && cd InterSubMod

# 1. 編譯（repo 內的執行檔是 STALE 的，務必重新編譯）
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)

# 2. 驗證   -> 以此 commit 實際輸出與退出碼為準，不比對硬編測試數
./build/bin/run_tests
```

完整六步驟與每步的驗收條件見 **[How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)**。

---

## 一定要先知道的三件事

> **① `tree.nwk` 的葉子是 read，不是 clone。**
> 它是「read 依甲基化相似度」的階層分群樹，**不是**亞群演化譜系樹。這是最常被誤讀的輸出。

> **② 88.26% 不代表「腫瘤演化史解出 88%」。**
> 它的分母是「已經可排序的 71,955 個單元」。另有 **170,131 / 255,752（66.52%）個 strict read-linkage components 為 k=1**；66.52% 的分母是 strict components，不能說成「全部突變的三分之二」。
> 88.26% 是 local、recurrence-allowed、model-conditional 的數學圖形統計。

> **③ LongLineage 的「BLOCKED」不是程式沒寫。**
> 指的是對照驗證的證據尚未存在。M1／M2／topology 的核心都已實作，也被實際執行過。
> 這是刻意的 fail-closed 設計。

---

## 這些內容是怎麼來的

Wiki 的所有內容源自 2026-08-06 對兩個 repo 的一次系統性實測收集，
再由獨立的對抗式驗證檢查過：

- **334 個「檔案:行號」宣稱 + 111 個路徑宣稱**經機械重驗 → **0 捏造、0 行號越界**
- 標示「可跑」的部件，都是**實際執行並檢查 exit code** 的結果
- funnel 各層數字取自凍結的 canonical 輸出，且**已驗證各層加總自洽**

同樣的內容另有排版更豐富的 **standalone HTML 版本**（含可展開的細節區塊）位於
repo 的 `docs/explain/`，clone 後可離線開啟。

### 誠實標註的已知缺口

- 「純 parsimony 單一拓撲率」的**下界 64.89%** 在 repo 內找不到原始出處，對外引用前需補齊（上界 88.26% 有佐證）
- LongLineage 的 **7 樣本執行時間與記憶體上界從未實測**，其自身文件明文禁止由單一樣本外推
- 拷貝數（CN／LOH）目前狀態為 `NOT_INTEGRATED`，尚未接入主線 —— 因此現有結論皆為「未經拷貝數校正」
