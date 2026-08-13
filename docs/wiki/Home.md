# InterSubMod Wiki

**Local mutation-state candidate reconstruction from single-molecule sSNV co-occurrence。**

本專案重建的是**局部、允許遞迴、最小的突變狀態候選 arborescence／topology**；它不是已確認的細胞 clone、subclone 或腫瘤譜系。

這個 Wiki 是專案的完整說明：**有哪些部件、各自吃什麼吐什麼、現在哪些真的跑得動、以及怎麼動手跑。**
給指導教授、實驗室同學，以及任何第一次接觸這套系統的人。

---

## 先看這張圖

![系統架構](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/a34b0cb96a8ef247c5a6f423d46b2920c7e541aa/docs/images/architecture-overview.png)

由下而上五層。方框顏色代表現在的可用狀態：綠＝可跑、黃＝有限制、紅＝被鎖住。

---

## 這個專案在解決什麼問題

一顆腫瘤不是單一細胞群，而是好幾群帶著不同突變組合的細胞混在一起。
要理解抗藥性與病程發展，就需要知道**哪個突變先發生**、**哪些突變住在同一群細胞裡**。

短讀長測序只看得到每個位點各自的變異等位頻率；若只靠 per-locus marginal VAF、
且沒有 linkage 或額外模型假設，聯合結構在此逆問題下不可識別。

ONT 長讀長改變了可觀測資料：一條分子可以同時跨過好幾個體細胞突變，
因此可以**直接觀測同一條物理分子上的共現**。兩個突變是否屬於同一細胞、哪個先發生或是否形成細胞譜系，仍是受模型與資料限制的推論，不能由一條 read 直接確認。

---

## 三分鐘掌握重點

| # | 重點 | 一句話 |
|---|---|---|
| 1 | **骨幹是什麼** | 同一條物理分子上的候選體細胞 allele 共現。只在不使用甲基化衍生標籤的限定下非循環；仍依賴 variant calling、alignment、basecalling 與 haplotag 假設。 |
| 2 | **甲基化的角色** | **bounded-auxiliary（有界輔助）**。候選突變狀態拓撲定好之後才計算，只做 association-only 註記，**動不了任何一條邊**。 |
| 3 | **為什麼不能用甲基化重建** | 單一 bulk 無法區分四種成因（germline ASM / LOH 解遮蔽 / 拷貝數劑量 / 真譜系差異）→ 循環論證。 |
| 4 | **實際產出** | 7 個 technical datasets／6 個 biological IDs、chr1–22：frozen model 對 63,506 / 71,955 個 rankable candidate units 指派單一 rooted-unlabeled candidate-shape signature（88.2579%）；不是 biological topology、accuracy 或 prevalence。 |
| 5 | **能力天花板** | 在目前單一 bulk observation／model 且未整合 CN／LOH 下，**確認的細胞亞群 = 0、確認的線性祖先關係 = 0**。single-cell、multi-region、orthogonal CN／purity 等獨立證據可能提高識別性；不宣稱只有某一種方法可行。 |
| 6 | **哪條線產生 exact-PS funnel** | 獨立的 research `exact_ps_topology_af` C++ solver ＋ Python runners。`inter_sub_mod` 本身產生 per-region 甲基化／統計輸出；**不是** LongLineage 主線（見下）。 |

---

## 各頁導覽

| 頁面 | 你會得到什麼 |
|---|---|
| **[System Overview 系統全景](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview)** | 五層架構、兩 repo 關係、**誠實狀態表**、7 technical datasets／6 biological IDs funnel、能與不能回答什麼 |
| **[InterSubMod Engine](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine)** | 3 個必填輸入、8 個內部階段、17 種輸出檔（附真實 header）、3 條實跑指令、9 條陷阱 |
| **[LongLineage Engine](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine)** | PRIVATE research-preview candidate `b9aaa12` 的子命令狀態、M1→M2→topology 鏈與 artefact 契約；frozen HCC1395 receipt 為 **0 candidate topology units**，不是「公開版解出 0 棵譜系樹」 |
| **[Upstream & Data](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data)** | Dorado / ClairS / LongPhase-S / SAVANA、sidecar 格式、7 technical datasets／6 biological IDs、三個算錯數字的陷阱 |
| **[Analysis & Presentation](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation)** | 該用哪些 Python 腳本、拒絕渲染如何防止 spec 宣告的必填欄位被靜默省略（不驗證來源真實、spec 完整或科學正確）、這層的四個坑 |
| **[How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)** | 六個步驟，每步附驗收條件與真實輸出；常見狀況排除表 |

---

## 最快的上手路徑

```bash
git clone https://github.com/liaoyoyo/InterSubMod.git && cd InterSubMod
HANDOFF_COMMIT="<IMMUTABLE_HANDOFF_COMMIT_SHA>"
git checkout --detach "$HANDOFF_COMMIT"
test "$(git rev-parse HEAD)" = "$HANDOFF_COMMIT"
test -z "$(git status --porcelain)"

# 1. repo 外 clean build；build output 不進版本控制
REPO_ROOT="$(pwd -P)"
BUILD_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/ism-build.XXXXXXXX")"
cmake -S "$REPO_ROOT" -B "$BUILD_ROOT" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_ROOT" -j$(nproc)
test -z "$(git -C "$REPO_ROOT" status --porcelain)"

# 2. 驗證   -> 以本次 CTest/run_tests 輸出動態取得 test/suite count，且 failure=0
"$BUILD_ROOT/bin/run_tests"
ctest --test-dir "$BUILD_ROOT" --output-on-failure
```

完整六步驟與每步的驗收條件見 **[How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)**。

---

## 一定要先知道的三件事

> **① `tree.nwk` 的葉子是 read，不是 clone。**
> 它是「read 依甲基化相似度」的階層分群樹，**不是**亞群演化譜系樹。這是最常被誤讀的輸出。

> **② 88.26% 不代表「腫瘤演化史解出 88%」。**
> 它的分母是「已經可排序的 71,955 個單元」。另有 **170,131 / 255,752（66.52%）個 strict read-linkage components 為 k=1**；66.52% 的分母是 strict components，不能說成「全部突變的三分之二」。
> 88.26% 是 local、recurrence-allowed、model-conditional 的數學圖形統計。

> **③ LongLineage 的「BLOCKED」是具名 gates，不是單一 parity 總結。**
> `P3/P4/P5/P7/P8` 尚未通過，範圍涵蓋 parity／validation、source-origin、license、dependency 與 release-safety。
> 部分 M1／M2／topology kernels 曾由具名 dataset-gate 執行，不等於全部核心、production entry 或公開資格完成。

---

## 這些內容是怎麼來的

Wiki 的 editorial 初稿源自 2026-08-06 對兩個 repo 的系統性實測收集；科學數值權威固定為
2026-08-01 authority bundle，公開 claim inventory 鎖於 2026-08-12，本輪 source corrections
日期為 2026-08-13。這些時間層不可合併成「同一次全部驗證」：

- **歷史自述、目前 `UNVERIFIED`（ALG-023）**：2026-08-06 文件曾記錄 334 個 source refs
  ＋111 個 paths 為 0 missing／0 out-of-bounds；公開 repo 缺 commit-pinned inventory、commands、
  hashes 與 replay receipt，因此不能把它當成 current blanket guarantee
- 標示「可跑」的部件，都是**實際執行並檢查 exit code** 的結果
- funnel 各層數字取自凍結的 canonical 輸出，且**已驗證各層加總自洽**

同樣的內容另有排版更豐富的 **standalone HTML 版本**（含可展開的細節區塊）位於
repo 的 `docs/explain/`，clone 後可離線開啟。

### 誠實標註的已知缺口

- 「純 parsimony 單一拓撲率」的**下界 64.89%** 在 repo 內找不到原始出處，對外引用前需補齊（上界 88.26% 有佐證）
- LongLineage 對 7 technical datasets／6 biological IDs cohort 的**執行時間與記憶體上界從未實測**，其自身文件明文禁止由單一 dataset 外推
- 拷貝數（CN／LOH）目前狀態為 `NOT_INTEGRATED`，尚未接入主線 —— 因此現有結論皆為「未經拷貝數校正」
