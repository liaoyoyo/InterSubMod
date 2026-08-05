# V2 · unphase / HP3 分樹決策 — 驗證 digest

**軌**: V2 (unphase/HP3 tree-splitting decision)
**樣本**: HCC1395 單樣本 (⭐3 上限), tagged tumor BAM
**狀態**: ✅ COMPLETED — 決策有實際 BAM 掃描數字支持
**日期**: 2026-07-04

---

## 一句話結論

> **HP3 read 應併入「帶對應 somatic 突變的 lineage 樹」,不另立第三棵;unphase read 是定相失敗背景,不進任何 somatic lineage 樹(獨立/丟棄)。**

三條獨立證據軸一致指向此結論,且與 PI 2026-07-04 修正的 HP-tag 語意(判別軸=是否經過 ALT)完全吻合。

---

## 決策依據:unphase vs HP3 招募力對比表

全部數字來自實際 BAM 掃描 (`docs/methodology/_assets/20260704_V_verification/V2/sm_read_recruitment_v2_split.json` + `_perregion.tsv`),
掃描 754 個 LOH/neutral 目標區、51 個具 germline-ASM 甲基錨點的區。

| 指標 | **HP3**(HP tag=3) | **unphase**(無 HP tag) | 判讀 |
|---|---|---|---|
| **ALT-carriage(定義軸)**anchored | **112/112 = 100%** | 3/4386 = **0.07%** | HP3 定義=經過 ALT;unphase 幾乎不帶 somatic 訊號 |
| ALT-carriage 全 754 區 | 1753/1798 = **97.5%** | 57/27246 = **0.21%** | genome-wide 同結論 |
| 對照:HP1 / HP2 ALT-rate | HP1 43%, HP2 44%(全區) | — | HP1/HP2 為混合(germline 家族含帶/不帶突變的 read) |
| **可招募池**(覆蓋≥3 錨點 CpG) | 16 | 318 | unphase 池大(定相失敗多),HP3 池小(reads 少+多區無 germline 錨) |
| 可招募且帶 ALT | **16/16 = 100%** | **0/318 = 0%** | 招募到的 HP3 全帶突變;unphase 全不帶 |
| **招募指派**(最近甲基質心→HP1/HP2) | →HP1=**15**, →HP2=**1**(15:1) | →HP1=161, →HP2=157(≈50/50) | HP3 強偏一條 lineage;unphase 是隨機背景 |
| germline-HP 招募力(HP1/HP2 LOO) | overall **96.3%**(42 區, 5147 reads;per-region median 98.6%) | 同一基線 | 甲基錨點能到 haplotype 層(既有 96.3% 細分,已複現) |
| within-HP subclone 甲基可分性 | HP3 內: **0/0 區**(無 HP3 多基因型子群可測) | HP1/HP2 內: **2/19 = 10.5%** | 甲基在 subclone 層招募力弱(既有 10.5% 細分,已複現) |

### 交叉驗證(證明重寫忠實)
拆分後各池加總 = 原始合併數字,逐項吻合:
- unphase(4386) + HP3(112) = 4498 = 原 `total_hp0_in_anchored` ✔
- HP1(7568) + HP2(5030) = 12598 = 原 `total_hp_tagged_in_anchored` ✔
- 可招募 318 + 16 = **334** = 原 `total_unphased_hp3_recruitable` ✔
- LOO 42 區 / 5147 reads / 4956 correct / **0.9629** — 與原始逐位元相同 ✔

---

## 分樹決策(逐條)

### 1) HP3 → **併入帶該突變的 lineage 樹**(不另立第三棵)
- **定義證據**: HP3 的判別軸就是「經過 ALT」— anchored 100%、全區 97.5% 帶 somatic ALT。HP3 read *本身攜帶* somatic 訊號,不是無訊號的雜項。
- **歸屬證據**: 16 個可招募 HP3 中 15 個甲基質心最近 HP1、1 個最近 HP2(15:1),明確偏向*特定一條* germline-haplotype lineage,而非等距於兩者。
- **甲基定相證據**(既有 `h3_methyl_phasing.json`,93 區): 10/10 可測 H3 區判為「近HP(定相失敗)」,0 真第三群 — HP3 甲基輪廓總是貼近既有 HP1 或 HP2,從不自成一群。
- **含意**: HP3 = 「帶 ALT 但 germline 定相未確認/HP 衝突」的 read。它屬於某條已知 lineage(帶該 ALT 的那條),應被招募進去,**不構成新的克隆分支**。單樣本無法對 1317 個 hp3-ungateable somatic 逐一 germline-gate(無錨),故此為**候選**級結論:方向確定(併樹),但每條 HP3 read 具體歸哪棵,需 normal-same-HP 基線或更深錨點逐區確認。

### 2) unphase → **不進任何 somatic lineage 樹**(定相失敗背景;獨立看待或丟棄)
- ALT-carriage ≈ 0%(anchored 0.07%、全區 0.21%): unphase read 幾乎不帶 somatic 突變 → 無 lineage 訊號。
- 招募指派 ≈ 50/50(161:157): 甲基把它們隨機丟向 HP1/HP2,無方向性 → 確認是背景,不是隱藏的第三群。
- **含意**: unphase 是「REF-only、germline 未碰或衝突」的定相失敗 read。可作背景池(甲基招募回 HP 層做覆蓋補強),但**不應以 unphase 身分獨立成一棵 somatic 樹**,也不攜帶可定義新分支的突變。

### 為何不設「第三棵 HP3 樹」
若 HP3 各自獨立成樹,等於假設 HP3 是與 HP1/HP2 平行的第三個 germline haplotype 群 —
但 (a) 甲基輪廓 0/10 支持獨立群,(b) 招募 15:1 偏向既有 HP,(c) HP3 定義本就是「帶 ALT 但定相未定」,
其突變本身即指向某條既有 lineage。第三棵樹會把同一 lineage 的證據拆散、虛增分支。

---

## 方法 / 檔案

- **重寫腳本**: `docs/methodology/_assets/20260704_V_verification/V2/sm_read_recruitment_v2_split.py`
  (改自 `docs/methodology/_assets/20260618_subcluster_pilot/scripts/sm_read_recruitment.py`)
- **核心修正**: 原 `hp_germ()` 把 HP tag `"3"` 與無 tag(unphase)都映射成 `0`,`recruitable_hp0=334` 混淆兩者。
  本版 `hp_class()` 拆成 hp1/hp2/**hp3**/**unphase** 四類,並對每條 read 記錄 somatic-ALT carriage(PI 判別軸)。
- **輸出**: `docs/methodology/_assets/20260704_V_verification/V2/sm_read_recruitment_v2_split.json`, `..._perregion.tsv`(754 區)
- **既有輸入**: `h3_methyl_phasing.json`(93 H3 區甲基定相)、`sm_hp_contribution.json`(hp3_ungateable_somatic=1317)
- **執行環境**: `whatshap` conda env (pysam 0.22.1 / scipy 1.14.1 / numpy 2.1.1), BAM = HCC1395 tagged tumor, MAPQ≥20, dbeta≥0.2, anchor_min=3。掃描 ~28 min 單執行緒。
- **甲基定位**: 全程 Tier-3,不 override genetic/HP;此分析是甲基*招募力測試*,不改變 genetic 定樹。

## 限界(⭐3 單樣本,不 overclaim)
- **已驗證**: 拆分數字忠實(交叉核對通過);ALT-carriage 對比(100% vs 0.07%)為定義層事實;招募指派方向(15:1 vs 50/50)。
- **候選/待驗**: HP3 併入*哪一棵*具體 lineage 的逐 read 指派(需 germline 錨/normal-same-HP 基線;1317 hp3-ungateable 無法單樣本逐一確認);跨樣本重現性。
