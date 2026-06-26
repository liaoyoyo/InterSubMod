<!--
建立時間: 2026-06-27（2026-06-27 修正 adversarial F1: baseline 0.443→0.500 gateable-only）
類型: L1 HP-tag 貢獻 — 單變量 ablation（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3）
data_sources: data/sm_hp_contribution.json
-->

# L1 — longphase-S HP tag 的貢獻（單變量驗證；已修對抗稽核）

> 問題：HP tag 對 subclone 判別有沒有 data-supported 的貢獻？做 ablation + **正確的 chance baseline** 對照。
> 圖：`figures/02_hp_contribution.png`。資料：`data/sm_hp_contribution.json`。

## §1 方法（修正後）
- **chance baseline（修正）**：observed same/diff-HP 只在 **gateable 對（HP1-1 vs HP2-1）**上算（HP3 已排除於投票）。故 baseline 必須**只用 gateable HP**：p(1-1)=10006/19497、p(2-1)=9491/19497 → **chance same-HP = 0.500**。（先前用全類別含 HP3 得 0.443，會虛增 enrichment ~13% — 已修。）
- **ablation**：mutual_excl 在「無 HP 資訊」全當 sibling vs「有 HP gate」只留 same-HP。

## §2 數據觀察（verified；powered+somatic；baseline=0.500）

| 關係 | same-HP | diff-HP | same-HP 率 | vs chance(0.500) |
|---|---|---|---|---|
| independent | 5,198 | 352 | **93.7%** | **1.87×（最高）** |
| nested b⊂a | 5,509 | 518 | 91.4% | 1.83× |
| co_linked | 9,376 | 1,265 | 88.1% | 1.76× |
| nested a⊂b | 5,334 | 924 | 85.2% | 1.70× |
| **mutual_excl** | **3,949** | **5,238** | **43.0%** | **0.86×（DEPLETED 低於背景）** |

## §3 正確判讀（修對抗稽核：same-HP 高 ≠ 克隆證據）

🔑 **同-HP 高是「背景」，不是克隆連鎖的特異證據**：所有「正共現」關係（nested/co_linked/**independent**）same-HP 都 1.7–1.87×，而且 **independent 最高（1.87×）** —— 若 independent 比 nested 更 same-HP，代表「高 same-HP」只是「能被共讀分析的 sSNV 對本來就在同一單倍型上」的**區域背景屬性**，不是 nested 特有的克隆證據。**克隆連鎖的證據是 sSNV 的「共現結構（2×2）」本身，不是 same-HP 率。**

🔑 **唯一偏離背景的是 mutual_excl（DEPLETED, 0.86×）= HP 的真正診斷價值**：互斥**低於**背景 same-HP → 因為互斥同時來自 sibling（同 HP）**與** allelic（異 HP，不同染色體，兩 ALT 永不共讀）。**HP gate 的核心 = 從互斥中移除 allelic（異 HP）**：
- 無 HP 資訊：9,187 對互斥（vs 此 naive baseline）全可能被當 sibling。
- 有 HP gate：只留 3,949 same-HP（真 sibling 候選）。
- **HP 移除 5,238 diff-HP（57%）= 過半是 allelic（不同染色體）非 subclone**。
- ⚠ 此「57%」是對比「**完全無 HP/allele 資訊**」的 baseline；非對比某個更聰明的方法。意義 = HP 提供的「同細胞 vs 不同染色體」資訊不可由共現結構單獨取得。

## §4 限制
- **HP3 = 1,317 linked-somatic（gateable 中 6.3%）無 germline 錨** → 無法 HP-gate，這群不判 sibling/allelic。
- HP=germline-somatic 軸（HP1 vs HP1-1 = 是否帶 phased somatic），非「同單倍型兩拷貝」。同-HP = 同 germline 單倍型，**必要非充分**於同克隆（見 `00b` §3）。
- HP tag 本身有 longphase phasing 誤差 + problem PS block（L3 已標記排除）。

## §5 對論文的修正後意義
HP tag 的 data-supported 貢獻 = **把互斥中的 allelic（57%，異染色體）從 sibling 判定移除**（這是共現結構單獨做不到的）。它**不是**用「same-HP 率高」去確認 nested/co_linked（那是背景）。→ somatic haplotagging 的價值 = **sibling vs allelic 的鑑別器**，精確而有界。
