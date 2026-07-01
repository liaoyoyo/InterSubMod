<!--
建立時間: 2026-07-01
類型: L10 甲基「兩層」大規模驗證 + read 招募力（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3）
data_sources: data/sm_two_layer_validation.json, data/sm_two_layer_validation_perregion.tsv, data/sm_read_recruitment.json
build_branch: feat/summary-nreadsvalid
-->

# L10 — 甲基「兩層」大規模驗證 + read 招募力

> **緣起（使用者）**：甲基資訊分「haplotype 層（germline ASM）」與「subclone 層（somatic-specific）」的假設，如何**大規模**驗證（非只 49 corroborated 小 N）？如何判斷資訊是否有用？甲基錨點能否招募更多 read（含未穿 sSNV / unphased / HP3）歸群，串接更遠的 read、補上無法建樹的位點？
> 本層用**兩個新全基因組分析**（754 區）逐一數據驗證。HP 比對用 `str()`（延續 2026-06-29 str/int bug 修正）。

## §0 一句話總答

🔴 **兩層假設在全基因組成立且非對稱**：**haplotype 層（germline ASM）真實、強、普遍**（neutral 區 62%、Δβ 0.97、非循環）→ 對 **phasing/read 招募**有用；**subclone 層（somatic-specific）稀少**（全基因組 7–9 區）→ 需 matched-normal 才非循環。**甲基錨點能以高準確度把 read 招募到 HP（haplotype），但招募不到 subclone。**

## §1 大規模驗證設計（三軸，754 區單次 fetch）

| 軸 | 測什麼 | 循環性 |
|---|---|---|
| **(A) HP-axis germline ASM** | HP1{1,1-1} vs HP2{2,2-1}，所有 CpG，MWU+BH-FDR | ✅ **完全非循環**（純 germline 軸，不碰 somatic）|
| **(B) genotype-axis somatic** | 兩大基因型 population（= 重現 corroboration）| — |
| **(C) 兩軸交叉分解** | genotype-sig CpG 中，HP 也 sig（germline-cis）vs 殘差（subclone 候選）| — |

## §2 [A] Haplotype 層 — ✅ 真實、強、普遍（大規模非循環證實）

| 指標 | 值 |
|---|---|
| HP 可評估區 | **175 / 754** |
| germline ASM 區 | **51（29.1% of 可評估）** |
| **neutral（germline ASM 本可能）** | **39 / 63 = 62%** ⭐ |
| LOH（單倍型，本不該有）| 12 / 112 = 11% |
| effect size | max Δβ median **0.972**、**100% ASM 區 Δβ>0.5**、median **15** sig CpG/區 |

→ **在 germline ASM 結構上可能的 neutral 區，62% 真有，且效應極強（Δβ≈0.97 all-or-nothing）。** 非循環鐵證：haplotype 層是普遍現象非特例。**［L1 全基因組重算，`sm_two_layer_validation.json`］**

## §3 [B] Somatic 軸 — corroboration 完美重現

- **49 corroborated**（重現 ✓）；其中 17 HP 可評估、**8 germline-explained（47%，此判準）**。
- ⚠ 與 L8「11/14」差異＝判準不同（此處看 geno-sig∩HP-sig CpG 重疊比例；L8 看 region 級 HP1-vs-HP2 MWU）；方向一致（約半至多數為 germline-cis），精確比例隨判準 47–79%。

## §4 [C] 兩層交叉表 — 決定性圖像

| 類別 | 區數 | 意義 |
|---|---|---|
| **純 germline ASM（有 ASM、無 somatic 對齊）** | **41** | haplotype 層獨立存在，不靠 somatic 驅動 |
| germline ASM + corroborated | 10 | germline-cis 重疊（「對齊」其實是等位特異）|
| **corroborated 但無 germline ASM** | **7** | somatic 對齊但非等位特異 → **最可能真 somatic/subclone** |
| subclone 候選（殘差判準）| **9**（LOH 6 + neutral 3）| 隨判準 3–9 |

## §5 逐步確認：7 個「真 somatic 候選」

- **6 個 LOH**（chr5/9/1/10/14/19，`linear_nested`）：HP 讀數極度偏斜（119/6、115/4、5/124…）→ 符合 LOH 單倍型；其 somatic 甲基差異最可能 **somatic-cis（突變直接效應）**，需 normal 排除。
- **1 個 neutral：`chr21:30520863`（co_linked，HP 11/26 平衡，5 sig CpG）** ← **最乾淨單一候選**：CN neutral（germline ASM 本可能）、HP 平衡、corroborated、但 HP1≠HP2 沒解釋 → 甲基差異是**基因型連鎖但非等位特異** = 真 somatic-specific 甲基化最強單例。**［L1 `sm_two_layer_validation_perregion.tsv`］**

## §6 判斷「資訊是否有用」的框架

**判準 = 普及率 × effect size × 非循環性 × 判別力**：

| 層 | 普及率 | effect | 非循環 | 判別力 | 有用？ |
|---|---|---|---|---|---|
| **Haplotype（germline ASM）** | **62%**(neutral) | Δβ 0.97 | ✅ | 強 | ✅ **有用 → phasing/招募先驗** |
| **Subclone（somatic-specific）** | **7–9 區/754** | 高但少 | 需 normal | 未證 | 🟡 **有限 → 需 matched-normal** |

## §7 read 招募力（甲基錨點能招募多少 read）— 全基因組（`sm_read_recruitment.json`）

51 個有 germline-ASM 錨點的區域，測甲基錨點的招募力：

| 招募對象 | 方法 | 結果 | 判定 |
|---|---|---|---|
| **① 招募到 HP（haplotype 層）** | 錨點 CpG 甲基 leave-one-out 最近質心分類 HP-tagged read | **4,956 / 5,147 = 96.3%**（42 區，per-region median **98.6%**）| ✅ **高準確度可行** |
| **② 可招募 unphased/HP3 池** | hp0 read 覆蓋 ≥3 錨點 CpG | **334 read**（占 anchored 區 hp0 總 4,498 的 **7.4%**；對比 HP-tagged 12,598）| 🟡 **真但小**（覆蓋限制）|
| **③ 招募到 subclone（within-HP 子群）** | 同一 HP 內 ≥2 somatic 子群，甲基能否分開 | **2 / 19 = 10.5%** | ❌ **多數失敗（data-starved）** |

→ **決定性**：甲基錨點**能以 96% 準確度把 read 招募到 HP（haplotype）**，可救回 **334** 個 unphased/HP3 read 給 phasing；但**招募到 subclone（within-HP 子群）只 10.5%** → **甲基串接 read 是 phasing 層的擴充，不是 subclone 層的擴充**。**［L1 全基因組 LOO，`sm_read_recruitment_perregion.tsv`］**

> **對使用者「串接更遠 read → 更準演化樹」的直接回答**：
> - ✅ **能救 unphased/HP3 → HP**（334 read、96% 準確）→ 擴充**已 phase 的骨幹覆蓋**，間接讓更多 read 可用。
> - ❌ **不能把「只被 germline-phase 標到 HP1、未穿 sSNV」的 read 可靠地判成某 subclone**（within-HP 甲基分群 10.5%）→ 這條鏈**斷在「HP 內再分 subclone」**。
> - 淨效果：甲基招募**增加 phasing 覆蓋與 read 利用率**，但**不新建 subclone 樹**；補「無法建樹位點」靠的是遺傳（更多 sSNV / 更深 read）非甲基。

## §8 結論（給研究決策）

1. **甲基主要訊號 = germline 等位特異甲基化（haplotype 層）**：普遍（62% neutral）、強（Δβ 0.97）、非循環 → 對 **phasing/read 招募**有用（甲基真正的研究價值），**非 subclone 偵測**。
2. **真 somatic-specific 甲基全基因組僅 7–9 區**（6/7 LOH 最可能 somatic-cis，唯 chr21 neutral 乾淨候選）→ 甲基當 subclone 判別器上限就這幾區，都需 matched-normal 收尾。
3. **紅線**：甲基錨點招募到的是 **haplotype 非 subclone**；把「招募到 HP」誤稱「招募到 subclone」= baseline-dependence 循環。
