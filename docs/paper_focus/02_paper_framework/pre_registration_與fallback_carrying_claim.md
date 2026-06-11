<!--
建立時間: 2026-06-09
狀態: pre-registration (carrying claim 雙主軸 + #21/#22 成功判準 + fallback decision tree)
報告類型: paper_focus_preregistration_fallback
受眾: 廖子游 · PI · 執行 #21/#22 的其他 session
framework: scientific-rigor §7 Pre-registration（跑前寫死判準，防 HARKing）
provenance_note: de-confound 主軸數字已 🟢（catalog wf_5644ed77-082 6/6 驗 + 本 session grep）；co-validation 軸 contingent on #21/#22。
-->
<!-- provenance-verified: catalog 數字本 session grep 確認(A=1/C=12868/G=291518/total 332705)；#21/#22 為未跑,判準跑前寫死。 -->

# Pre-registration + Fallback — carrying claim（雙主軸）

> **L0 一眼結論**：carrying claim = **雙主軸＝「規則 + 例外」同一故事**：①**規則（de-confound，已證 🟢）**：read-level 甲基分群真實（12,868 reliable）但**壓倒性是 germline-allelic**，全基因組僅 **1/332,705** 乾淨 somatic cis、甲基不能 filter；②**例外（co-validation，contingent）**：在那唯一位點（chr17）+ phasing 軸上，somatic／haplotype／甲基**正交佐證**同一事件。**de-confound 是 robust floor（不依賴未跑工作）；co-validation 是 contingent ceiling（押在 #21/#22）。**
>
> **L1 重點邏輯**：
> ① **雙主軸不是兩篇** —— 是「**規則（甲基跟 germline 不跟 somatic）+ 例外（罕見的 somatic 真的 co-opt 甲基 = chr17）**」。例外讓規則更有力（rule 因有 well-characterized exception 而完整）。
> ② **floor 已穩**：de-confound 用 catalog（🟢 6/6 驗），#21/#22 任一失敗論文不垮。
> ③ **跑前寫死判準**（§7 防 HARKing）：#21 R-SELFREF + #22-reframed 成功/失敗條件如下。
> ④ **每個 outcome 有 decision rule** → fallback tree。

---

## §A 雙主軸整合敘述（rule + exception，防「散成兩篇」）

> **整合 one-liner**：*在癌症長讀 read-level，**規則**是甲基 allelic 分群跟隨 germline haplotype 而非 somatic 結構（12,868 reliable 分群全 germline-allelic，僅 1/332,705 乾淨 somatic cis，甲基不能 filter 變異）；**例外**是當 somatic 事件罕見地 co-opt allelic 甲基時（chr17/TBC1D16），somatic／haplotype／甲基三訊號正交佐證同一 LOH 重塑事件。*

| 軸 | 角色 | 證據 | 強度 |
|----|------|------|------|
| **de-confound（規則）** | 主 floor | catalog 332,705：C=12,868 germline / A=1 somatic cis / filter 死四道 / dosage-refuted | 🟢 已證、robust |
| **co-validation（例外）** | contingent ceiling | chr17 三訊號對齊（🟡 perm p=0.001）+ phasing NG=2 7/7（contingent on #21）| 🟡 押 #21/#22 |

⚠ **誠實 billing**：兩軸**證據不等強**。de-confound 已證；co-validation 是 1 個 locus（甲基）+ phasing（待 R-SELFREF）。**投稿前若 #21/#22 不如預期 → co-validation 降為 de-confound 的「well-characterized exception」一段，不佔 co-equal headline。**

---

## §B Pre-registration（#21 / #22 — 跑前寫死，§7 防 HARKing）

### PR-1 ｜#21 R-SELFREF circularity 對照（phasing 軸）
- **H（pre-reg）**：NG=2 Inner same-HP1 > Outer 的方向，在**移除/控制 somatic-attributed bucket 的循環成分**後仍成立。
- **度量**：全 7 樣本 flag-on 重跑後，per-sample Inner-vs-Outer gap 的方向 + signed-rank p。
- **🟢 成功判準（pre-reg）**：控制後方向一致 **≥5/7** 樣本 **且** signed-rank p<0.05 **且** median gap 仍 >0（非被 bucket 定義 forced）。→ phasing 當 positive（co-validation 軸成立）。
- **🔴 失敗判準**：控制後方向 <5/7 或 p≥0.05 或 gap 塌到 0。→ **phasing 降為 characterization**，co-validation 軸只剩 chr17。
- **狀態**：未跑（task #21，~25-50hr C++ flag-on，給其他 session）。

### PR-2 ｜#22-reframed（co-validation 例外 + germline 背景）
- **H（pre-reg）**：chr17/TBC1D16 是**唯一**全三訊號（somatic×haplotype×甲基）對齊位點；其餘 reliable 分群為 germline-allelic 背景。
- **🟢 成功判準**：concordance 分析後 TAG-A 仍 =1（chr17 唯一）**且** 12,868 germline-allelic 分群結構可量化（confirm de-confound）。→ 例外+規則都成立。
- **🟡 意外判準**：若出現 chr17 以外的乾淨 somatic co-val 位點（TAG-A>1）→ co-validation 軸升級（更多 exemplar），記為**正向意外**（非 HARKing，因 pre-reg 已說「若 >1 則升級」）。
- **狀態**：confirmatory（catalog 已預示 N=1，低風險），task #22。

---

## §C Fallback decision tree（每個 outcome 論文長怎樣）

```
de-confound floor（catalog 🟢，已證）── 永遠成立，論文不垮
        │
        ├─ #21 成功 + #22 A=1 → 雙主軸：規則(de-confound) + 例外(chr17 co-val + phasing positive)  ★最強
        ├─ #21 失敗 + #22 A=1 → 單主軸：de-confound + chr17 例外；phasing 降 characterization（不當 positive）
        ├─ #21 成功 + #22 A>1 → 雙主軸加強：多個 co-val exemplar + phasing positive  ★★（正向意外）
        └─ #21 失敗 + #22 A=0(chr17 也掉) → 純 de-confound 論文（「甲基全 germline-allelic、零乾淨 somatic cis、不能 filter」）仍可投 GB
```

> **關鍵**：**任一分支論文都站得住**（de-confound floor）。這就是「平行寫 fallback」的意義——主文先寫 de-confound（不等 #21/#22），co-validation 段依 outcome 補/降。

---

## §D 整合 title 候選（雙主軸 = rule+exception）

1. *"Allele-specific methylation in long-read cancer genomes tracks germline, not somatic, haplotypes: a read-level de-confounding with a rare somatic exception"* （rule-led，誠實，推薦）
2. *"What read-level methylation co-validates and cannot filter in cancer: LOH-driven haplotype structure and the rarity of somatic cis-ASM"*
3. *"Read-level co-validation and de-confounding of somatic mutation, haplotype and methylation in cancer long-read sequencing"*（co-val 並列，需 #21/#22 撐）

---

## L1 — 對映 + 下一步
- **floor 可即寫**：de-confound 主文（catalog R6 + filter 死四道 + dosage-refuted）今天就能寫，不等 #21/#22。
- **contingent 待跑**：#21（其他 session，~25-50hr）+ #22-reframed（confirmatory）→ 回來依 §C tree 補 co-validation 段。
- **必修**：same-hap 口徑 3/6 非 6/6（T-PROV）→ phasing 段同步改。
- **對映**：論文架構 R2(phasing)/R3(ASM)/R6(catalog) + 任務 #21/#22。
