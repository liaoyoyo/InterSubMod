<!--
建立時間: 2026-06-11
狀態: addendum (方法比較研究 → 論文整合 + 開論文就緒確認; 非新主軸 SoT)
報告類型: paper_method_positioning_integration + readiness_check
受眾: 廖子游 · PI · 後續論文 session
framework: Verdict-Pyramid + 就緒 checklist
title: "Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing"
provenance_note: 本檔不重複既有 2 份新主軸 SoT; 只折進本 session 獨有的 method_comparison 研究(00-12) + 確認就緒。framing 對齊既有整合篇章(reconstruction=genetic-haplotagging 驅動, 非甲基)。
-->
<!-- provenance-verified: 既有 SoT 數字引 foundation/整合篇章/convergence; method_comparison 引本 session commit 891e04b…8bb9d4e; 本檔導航非新分析。 -->

# 📎 方法比較研究 → 論文整合 + 開論文就緒確認

> ⚠ **本檔不是主軸 SoT，是 addendum**。新主軸論文的**權威 SoT = 既有 2 份**（2026-06-11 已建）：
> 1. **`InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`**（甲基-phasing-assist 面：V1-V12 + 6 樣本×3 癌種資產 + G-A~G-E gap）— **新主軸 SoT**。
> 2. **`InterSubMod/docs/concepts/2026/06/20260611_Subclonal_Reconstruction_Paper_Focus_整合篇章_01.md`**（ASM-characterization + 四道 NEGATIVE + LOH-phasing 脊柱 + 章節映射 + 框架紀律）。
>
> **本檔唯一職責**：把**本 session 獨有**的「業界方法比較研究」（`docs/method_comparison/.../00-12`）折進論文，並做開論文就緒確認。

---

## L0 — 框架已定（更正我先前的誤判）

> 我先前把「reconstruction vs characterization」當未定決策 —— **錯了，既有整合篇章已精準解決**：
> **「subclonal reconstruction」由 somatic haplotagging（genetic 引擎，Grade B+ ⭐3）驅動；甲基 = characterization（不獨立重建子克隆）。** 標題的 reconstruction **站得住**（是 genetic 產品）；validated KB 04 說「缺 phylogeny/CAMDAC/單細胞」是針對**甲基-based** reconstruction —— **論文不那樣 claim**，所以那 3 支柱是 Discussion roadmap 不是阻擋。
> 🔴 **真正唯一卡關仍是 HD-1**（phasing by-construction 循環 → 跑 R-SELFREF 對照保 positive-spine，或降為 characterization），**你決定**。

---

## L1 — 本 session method_comparison（00-12）→ 論文章節折入

| method_comparison 資產 | → 論文章節 | 用途 |
|----------------------|-----------|------|
| **`02`/`12` 業界 4 軸地景 + 8 共同難點** | **§Related Work / Background** | 定位：業界主流是「軸A tag 比 HP1-vs-HP2 平均 = **germline 軸**」；ISM 站「軸C read-read 距離 + **PERMANOVA 結構檢定** + normal-anchored cis」= 與主流本質不同 |
| **`01` ISM 方法源碼基準 + `03` 對照矩陣** | **§Methods** | ISM 6 核心精確描述 + 與 modkit/DSS/cvlr/DAMEfinder 逐項差異 |
| **`08` 外部程式稽核**（modkit 驗 ISM 解析 **r=0.98**；Fisher 不偏樂觀）| **§Methods validation / reviewer 防守** | 用獨立工具證 ISM 甲基提取正確 + per-CpG 統計不偏樂觀 |
| **`11` 816-loci 大規模驗證**（normal-HP germline 基線 0.237>tumor 0.181 值得加；結構顯著率測弱 8.5%）| **§Methods / Discussion** | 量化「率 vs 結構」互補 + normal-HP 參考價值（對應 §Methods-Neg-4 的 germline-allelic confound） |
| **`04` 結果×文獻交叉** | **§Discussion 防守** | 我們結果 vs Do2020/Martin-Trujillo/POG 的口徑釐清 |
| **`05`/`10` 改進執行計畫 + `12` CAMDAC/MethPhaser** | **§Discussion / future-work + roadmap** | CAMDAC(CN-purity 反卷) + EVOFLUx(phylogeny) = 走向甲基-reconstruction 的 roadmap（補 3 支柱）|
| **修 DAMEfinder 引用**（891e04b）| §References | 引用誠信 |

> 🔑 **最該寫進論文的單一定位句**（本 session 確認）：
> *「業界甲基 haplotype 分析主流是把每 haplotype 塌成 per-CpG 平均、比 HP1 vs HP2（germline 軸）；我們改在 read-level 用 read-read 距離 + PERMANOVA 結構檢定 + normal-anchored somatic cis-test，在 **somatic HP1-vs-HP1-1 軸**量結構而非率差 —— 這是唯一結構性新穎，外部稽核(modkit r=0.98)證明我們的甲基提取與既有工具一致。」*

---

## L2 — 開論文就緒確認（回答你「可以開啟重點論文目標了嗎」）

✅ **可以開** —— 五大支柱齊備：
1. **NEGATIVE 底座（最強）** 🟢 READY — 死四道（整合篇章 §3）。
2. **phasing 脊柱** 🟡 Grade B+ ⭐3 — NG=2 7/7 p=0.0078（HD-1 caveat）。
3. **ASM characterization** 🔵 ⭐3 — chr17 cis / 6/6 excess / BRCA2 退役。
4. **方法定位（本 session）** 🟢 READY — 4 軸地景 + 外部程式稽核 + 816-loci 驗證。
5. **6 樣本×3 癌種資產** — 可衝 ⭐4（foundation §資產）。

🔴 **唯二 gate**（不阻擋開寫，阻擋論文最終強度）：
- **HD-1**：R-SELFREF 循環對照跑 or 降 characterization（你決定）。
- **G-B**：subclone-甲基 somatic-specific vs germline-allelic 的**正確 null（within-haplotype 非 germline-het）**未跑 → 甲基-subclone 貢獻強度押這（整合篇章 §reconcile）。

⏸ **整理放緩（已完成）**：所有現有任務 park；method_comparison Phase B（modkit/DSS anchor 已跑，擴 panel 放緩）；ASM 全基因組 survey 背景。

---

## 下一步（你決策後，後續 session）
1. **HD-1 決策** → 定 phasing 脊柱 grade。
2. 把 method_comparison 定位句 + 4 軸 Related Work 編入 `02_paper_framework/論文架構_正式學術版`。
3. **G-B 正確 null 對照**（within-haplotype）→ 定甲基-subclone 貢獻。
4. （選）CAMDAC-style CN-purity 反卷積 pilot → 補 copy-confound 控制（improvement #13 延伸）。

## Provenance
- 主軸 SoT：foundation_01 + 整合篇章_01（既有，2026-06-11）+ `knowledge/11_external_literature/10`。
- 本 session method_comparison：`docs/method_comparison/20260609_ism_vs_external_methylation_tools/00-12`（commit 891e04b…8bb9d4e）。
- CURRENT_FOCUS pinned 已指向新主軸（2026-06-11）。
- ⚠ 我先前「reconstruction 框架未定」是誤判，已對齊既有整合篇章（reconstruction=genetic 驅動）。
