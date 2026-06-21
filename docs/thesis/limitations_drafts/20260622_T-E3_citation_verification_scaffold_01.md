---
title: "T-E3 citation-verification scaffold（投稿前 web 核對工作清單）"
date: 2026-06-22
status: draft
task: T-E3
audience: 論文作者 + 未來 web-enabled session
build_branch: docs/limitations
note: >
  ⚠ 本 session **沒有 web 工具**（WebSearch/WebFetch/Scholar 不可用）。
  因此本檔**只建立 scaffold + 清單**，明確標出「哪些條目需要 web-enabled session 用
  /citation-verification（WebSearch + Google Scholar）核對」。
  **本 session 未驗證任何 citation，不假裝驗證過。**
data_sources:
  - knowledge/11_external_literature/09_references_and_provenance.md
  - docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2/refs_v2.bib
  - docs/paper_focus/03_references/00_外部文獻索引.md
related_memory:
  - project_external_validation_library
---

# T-E3 — citation-verification scaffold

> **狀態宣告（誠實）**：此 session 為文件軌、**無 web 工具**。本檔產出 = (1) 需驗證項清單；(2) `.bib` 骨架（同目錄 `20260622_thesis_references_skeleton.bib`）；(3) 標明哪些**必須**由 web-enabled session 用 `/citation-verification` 核對。**沒有任何一條在本 session 被「驗證」**。

---

## E3.1 既有 provenance SoT（不重做，直接接續）

論文外部文獻的權威 provenance SoT 已存在：

- **`InterSubMod/knowledge/11_external_literature/09_references_and_provenance.md`** — ~38 篇經 workflow wf_37b2cc97-663（2026-06-05）的 verify agent **實際 WebFetch** 取得真實 PMID/DOI；無法 fetch 者標 `UNVERIFIED`（未編造識別碼）。
- **`InterSubMod/docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2/refs_v2.bib`** — 24 條 deck 用 bib（2026-06-03 /citation-verification 跑過；8 條 CORRECTED）。
- **`/big7_disk/liaoyoyo2001/external_validation/`**（repo 外）— 74 源 CONTEXT 卡（memory `project_external_validation_library`）。

> **關鍵原則（檔 09 自宣告）**：「所有外部文獻**正式投稿引用前，仍須用 /citation-verification 對每篇再核一次**。」→ 本 scaffold 的存在不取代投稿前的逐條 web 複核。

---

## E3.2 🔴 必須 web 核對清單（給 web-enabled session）

> 以下三類**必須**用 WebSearch + Google Scholar（PubMed + 出版商頁面）逐條 4 欄對照（標題 / 第一作者 / 年 / 期刊 + PMID/DOI）才能寫入最終 `.bib`。**本 session 一律不填這些欄位的「verified」狀態**。

### (a) §B 三個跨-agent 識別碼衝突（最高優先 — 檔 09 §B 已標）
| # | 條目 | 衝突內容 | 核對方法 |
|---|------|---------|---------|
| B1 | O'Neill / POG long-read cancer cohort (*Cell Genomics* 2024) | 3 個 article-number（100674 / 100693 / 100538）+ 2 種署名（O'Neill / Lin）| 用 **PMID 39406235** 在 PubMed 核確切 article number + 第一作者；採 anchor PMC11605692 |
| B2 | epihet (*Sci Rep* 2021) | 作者（Chen X / Pan? / Ashoor）+ PMID 不一致 | DOI 10.1038/s41598-020-79627-x / PMC7801679 反查 PubMed 取 canonical PMID + author list |
| B3 | HiCancer (*Sci Rep* 2021) | venue 標註存疑 | DOI 10.1038/s41598-021-86104-6 / PMID 33758310 反查確認 venue |

### (b) UNVERIFIED-PMID 條目（檔 09 已標 `UNVERIFIED`，需補真 PMID）
- Metheor (PLOS Comp Biol 2023, DOI 10.1371/journal.pcbi.1010946)
- ClairS-TO (Nat Commun 2025, DOI 10.1038/s41467-025-64547-z)
- Krishnamachari cfDNA ML FP filter (Sci Rep 2023, DOI 10.1038/s41598-023-37409-1 / PMC10300101)
- Hesson BRCA1 PDAC (Mol Diagn Ther 2022, DOI 10.1007/s40291-022-00614-1 / PMC9626413)
- Herman & Baylin VHL review (J Clin Invest 2001, DOI 10.1172/JCI9462 / PMC289180)
- NanoEM (Nucleic Acids Res 2021, DOI 10.1093/nar/gkab397)
- Soneson batch-class confounding (PLOS One 2014, DOI 10.1371/journal.pone.0100335; PMID 24967636 前 wf 已驗、本輪 agent 標 UNVERIFIED → 重核)
- pycoMeth (Genome Biology 2023, DOI 10.1186/s13059-023-02917-w; PMID 37076875 未直接 fetch)

### (c) 限制段（本批 T-C-UNMASK / T-METHOD-PSLIMIT）新引入或需確認身分的條目
| 條目 | 用在 | 目前狀態 | web 動作 |
|------|------|---------|---------|
| Martin-Trujillo 2017 (PMID 28883545) | T-C-UNMASK | 檔 09 標 high/verified | 確認 PMID + **scope = imprinted DMR**（避免泛化 over-claim）|
| Chase 2015 (PMID 26114957) | T-C-UNMASK | **L3，未在檔 09 清單** | 全新核：標題/作者/年/期刊 + 確認 r≈0.76 與 aUPD 的關係 |
| Rosenski 2025 (PMID 40069157) | T-C-UNMASK | **L3，未在檔 09 清單** | 全新核：確認 460 parental-ASM + 34,426 SD-ASM 座標可作 blocklist |
| MethPhaser / Fu 2024 (PMID 38909018, DOI 10.1038/s41467-024-49588-0) | T-METHOD-PSLIMIT | 檔 09 #2 high；refs_v2.bib `fu2024methphaser` **缺 pmid 欄** | **身分一致性核對**：確認 refs_v2.bib 的 fu2024methphaser（無 pmid）= 檔 09 的 PMID 38909018 是同一篇（Nat Commun 15:5327）→ 合併時補上 pmid |
| Do 2020 (PMID 32594908) | cis 口徑 | 檔 09 #33 high | 確認「+5× hypo-dominant，**口徑差非矛盾**」措辭依據 |
| Tarabichi 2021 (subclonal recon. consensus) | Q2 motivation / reconstruction 用詞 | memory 標「投稿前核 Tarabichi 原句」| 取原文驗 multiplicity 歧義 / 單樣本=characterization 非 confirmation 的原句 |

---

## E3.3 .bib 骨架（同目錄檔）

骨架檔：`20260622_thesis_references_skeleton.bib`

骨架規則（**反假引用，對齊 /citation-verification + §13**）：
1. 每條 `@article` 帶 `note = {STATUS=...}`：
   - `STATUS=VERIFIED_PROVENANCE_SOT`：識別碼來自檔 09（已 WebFetch），但**投稿前仍須 web 複核一次**。
   - `STATUS=NEEDS_WEB_VERIFY`：L3 二手錨 / UNVERIFIED-PMID / §B 衝突 → **web-enabled session 必核才可定稿**。
2. **缺 PMID 的欄位留 `pmid = {NEEDS_WEB_VERIFY}` 佔位，不編造數字**（§13 反捏造）。
3. 骨架只放本批草稿（T-C-UNMASK / T-METHOD-PSLIMIT / T-P2）直接用到的條目 + §B 衝突項；**全論文 bib 應由投稿前一次性 /citation-verification 從檔 09 + external_validation 庫生成**，本骨架不冒充完整 bib。

---

## E3.4 web-enabled session 的下一步（交接指令）

1. 跑 `/citation-verification`（WebSearch + Google Scholar），逐條核 E3.2 (a)(b)(c)。
2. 把 §B 三衝突解決為單一 canonical 識別碼（不並列衝突）。
3. 補齊骨架 .bib 的 `NEEDS_WEB_VERIFY` 佔位，把 `STATUS` 升為 `VERIFIED`。
4. 從檔 09（38 篇）+ external_validation 庫（74 源）生成**全論文** .bib（本骨架只是子集）。
5. 任何 web 查不到的條目 → 標 `NOT_FOUND`，**不寫入正文引用**。

---

## 待填 / 缺資訊

- `{{待填}}` 全論文最終引用清單（須投稿前一次性生成；本檔只列本批用到的子集）。
- `{{待填}}` Chase 2015 / Rosenski 2025 / Tarabichi 2021 的完整 4 欄識別碼（本 session 無 web 不可填）。
