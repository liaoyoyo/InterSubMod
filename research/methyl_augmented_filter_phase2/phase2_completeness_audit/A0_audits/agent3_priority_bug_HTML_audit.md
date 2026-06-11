# Agent 3 Audit: Priority Bug Standalone HTML (2 files)

**Auditor**: Claude Opus 4.7 (read-only)
**Date**: 2026-05-21
**Scope**: 2 standalone HTML reports → numerical claim cross-verification against source MD / TSV / evidence_ledger / C++ source line refs.

**Files audited**:
1. `InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_17to1_rootcause_explained_01.standalone.html` (898 lines, 41 KB)
2. `InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_priority_bug_engineering_report_01.standalone.html` (1373 lines, 94 KB)

**Note**: file (2) header banner declares itself **SUPERSEDED** by `20260514_Self_Phasing_4Layer_Synthesis_Engineering_Report_01.html`; retained as amendment-history audit trail. Numbers below still audited as the file lives in `validated/`.

---

## Document 1: 20260514_HP_tag_17to1_rootcause_explained_01.standalone.html

Source path used for each claim: parent synthesis = `docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`; paired audit = `research/paired_priority_bug_audit/{00,01,06,07,08,09}_*.md`; evidence_ledger entries `PROV-V5-002 / PROV-V5-003 / PROV-V5-SYNTH-001`.

| # | Claim | Value | Section | Source Line | Verdict |
|---|---|---|---|---|---|
| 1.1 | HP1 系列 reads count | ~614 K | §1 stat-grid | 00_audit_report.md:102 "614,000:35,500" | ✓ verified |
| 1.2 | HP2 系列 reads count | ~35 K | §1 stat-grid | 00_audit_report.md:102 | ✓ verified |
| 1.3 | HP1:HP2 全基因組 ratio | 17.3 : 1 | §1 stat / TL;DR | 20260508 §2.1 line 101; 00_audit_report.md:32 | ✓ verified |
| 1.4 | HP1 占 (HP1+HP2) % | 94.6% | §1 / TL;DR / §5 arithmetic | 20260508 §2.1 line 102; arithmetic 17.3/18.3=0.9454 | ✓ verified |
| 1.5 | chr19 dump rows (each version) | 549,206 | §1 / §8 details | ledger PROV-V5-002 "549,206 vote dumps × 3 binary versions" | ✓ verified |
| 1.6 | 3-way merged events | 1,069,832 | §1 table | ledger PROV-V5-002 "1,069,832 3-way merged" | ✓ verified |
| 1.7 | chr19 priority bug victims | 752 | §1 / §8 / TL;DR | ledger PROV-V5-002 + 06_V3F_vs_V5:32; 20260508 §4.2 line 267 | ✓ verified |
| 1.8 | Genome-wide victims | 34,855 | §1 / §8 | ledger PROV-V5-003 "全基因組 34,855 priority bug victims"; 06_V3F_vs_V5:32 | ✓ verified |
| 1.9 | V6 修正率 | 100.00% | §1 / §8 | 06_V3F_vs_V5:32 "100% ... 修正"; 20260508 amend banner | ✓ verified |
| 1.10 | 反向案例 (baseline=21 → V6=11) | 0 | §1 / §8 | 20260508 §4.2 line 267 "全 752 條 ... 無一條反向" | ✓ verified |
| 1.11 | VCF GT2 偏 `1\|.` 比例 | ~95% | TL;DR / §4 | 4Layer-synthesis §1 "1.77:1 偏 HP1" (precise = 26,436:14,931 → 63.9%); HTML wording "約 95%" derives from old framing | ⚠ **IMPRECISE — pre-erratum framing** (correct value = 1.77:1 = 63.9%; 95% conflates GT2 偏向 with downstream HP1 占比 94.6%) |
| 1.12 | chr19 占 priority bug % | 2.16% | §5 | ledger 20260508_T1_2_F1 "chr19 占 2.16%"; 20260508 line 75 | ✓ verified |
| 1.13 | chr8 enrichment | 0.34× genome avg | §5 | 20260508 line 314 "chr8 enrichment 0.34× genome avg" | ✓ verified |
| 1.14 | chr8 LOH+HPSig FP enrichment | 7.4× | §5 / cross-ref | MEMORY (`project_hcc1395_chr8_hotspot.md`) "LOH+HPSig 7.4× FP enrichment"; flagged "不同 layer" — correct | ✓ verified |
| 1.15 | SP1 chr19:17,565,944 HP1:HP2 | 113 : 0 | §5 | 20260508 §2.2 line 120; 00_audit_report.md:147 | ✓ verified |
| 1.16 | TO mode HP1:HP2 vs Paired | 17.3:1 vs 1:1.275 | §6 table | 00_audit_report.md:32 + line 97 | ✓ verified |
| 1.17 | Cross-sample ratio range (4 samples) | [0.611, 1.243] | §7 / §8 | 08_phaseD §117 "1.243 / 0.901 / 0.958 / 0.611" → range OK | ✓ verified |
| 1.18 | V5 germline-absent 4.19:1 | 4.19:1 | §7 trade-off | 01_step_D_germline_absent:75-77,93 "3,313:790 = 4.19:1" | ✓ verified |
| 1.19 | Caller F1 (HCC1395 0.93 / 0.6) | 0.7166 / 0.6273 | §7 / §8 gate 5 | 09_V6_caller_F1_verification:24,70-77 (TP/FP/FN identity-preserved) | ✓ verified |
| 1.20 | Source lines: baseline 506-530 / V6 512-559 | as cited | §2/§7/footer | 20260514 .md frontmatter inputs §; line numbers consistent across all cross-refs | ✓ verified |
| 1.21 | Util.h Haplotype enum lines | 19-26 | §2 table | 20260514 .md frontmatter "Util.h (lines 19-26 Haplotype enum)" | ✓ verified |
| 1.22 | Arithmetic 17.3/18.3 = 0.9454 ≈ 94.6% | identity | §5 | self-consistent + 20260508 §4.4 line 331 | ✓ verified |

**Doc 1 issue (1 SUSPECT)**: Claim 1.11 "VCF GT2 約 95% 偏 1|." in TL;DR §4 is **pre-erratum framing**. The 4-Layer Synthesis report (which supersedes this) decomposes the bias as L1 assignment 1.77:1 (63.9%) × L2 priority bug 9.8× → 17.3:1 (94.6%). The HTML conflates `94.6% HP1 系列 BAM 統計` (downstream outcome) with `~95% GT2 偏向` (upstream cause). Source data (4Layer §1, ASCII flow) shows GT2 偏 HP1 is only `26,436:14,931 = 63.9%`, NOT 95%. This is a holdover from the original (pre-erratum) causal framing.

---

## Document 2: 20260514_HP_tag_priority_bug_engineering_report_01.standalone.html

This doc has SUPERSEDED banner but is still inside `validated/`. The .md companion was already audited by Agent 2; this audit focuses on HTML-specific claims (chart inline data, stat-box numbers, embedded SVG labels).

| # | Claim | Value | Section | Source Line | Verdict |
|---|---|---|---|---|---|
| 2.1 | TL;DR stat-box HP1:HP2 | 17.3:1 | §0 | 00_audit_report.md:32,102 | ✓ verified |
| 2.2 | TL;DR victims | 34,855 | §0 | ledger PROV-V5-003 | ✓ verified |
| 2.3 | TL;DR V6 修正率 | 100% | §0 | 06_V3F_vs_V5:32; ledger PROV-V5-002,003 | ✓ verified |
| 2.4 | 5 / 5 verification gates | 5 / 5 | §0 | 07_V6_validation_findings + 20260508 §6.2 | ✓ verified |
| 2.5 | §0 erratum banner: 17.3 = 1.77 × 9.8 decomposition | 1.77:1 × 9.8× | §0 erratum / §5.4 | 4Layer-synthesis §1 line 278 "1.77:1 × 9.8×"; line 488 "26,436/14,931 = 1.771"; line 500 "17.3/1.77 = 9.8" | ✓ verified |
| 2.6 | §4.1 HP1 / HP2 series counts | ~614 K / ~35 K | §4.1 table | 00_audit_report.md:102 | ✓ verified |
| 2.7 | §4.1 paired counts (HP:Z:1 / 2 / 1-1 / 2-1) | 143,760 / 183,309 / 12,401 / 14,504 | §4.1 table | 00_audit_report.md:90-93 (exact match) | ✓ verified |
| 2.8 | §4.1 paired HP1:HP2 series ratio | 1 : 1.275 | §4.1 table | 00_audit_report.md:99 "1:1.267"(family) / line 32 "1:1.275"(整體) | ✓ verified |
| 2.9 | §4.2 chr19 dump rows / merged events | 549,206 / 1,069,832 | §4.2 | ledger PROV-V5-002 | ✓ verified |
| 2.10 | §4.2 chr19 victims / genome victims | 752 / 34,855 | §4.2 | ledger PROV-V5-002,003 | ✓ verified |
| 2.11 | §4.2 V3F 修正率 / V6 修正率 / 反向 0 | 100% (752/752) / 100% (34,855/34,855) / 0 | §4.2 | 20260514 .md:165-172; 06_V3F_vs_V5:32 | ✓ verified |
| 2.12 | §4.3 SP1/SP2/SP3 baseline | 113:0 / 109:1 / 108:0 | §4.3 table | 20260508 §2.2 line 120 + §6 line 585; 00_audit_report.md:147 | ✓ verified |
| 2.13 | §4.3 SP1 paired som_ratio | 0.500 (265:265) | §4.3 | 00_audit_report.md:126,147 | ✓ verified |
| 2.14 | §5.4 L1 assignment HAP1+HAP1_1 | 26,436 | §5.4 code block | 4Layer-synthesis line 488 "26,436 / 14,931 = 1.771" | ✓ verified |
| 2.15 | §5.4 L1 assignment HAP2+HAP2_1 | 14,931 | §5.4 | 4Layer-synthesis line 488 | ✓ verified |
| 2.16 | §5.4 L2 amplification | 9.8× (17.3/1.77) | §5.4 | 4Layer-synthesis line 500 | ✓ verified |
| 2.17 | §5.4 PON-only V5 assignment | 2.03:1 | §5.4 / §5.3 table | 4Layer-synthesis line 389 "baseline 1.77:1 / V5 2.03:1" | ✓ verified |
| 2.18 | §5.4 v3f_no_pononly chr19 result | 1.21:1 | §5.4 / 4Layer §3.1 | 4Layer-synthesis line 673 "chr19 best case 1.21:1" | ✓ verified |
| 2.19 | §8.1 marker coverage V3F / V5 / V6 | 21,997 / 18,382 / 23,980 | §8.1 table | 07_V6_validation:179-181; 20260508:1214-1216 | ✓ verified |
| 2.20 | §8.1 V5 germline-absent 4.19:1 | 4.19:1 | §8.1 row | 01_step_D_germline_absent:75-77 | ✓ verified |
| 2.21 | §8.1 V5 全基因組 ~1.86:1 | ~1.86:1 | §8.1 row | 08_phaseD line 14,96,131 "V5 1.86" (rounded from 1.838 in 07_V6_validation:111,162 "V5 1.863") | ⚠ **rounding-only discrepancy** — 1.863 → "1.86" OK |
| 2.22 | §7.2 commit chain hashes (4) | 8b8c1fd / 41ff147 / 380e8d2 / d0bcd8c | §7.2 + SVG fig 4 | 20260508 line 16 commit list; 07_V6_validation:35,231,247,306 | ✓ verified |
| 2.23 | §7.2 V3F TP region count | n/a (not claimed) | — | — | — |

**Doc 2 issues**:
- **No fabrication found.** The 4-Layer decomposition (1.77 × 9.8) is rigorously source-traced to `20260514_Self_Phasing_4Layer_Synthesis_Engineering_Report_01.html`. The HTML is **internally consistent** with that synthesis.
- **Minor rounding (claim 2.21)**: "V5 baseline ~1.86" rounds 1.863 (07_V6_validation:111) / 1.838 (08_phaseD:117 for HCC1395). Acceptable (<2% delta).
- **SUPERSEDED banner** on top is honest disclosure; HTML retained for audit trail.

---

## Cross-document consistency

| Cross-claim | Doc 1 | Doc 2 | Same source | Verdict |
|---|---|---|---|---|
| HP1:HP2 ratio | 17.3:1 | 17.3:1 | 20260508 §2.1 | ✓ consistent |
| 752 / 34,855 victims | 752 / 34,855 | 752 / 34,855 | ledger PROV-V5-002,003 | ✓ consistent |
| 549,206 dump rows | 549,206 | 549,206 | ledger PROV-V5-002 | ✓ consistent |
| 113:0 SP1 baseline | 113:0 | 113:0 | 20260508 §2.2 | ✓ consistent |
| V6 cross-sample range | [0.611, 1.243] | [0.611, 1.243] | 08_phaseD line 117 | ✓ consistent |
| Caller F1 0.7166 / 0.6273 | both cite | both cite | 09_V6_caller_F1 | ✓ consistent |

No inter-document contradictions.

---

## Overall Verdict

**Doc 1 (17to1_rootcause)**: 21 / 22 claims verified. **1 SUSPECT** = §4 TL;DR "VCF GT2 約 95% 偏 1|." conflates upstream L1 assignment ratio (correct 1.77:1 = 63.9%) with downstream BAM HP1 占比 (94.6%). Pre-erratum framing leak; technically wrong but the headline number 17.3:1 / 94.6% / 752 / 34,855 are all correct. Recommend PI footnote or §4 sentence patch.

**Doc 2 (engineering_report)**: 22 / 22 claims verified (1 minor rounding 1.86 ≈ 1.863 acceptable). HTML is internally consistent with both its source .md and the 4-Layer Synthesis supersedure. SUPERSEDED banner is honest disclosure.

**No fabrication detected** (i.e., no integer-rounded "LLM smoothed" numbers like 1000 / 50% etc. — all precise figures trace to TSV / .md / ledger). Both files are PI-safe to cite, with caveat on Doc 1 claim 1.11 phrasing.

**Risk level**: LOW. No 2026-05-20 +0.057-fabrication-class issues present.

**Read-only confirmation**: no files modified; only `agent3_priority_bug_HTML_audit.md` created (matches request).
