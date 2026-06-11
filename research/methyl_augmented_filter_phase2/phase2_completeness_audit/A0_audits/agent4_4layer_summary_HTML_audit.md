# Agent 4 — V6 vs baseline 4-Layer + Summary HTML Audit

**Date**: 2026-05-21
**Auditor**: Agent 4 (read-only)
**Scope**: 2 PI-facing HTML artifacts (4-Layer Synthesis Engineering Report + V6-vs-baseline summary)
**Methodology**: Read both HTML → extract numerical claims → grep source `.md` / `.tsv` / evidence_ledger → label `verified` / `suspect` / `unverifiable`

---

## Files audited

1. `/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/05/20260514_Self_Phasing_4Layer_Synthesis_Engineering_Report_01.html` (100 KB, 1,838 lines)
2. `/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/05/20260514_V6_vs_baseline_HTML_summary_01.html` (25 KB, 376 lines)

Both are companions of validated `.md` reports + are explicitly marked source-verified.

---

## Verdict table

| # | Claim (location) | Value | Source data file | Status |
|---|---|---|---|---|
| 1 | baseline HP1:HP2 family ratio (4Layer §0, §1.2, §1.3; Summary §3.1) | 17.3 : 1 | `20260508_Self_Phasing_完整觀察整合報告_01.md` §2.1 + PI report `20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` | **verified** |
| 2 | V6/V5 全基因組 ratio (4Layer §0, §5.2, §7; Summary §3.3) | 1.84 : 1 | `20260513_V6_Attribution_Errata_01.md` "V6 1.838:1, 64.8%" | **verified** |
| 3 | V3F-only 全基因組 ratio (4Layer §0, §4.4) | 1.85 : 1 | Errata §S5 "V3F-only 全基因組 BAM = 1.85:1" 5/15 empirical | **verified** |
| 4 | chr19 V3F-only best case (4Layer §0, §5.2, §7) | 1.21 : 1 (1.209) | Errata §S5; T1.2 chr19 audit corroborates | **verified** |
| 5 | L1(b) altHaplotype 26,436 : 14,931 = 1.771:1 (4Layer §0, §1.3, §3.2) | 1.77 : 1 | `20260513_V6_Attribution_Errata_01.md` "HAP1 family 26,436 / HAP2 family 14,931" + standalone v1.7-J | **verified** |
| 6 | L2 priority bug amplification (4Layer §0, §1.3) | 9.8× theoretical / 9.4× empirical | Errata "9.4× (17.3/1.85 empirical) / 9.8× (17.3/1.77 theoretical)" | **verified** (both values traceable) |
| 7 | 34,855 read-level victims (4Layer §4.4) | 34,855 | `research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md` "34,855 — PASS" + evidence_ledger cycle `20260508_T1_2_F1_genome_wide_audit` | **verified** |
| 8 | Phase block N50 4,061 → 8,109 (+99.7%) (Summary §3.1) | 4,061 → 8,109 | Errata quote + `20260511_V6_binary_complete_documentation_01.md` | **verified** |
| 9 | Phased rate 54.9% → 78.5% (+23.6pp) (Summary §3.1, §8.6.1) | +23.6pp | Same source as #8 | **verified** |
| 10 | caller F1 = 0.7166 invariant (Summary §3.5; 4Layer §6.1, §8.5) | 0.7166 (TP=28,509, FP=11,606, FN=10,938) | `research/paired_priority_bug_audit/09_V6_caller_F1_verification.md` — 5 stages identical | **verified strongly** |
| 11 | Phase D 4 樣本 ratio: H1437 1.243 / H2009 0.901 / HCC1954 0.958 / HCC1937 0.611 (Summary §3.3) | 0.611 - 1.243 | `research/paired_priority_bug_audit/v6_quantification_findings.md` quoting `phaseD_v6_5sample/v6_cross_sample_summary.tsv` | **verified** |
| 12 | Marker coverage NG≥3: V3F 21,997 → V6 23,980 (+9.0%) (Summary §3.4, 4Layer §8.6.1) | +9% | `20260514_HP_tag_priority_bug_engineering_report_01.md` table; `20260511_V6_binary_complete_documentation_01.md` | **verified** |
| 13 | hp=33 reads 132,060 → 138,317 (+4.7%) (Summary §3.4) | +4.7% | `20260511_V6_binary_complete_documentation_01.md` "V3F 132,060 / V6 138,317 (+6,257 +4.7%)" | **verified** |
| 14 | 15-Clean PS GT +13.3pp (Summary §3.2, §6 bottom line) | +13.3pp (74.9% → 88.2%) | PI report `20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` "74.9% → 88.2% +13.3pp" | **verified** |
| 15 | 17,404 V3F victim subset (Summary §5.2 V6 vs V5 caveat) | 17,404 | Mentioned in summary only; T1.2 全基因組 audit uses 34,855 — 17,404 likely chr19 or specific subset | **suspect** (need source citation; not found in surface grep) |
| 16 | LOH 31.2% self-phasing artifact (4Layer §1.2, §8.1) | 31.2% | `20260508_Self_Phasing_完整觀察整合報告_01.md` §6 | **verified** |
| 17 | Phase block N50 V5/V6 ~20,000+ (4Layer §4.1, §7) | ~20,000 | Summary uses 8,109 (verified); 4Layer says "~20,000+"; numbers do not align internally — 8,109 is genome N50, 20,000 might be sub-region | **inconsistent / suspect** |
| 18 | HP:i:33 ~571 reads (4Layer §4.5, §6.2) | 571 | 4Layer ASCII "HAP3 hp=33: 571" — appears to be derived from L1(b) altHaplotype assignment count, not actual BAM HP:i:33 output. Summary §3.4 says hp=33 = 138,317 reads. The 571 is altHaplotype assignment ≠ BAM tag count | **mislabeled / context risk** (two different "571" vs "138,317" populations use the same hp=33 label) |
| 19 | family-level paired_T alignment 51.95% (V6) / 51.05% (V3F) (Summary §5.1) | ~50% random | `paired_priority_bug_audit/01_step_D_germline_absent_finding.md` (chr19 paired:TO HP1:HP2 ratio) | **verified** |
| 20 | per-read paired_T 對齊 V3F 24.70% / V6 24.97% (Summary §5.1) | ~25% | Same source as #19 | **verified** |
| 21 | 5 commits chain (8b8c1fd / 41ff147 / 380e8d2 / d0bcd8c / 938f0df) (both reports) | 5 commits | evidence_ledger cycle `20260508_self_phasing_synthesis_report` confirms; `20260511_V6_binary_complete_documentation_01.md` documents each | **verified** |
| 22 | Pipeline speed +36% (1× → 1.36×) (Summary §3.1) | +36% | Not found in grep of integration report / errata / V6 binary doc — possibly from internal benchmark not yet published as `.tsv` | **unverifiable** (within audit scope; no source citation in summary) |
| 23 | HCC1395 5kHz purity 0.93 (header) | 0.93 | `paired_priority_bug_audit/09_V6_caller_F1_verification.md` "@ 0.93 purity F1 = 0.7166"; `20260505_self_phasing_V5_data_provenance_audit_01.md` | **verified** |
| 24 | "V5 BAM = V6 BAM 0 diff" / "17,404 V3F victim subset read-by-read 100% identical" (Summary §5.2) | 0 diff | Summary cites internal extract; not located in audit dirs | **unverifiable** in surface grep |
| 25 | 4Layer §0 + §6 5-gate "FILTER write" = 0 matches; diff = 0 lines | 0 / 0 | `paired_priority_bug_audit/09_V6_caller_F1_verification.md` confirms FILTER counts identical 47,798 PASS across versions | **verified** |

---

## Layer-attribution audit (4-layer report specific)

4Layer report's central decomposition `17.3 = 1.77 × 9.8` is internally consistent + matched by both Errata §S5 (9.8× theoretical / 9.4× empirical) and `T1_2_F1_genome_wide_audit.md` (34,855 victims at genome level).

L1(b) out-of-scope framing is defensible — `HaplotagProcess.cpp:161-194` parser GT mapping is in upstream codebase (baseline + V5 share); claim "9 年從未輸出 (HP:i:33 dead code)" rests on enum dual-value collision (HAPLOTYPE3=2 == HAPLOTYPE2=2) — code-level claim cross-checked with `Util.h` claim in §4.5.

Per-layer commit attribution:
- L0 `8b8c1fd` PON-only — verified in evidence_ledger commit hash + V6 binary doc
- L1(a) `d0bcd8c` ploidyRatio — verified, dated 2026-04-30
- L2 `41ff147` V3F two-layer — verified across multiple reports
- L2 minor `380e8d2` INDEL guard — verified
- L1(a) upstream `938f0df` purity threshold — verified (zhenyu, 2025-10-02 cherry-pick)

---

## Cross-document consistency

Both HTML reports are mutually consistent on the headline numbers (17.3:1, 1.84:1, 1.85:1, 1.21:1 chr19, 0.7166, +13.3pp, +9% marker, +99.7% N50, +23.6pp Phased). No contradictory numbers across the two artifacts.

Summary HTML reuses 4-Layer's L1(b) framing in §5.2 caveat (V6 vs V3F trade-off). Internal divergence: Summary §3.4 frames "+9%" as `V6 vs V3F` while 4Layer §8.6.5 calls it "+30.5% abs / +9% rel" (V5 18,382 → V6 23,980) — both numbers correct but anchors differ; reader could mis-attribute baseline.

---

## Findings summary

**Verified (high confidence, source TSV/MD found)**: 20 of 25 numerical claims
**Suspect/inconsistent**: 3 (#15 17,404 subset undocumented; #17 N50 ~20,000+ vs 8,109 inconsistency; #18 hp=33 dual usage 571 vs 138,317)
**Unverifiable in surface grep**: 2 (#22 pipeline speed +36%; #24 V5=V6 BAM 0-diff read-by-read)

**No fabricated numbers detected.** All headline values trace to evidence_ledger cycles `20260508_self_phasing_synthesis_report`, `20260508_T1_2_F1_genome_wide_audit`, `20260507_T1_2_priority_bug_mechanism`, `20260510_V6_proposal_evaluation` + companion `.md` reports. Mermaid/SVG embedded numbers (17.3:1, 1.77:1, 9.8×, 1.84:1, 1.85:1) all reproduce text claims — no hidden discrepancies.

**Recommendations for PI report**:
1. Clarify Summary §3.4 "hp=33 reads 132,060 → 138,317" — this is BAM HP:i:33 tagged-read count, not the 571 altHaplotype HAP3 assignment count in 4Layer §1.3. Two are different populations sharing the "hp=33" label.
2. 4Layer §4.1 / §7 "Phase block N50 ~20,000+" should align with the verified concrete `8,109` value in Summary §3.1 + Errata — drop the "~20,000+" hand-wave.
3. Summary §5.2 V5=V6 "0 diff" claim and "17,404 V3F victim subset" lack a citation in the audit dirs surface grep — add a pointer (likely in `T1_2_read_level_audit/step_B_read_level_dump/` or a new comparison script log).
4. Pipeline speed "+36%" (Summary §3.1) lacks a benchmark `.tsv` — locate or remove.

Overall both HTML artifacts are high-fidelity; PI-facing presentation does not contain fabricated metrics. Main caveats are framing precision (e.g. point #1 hp=33 dual usage) not numerical error.
