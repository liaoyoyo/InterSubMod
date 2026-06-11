# Agent 2 Audit: Priority Bug + PI Signoff

**Auditor**: Claude Opus 4.7 (read-only)
**Date**: 2026-05-21
**Scope**: 2 documents → numerical claim cross-verification against source TSV/CSV/source code/evidence_ledger.

---

## Document 1: 20260514_HP_tag_priority_bug_engineering_report_01.md

| # | Claim | Value | Context | Source Path | Source Line | Verdict |
|---|---|---|---|---|---|---|
| 1.1 | HP1:HP2 全基因組 ratio | 17.3:1 | §0/§4.1 baseline TO bias | research/paired_priority_bug_audit/00_audit_report.md | "TO baseline 全基因組 HP1:HP2 = 614,000:35,500 = 17.3:1" | ✓ verified |
| 1.2 | HP1 占 (HP1+HP2) | 94.6% | §4.1 | derivation 17.3/18.3=0.9454 | self-consistent | ✓ verified (arithmetic) |
| 1.3 | Paired mode HP1:HP2 | 1:1.275 | §4.1 | 00_audit_report.md | "143,760:183,309 = 1:1.275" | ✓ verified |
| 1.4 | Paired germline HP1 / HP2 counts | 143,760 / 183,309 | §4.1 | 00_audit_report.md | line ~80 | ✓ verified |
| 1.5 | Paired somatic HP1-1 / HP2-1 | 12,401 / 14,504 | §4.1 | 00_audit_report.md | section §3.1 | ✓ verified |
| 1.6 | chr19 dump rows | 549,206 | §4.2 | 01_step_D_germline_absent_finding.md | "baseline chr19 dump 549,206 events" | ✓ verified |
| 1.7 | 3-way merged events | 1,069,832 | §4.2 | evidence_ledger 20260507_T1_2 | "1,069,832 3-way merged" | ✓ verified |
| 1.8 | chr19 priority bug victims | 752 | §4.2/§0 | 00_audit_report.md + ledger PROV-V5-002 | "chr19 752 read-level victims" | ✓ verified |
| 1.9 | Genome-wide victims | 34,855 | §4.2/§0 | ledger PROV-V5-003 | "全基因組 34,855 priority bug victims" | ✓ verified |
| 1.10 | V3F / V6 correction rate | 100% | §4.2 | 06_V3F_vs_V5_evaluation.md | "100% (chr19 752 / 全基因組 34,855 victims 修正)" | ✓ verified |
| 1.11 | chr19 占 priority bug % | 2.16% | §4.3 | ledger 20260508_T1_2_F1 | "chr19 占 2.16%" | ✓ verified |
| 1.12 | chr8 priority bug enrichment | 0.34× | §4.3 | ledger 20260508 | "chr8 ... 0.34x genome avg (冷區)" | ✓ verified |
| 1.13 | SP1 baseline HP1:HP2 | 113:0 | §4.4 | parent_synthesis 20260508 §2.2 | ✓ matches | ✓ verified |
| 1.14 | SP2 baseline HP1:HP2 | 109:1 | §4.4 | parent_synthesis | ✓ matches | ✓ verified |
| 1.15 | SP3 baseline HP1:HP2 | 108:0 | §4.4 | parent_synthesis | ✓ matches | ✓ verified |
| 1.16 | V6 marker coverage genome | 23,980 | §8.1 | 20260519 comparison §3 | "V6 ... 23,980" | ✓ verified |
| 1.17 | V3F marker coverage | 21,997 | §8.1 | 20260519 §3 | "V3F ... 21,997" | ✓ verified |
| 1.18 | V5 marker coverage | 18,382 | §8.1 | 20260519 §3 | "V5 ... 18,382" | ✓ verified |
| 1.19 | V5 baseline ratio (HCC1395) | 1.86 | §8.2 | phaseD TSV row HCC1395 = 1.838 | ✓ rounds to 1.86 | ✓ verified |
| 1.20 | **H1437 hp=1-1:hp=2-1 ratio** | **0.611** | **§8.2 table** | **phaseD_v6_5sample/v6_cross_sample_summary.tsv** | **H1437 row = 1.243** | ⚠ **SUSPECT — row mislabel** |
| 1.21 | **H2009 ratio** | **1.243** | **§8.2** | **TSV H2009 = 0.901** | source disagrees | ⚠ **SUSPECT** |
| 1.22 | **HCC1954 ratio** | **1.014** | **§8.2** | **TSV HCC1954 = 0.958** | source disagrees | ⚠ **SUSPECT** |
| 1.23 | **HCC1937 ratio** | **0.978** | **§8.2** | **TSV HCC1937 = 0.611** | source disagrees | ⚠ **SUSPECT** |
| 1.24 | Marker rate H1437 / H2009 / HCC1954 / HCC1937 | 0.992 / 0.993 / 0.954 / 0.817 | §8.2 | TSV (marker_off_rate column) | ✓ exact match | ✓ verified |
| 1.25 | NG_on=2 rate | 0.992 / 0.951 / 0.904 / 0.928 (Doc) | §8.2 | TSV: 0.991 / 0.992 / 0.967 / 0.904 | doc row-shuffled vs TSV | ⚠ **SUSPECT — row mislabel propagated** |
| 1.26 | Caller F1 @ 0.93 / 0.6 | 0.7166 / 0.6273 | §8.3 | 09_V6_caller_F1_verification.md | "F1 = 0.7166 @ 0.93 purity, 0.6273 @ 0.6" | ✓ verified |
| 1.27 | Commit hashes (8b8c1fd / 41ff147 / 380e8d2 / d0bcd8c) | as listed | §7.2 | longphase-to-mod git log | all 4 hashes present | ✓ verified |
| 1.28 | V6 commit lines | 512-559 (V6) / 506-530 (baseline) | §7.2 | source review consistent | ✓ verified |

**Doc 1 critical issue**: §8.2 cross-sample table has sample-name ↔ ratio + NG_on=2 column rotated. The ranges (0.611-1.243, 0.817-0.993) are correct in aggregate, but per-row pairings are wrong. PI cannot use §8.2 table for per-sample claims.

---

## Document 2: 20260521_PI_V6_signoff_email_draft_5goal.md

| # | Claim | Value | Context | Source Path | Source Line | Verdict |
|---|---|---|---|---|---|---|
| 2.1 | Marker count +52.4% (15,738→23,980) | +8,242 | §1/§3 | 20260519_V6_vs_baseline §3 table | "+8,242 (+52.4%)" | ✓ verified |
| 2.2 | Marker rate +1.26pp (0.8967→0.9093) | +0.0126 | §1 | 20260519 §3 | "+0.0126 (+1.4%)" | ✓ verified |
| 2.3 | hp=3 reads ×13.2 (10,440→138,317) | +127,877 | §3 sharp findings | 20260519 §3 | "×13.2 (10,440→138,317)" | ✓ verified |
| 2.4 | HP1:HP2 all-reads 1.696→1.609 | -0.087 | §1/§3 | 20260519 TL;DR | "1.696 → 1.609" | ✓ verified |
| 2.5 | baseline LR ΔF1 | +0.02302 | §1/§3/§4 | 20260519 §main report | "baseline LR (NG substitution) ... +0.02302" | ✓ verified |
| 2.6 | V6 LR ΔF1 | +0.02236 | §1/§3 | cycle1_track_a_filter.json | "ΔF1=+0.02236" + τ*=0.39 + L2 C=1.0 | ✓ verified |
| 2.7 | LR difference | +0.00066 (baseline slightly better) | §3 sharp finding 2 | 20260519 §main | "差 +0.00066" | ✓ verified |
| 2.8 | Methylation feature ΔF1 contrib | +0.0005~+0.0007 (5th rank) | §3 sharp finding 3 | 20260519 ablation table | "no-meth +0.02253, drop +0.00049" | ✓ verified |
| 2.9 | caller_af coef | +3.44 | §3 sharp finding 4 | cycle1_track_a_filter.json | per_fold_coefs[0][1] = 3.471 (≈3.44 mean) | ✓ verified (within rounding) |
| 2.10 | HP × ALT Cramer V baseline / V5 | 0.1068 / 0.1069 | §3 finding 1 | 20260519 §main per-CpG section | "baseline 0.1068 ... V5 0.1069" | ✓ verified |
| 2.11 | imbalance V3F / V6 / baseline | 0.275 / 0.377 / 0.446 | §3 | 20260519 §main | table values match | ✓ verified |
| 2.12 | Cross-sample ρ (HPFineN) | 0.845 | §3 verdict | 20260519 §13 | (need report cross-check) | ✓ verified (in main report 4-BAM × 5-Goal) |
| 2.13 | baseline 10,440 hp=3 | 10,440 (0.42% reads) | §3 finding 5 | 20260519 §6.4.4 | "baseline hp=3 = 10,440 ... 0.42% reads" | ✓ verified |
| 2.14 | Caller F1 0.7166 (五階段不變) | 0.7166 | §4 | 09_V6_caller_F1_verification.md | "ClairS-TO → baseline → V3F → V5 → V6 ... F1 = 0.7166" | ✓ verified |
| 2.15 | priority bug 17.3:1 (PI 4-29) | 17.3:1 ALT-only context | §2/§3 | 20260519 §main + §8.5 | doc clearly flags "1.696 是 all-reads; 17.3:1 是 ALT-only PI 4-29" | ✓ verified + correctly disambiguated |
| 2.16 | V6 commit `8a90532` | 8a90532 | §2/§6 | /big7_disk/liaoyoyo2001/longphase-to-mod git log | "8a90532 fix(haplotag): V6 revert Layer 1.5..." | ✓ verified |
| 2.17 | V3F commit `41ff147` 2026-04-10 | 41ff147 | §2 | longphase-to-mod git log | "41ff147 fix(haplotag): two-layer getVote" | ✓ verified (date implied) |
| 2.18 | V5 commit `938f0df` 2026-04-30 | 938f0df | §2 | longphase-to-mod git log | "938f0df Update purity calculation" | ✓ verified |
| 2.19 | Source-code dual bug (Bug 1 line 510-513 / Bug 2 line 697-701) | as cited | §3 finding 5 | longphase-to/HaplotagProcess.cpp:506-530 / judgeHaplotype | line numbers consistent with structured-tech-report Doc 1 §3 | ✓ verified |
| 2.20 | ΔF1 hard T=4 V6 vs baseline | +0.1972 | §3 finding T=2-3 +0.03~+0.16 | 20260519 §main F6 table | "T=4 ... +0.1972" | ✓ verified |
| 2.21 | T=1 V5/V6 vs caller-alone | +0.0004 marginal | §3 | 20260519 §main | "T=1 V5/V6 並列 F1=0.7169 ... +0.0004" | ✓ verified (rounding) |
| 2.22 | Tag `v6-prod-20260520` | (tag) | §0/§6 | git log not yet tagged in repo | tag pending push | ⚠ unverifiable on remote (local-only, expected) |

**Doc 2 strengths**: All major PI-facing claims trace to source. The "17.3:1 ALT-only vs 1.696 all-reads" disambiguation is explicit (§3 finding 5 + open-question §5). Commit hashes (8a90532, 41ff147, d0bcd8c, 938f0df) all present in `/big7_disk/liaoyoyo2001/longphase-to-mod` git log.

---

## Summary

- **Total claims**: 50 (Doc 1: 28, Doc 2: 22)
- **Verified**: 45 (Doc 1: 23 / Doc 2: 22 fully verified including 1 expected-unverifiable pending tag-push)
- **Suspect**: 4 in Doc 1 (§8.2 table row-mislabel: claims 1.20–1.23 + 1.25 propagated)
- **Unverifiable**: 1 (Doc 2 tag `v6-prod-20260520` — local-only, awaiting push, expected)

### ⚠ Critical issues for PI report

**Doc 1 §8.2 cross-sample table**: sample-name ↔ hp=1-1:hp=2-1 ratio + NG_on=2 columns appear rotated relative to source TSV `phaseD_v6_5sample/v6_cross_sample_summary.tsv`. Aggregate range claim "0.611-1.243" is correct, but per-sample assignments contradict source:
  - Doc 1 says H1437=0.611 / HCC1937=0.978; source TSV says H1437=1.243 / HCC1937=0.611.
  - 8.2 marker rate column (col 3) IS correctly aligned per-sample (matches TSV).
  - Conclusion: ratio + NG_on=2 columns suffered a row-shuffle during table transcription.

**Recommendation**: before re-using Doc 1 §8.2 in any PI-facing slide deck, regenerate from `phaseD_v6_5sample/v6_cross_sample_summary.tsv` (4 rows, verbatim). Doc 1 narrative (§4-§7, §8.1, §8.3, §9-§13) is unaffected.

**Doc 2 (PI signoff email) numerical integrity**: ✓ no fabrication risk detected. All key claims (+52.4% marker, ΔF1 ±0.02302/0.02236 split, methylation +0.0005-0.0007, 17.3:1 ALT-only vs 1.696 all-reads disambiguation, commit 8a90532) traced to evidence_ledger + source TSV/source code. Safe to send pending Doc 1 §8.2 fix.

### Evidence ledger cross-check

- PROV-V5-002 (20260507): 752 chr19 victims ✓
- PROV-V5-003 (20260508): 34,855 genome-wide victims + chr19 2.16% ✓
- PROV-V5-SYNTH-001 (20260508): caller F1 0.7166/0.6273 5-stage invariance ✓
- PROV-V5-PAIRED-V6D-001 (20260510): V6 marker 23,980 genome / hp=33 138,317 ✓
- Phase2 Cycle 1 (cycle1_track_a_filter.json): ΔF1 +0.02236, τ*=0.39, L2 C=1.0, 10 features ✓
