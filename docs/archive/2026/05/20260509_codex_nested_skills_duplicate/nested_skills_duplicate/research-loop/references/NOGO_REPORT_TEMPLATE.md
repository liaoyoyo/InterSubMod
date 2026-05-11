# NO-GO Report: [方向名稱]

**Date**: YYYY-MM-DD
**Verdict**: NO-GO
**Confidence**: [high / medium]

---

## One-Line Conclusion

[一句話總結為何此方向不可行]

---

## Verification Summary

| Level | Result | Key Number |
|-------|--------|------------|
| L1 AUC | AUC = X.XX | threshold: 0.58 |
| L2 Confound | [passed/failed] | [residualized AUC] |
| L3 Cross-sample | N/7 consistent | |
| L4 Beyond-AUC | [if reached] | |

---

## Exclusion Reason Category

- [ ] **Statistical Ineffective**: AUC < 0.58 after confound correction
- [ ] **Confound Dominated**: signal entirely from known confounders (AF, n_reads, LOH, HP tags)
- [ ] **Sample Specific**: works on <=2/7 samples only
- [ ] **Mechanism Implausible**: no biological rationale or contradicts established knowledge
- [ ] **Practical Constraint**: works in theory but safety constraints prevent application (e.g., TP loss > 2%)

---

## Key Evidence

1. [Most important piece of evidence — include specific numbers]
2. [Second most important]
3. [Third — if needed]

**Evidence paths**: [links to cycle artifacts, figures, reports]

---

## Prevent Re-investigation

### Variants Already Tested

| Variant | Result | Why Different from Base |
|---------|--------|----------------------|
| [variant 1] | [AUC/result] | [what was changed] |
| [variant 2] | [AUC/result] | [what was changed] |

### Exception Conditions (reopen if...)

1. [Specific condition that would invalidate the NO-GO, e.g., "new low-purity samples become available"]
2. [Another condition, if any]

**If none of these conditions are met, do not revisit this direction.**

---

## Impact

- **Closed hypothesis IDs**: [H_XXX, H_YYY, ...]
- **Related memory**: [memory file name]
- **Research direction adjustment**: [how this changes the strategy going forward]
