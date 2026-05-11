# Hypothesis H[XXX]: [名稱]

## Status: pending | testing | concluded
## Verdict: — | positive | negative | NO-GO | conditional
## Track: paired_full | paired_pileup | TO
## Priority: [1-100]

---

## Hypothesis Statement

[一句話假說陳述，包含可否證條件]

**Falsifiable condition**: [具體的推翻條件，含數值閾值]

---

## Verification Ladder

| Level | Date | Result | Key Metric | Evidence Path | Pass? |
|-------|------|--------|------------|---------------|-------|
| L1: AUC Screening | | AUC=X.XX | AUC > 0.58 | | |
| L2: Confound Check | | | residualized AUC > 0.58 | | |
| L3: Cross-sample | | N/7 consistent | >=5/7 same direction | | |
| L4: Beyond-AUC | | | mechanism + literature | | |

### L1 Details
- Dataset: [HCC1395_5kHz / ...]
- Baseline F1: [X.XXXX]
- Result F1: [X.XXXX]
- Delta: [+/-X.XXXX]
- TP/FP/FN changes: [+N/-N/+N]

### L2 Details
- Confounders tested: [AF, n_reads, LOH, ...]
- Method: [within-group OLS residualization / AF-bin stratification / ...]
- Raw AUC → Residualized AUC: [X.XX → X.XX]

### L3 Details

| Sample | AUC | Direction | Notes |
|--------|-----|-----------|-------|
| HCC1395 | | | |
| HCC1395_DORADO | | | |
| HCC1937 | | | |
| HCC1954 | | | |
| H2009 | | | |
| H1437 | | | |
| COLO829 | | | |

### L4 Details
- Mechanism explanation: [biological rationale]
- Literature support: [papers/citations]
- Alternative explanations ruled out: [list]

---

## Conclusion

**Verdict**: [positive / negative / NO-GO / conditional]

**Stability**: [high / medium / low] — [reopen conditions]

**Impact on research direction**: [how this changes the strategy]

**Excluded alternatives**: [what was ruled out and why]

---

## Evidence Ledger Entries

| Cycle ID | Date | Scale | Delta F1 | Verdict |
|----------|------|-------|----------|---------|
| | | | | |

---

## Related Hypotheses

- [H_YYY]: [relationship — supersedes / contradicts / extends]
