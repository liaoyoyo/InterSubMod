---
globs:
  - "scripts/**/*.sh"
  - "scripts/**/*.py"
---

# Common Workflow Commands

## 1. Full VCF Analysis (default)

```bash
./scripts/run_vcf_all_snv.sh --mode all-with-w5000
```

## 2. Batch Analysis (TP/FP comparison)

```bash
./scripts/run_batch_vcf_analysis.sh
```

## 3. Single-point Quick Verification

```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification
```
