# cis-ASM characterization pipeline — Makefile (run from pipeline/)
# characterization only; no TP/FP discriminator. Python-first; no C++ this phase.
PY := python3

.PHONY: help test features clean
help:
	@echo "make test      — BRCA2 canonical regression test (golden guard; run before any commit)"
	@echo "make features  — recompute the 10-locus feature matrix (needs BAM + annotations)"

test:                  ## golden guard — must PASS before committing any stage change (runs all tests; BAM/bed tests SKIP if data absent)
	@for t in tests/test_*.py; do echo "--- $$t ---"; $(PY) "$$t" || exit 1; done

features:              ## full feature matrix (Level-1 + BAM + annotations)
	cd .. && $(PY) scripts/35_feature_matrix_trial.py

clean:
	rm -rf lib/__pycache__ tests/__pycache__
