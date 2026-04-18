---
globs:
  - "output/**/*"
  - "results/**/*"
---

# Output File Structure

Each Region directory contains:

- `metadata.txt`: Region basic info and statistics
- `reads.tsv`: Read detailed information
- `methylation.csv`: Methylation matrix (Read × CpG)
- `distance_matrix_[METRIC].csv`: Read-Read distance matrix
- `significance_summary.csv`: Significance analysis summary
- `*.png`: Visualization charts
