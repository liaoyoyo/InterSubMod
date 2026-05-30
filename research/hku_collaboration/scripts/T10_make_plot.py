#!/usr/bin/env python3
"""T10 plot finalizer — reads complete 24-chr TSV and emits PNG.

Used after partial-run resume: reads the full TSV produced by main script
(possibly assembled from multiple --append-tsv runs) and produces the
dual-axis bar chart per spec.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Re-use the make_plot function from the main script
sys.path.insert(0, str(Path(__file__).parent))
from T10_hp3_tp_rate_24chr import make_plot, CHROMS  # noqa: E402


def main() -> int:
    tsv_path = Path(
        "/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/findings_5_24/T10_HP3_TP_rate_24chr.tsv"
    )
    png_path = Path(
        "/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration/figures/T10_HP3_TP_rate_24chr.png"
    )

    rows_by_chr: dict[str, dict] = {}
    with tsv_path.open() as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            d = dict(zip(header, parts))
            d["total_reads"] = int(d["total_reads"])
            d["hp3_reads"] = int(d["hp3_reads"])
            d["hp3_fraction"] = float(d["hp3_fraction"])
            d["hp3_tp_reads"] = int(d["hp3_tp_reads"])
            d["hp3_tp_rate"] = float(d["hp3_tp_rate"])
            d["seqc2_truth_count"] = int(d["seqc2_truth_count"])
            rows_by_chr[d["chromosome"]] = d

    # Order by canonical CHROMS list
    rows = [rows_by_chr[c] for c in CHROMS if c in rows_by_chr]
    missing = [c for c in CHROMS if c not in rows_by_chr]
    if missing:
        print(f"[T10-plot] WARN: missing rows for: {missing}", file=sys.stderr)
    print(f"[T10-plot] {len(rows)}/{len(CHROMS)} chrs in TSV", file=sys.stderr)

    make_plot(rows, png_path)
    print(f"[T10-plot] PNG written: {png_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
