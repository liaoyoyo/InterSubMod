"""
genomic_context.py — P2 LOH + CNV integration (the biggest gap from the analysis).

Treats LOH/CN as FIRST-CLASS observation + an axis-validity gate, not a text label:
  - cn_annotate: integer total CN from SEQC2 gain_cn/loss_cn beds (not just gain/loss labels).
  - loh_status: point-in-region lookup against the LOH bed.
  - axis_validity: HP-axis (HP1 vs HP1-1) is structurally valid ONLY in non-LOH (both
    haplotypes present). In LOH one haplotype is lost -> HP-axis is reconstructing against a
    copy-lost haplotype, so primary axis switches to ALLELE-axis (valid at VAF<1).
    This fixes the 43% (10,220/23,840) HP-axis-in-LOH violation found in the audit.

Guardrail: characterization annotation only; produces NO TP/FP discriminator score.
Data sources are SEQC2 matched-normal NGS benchmarks (whole-region, coarser than CpG-window).
"""
import os

# SEQC2 CNV / LOH (matched-normal NGS benchmark)
SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_CN = os.path.join(SEQC2_DIR, "ngs_benchmark_cnv_gain_cn.bed")    # col4 = integer CN (>2)
LOSS_CN = os.path.join(SEQC2_DIR, "ngs_benchmark_cnv_loss_cn.bed")    # col4 = integer CN (<2)


def _load_cn_bed(path):
    regions = []
    if not os.path.exists(path):
        return regions
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                regions.append((p[0], int(p[1]), int(p[2]), float(p[3])))
            except ValueError:
                continue
    return regions


def _load_loh_bed(path):
    regions = {}
    if not path or not os.path.exists(path):
        return regions
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                try:
                    regions.setdefault(p[0], []).append((int(p[1]), int(p[2])))
                except ValueError:
                    continue
    return regions


class GenomicContext:
    """Loads CN/LOH beds once; annotates loci. Pass loh_bed path (pipeline already has one)."""

    def __init__(self, loh_bed, gain_cn=GAIN_CN, loss_cn=LOSS_CN):
        self._gain = _load_cn_bed(gain_cn)
        self._loss = _load_cn_bed(loss_cn)
        self._loh = _load_loh_bed(loh_bed)

    def total_cn(self, chrom, pos):
        """integer total CN; 2.0 (neutral) if not in any gain/loss region."""
        for c, s, e, cn in self._gain:
            if c == chrom and s <= pos <= e:
                return cn
        for c, s, e, cn in self._loss:
            if c == chrom and s <= pos <= e:
                return cn
        return 2.0

    def loh_status(self, chrom, pos):
        for s, e in self._loh.get(chrom, []):
            if s <= pos <= e:
                return "LOH"
        return "nonLOH"

    def axis_validity(self, chrom, pos):
        """Which axis is the trustworthy 'somatic-controlled' primary at this locus.
        nonLOH -> HP-axis valid (germline-haplotype controlled).
        LOH    -> HP-axis INVALID (one haplotype lost); ALLELE-axis is primary (valid at VAF<1)."""
        loh = self.loh_status(chrom, pos)
        if loh == "LOH":
            return {"loh_status": "LOH", "primary_axis": "ALLELE", "hp_axis_valid": False,
                    "reason": "one haplotype lost in LOH -> HP1 vs HP1-1 not germline-vs-somatic"}
        return {"loh_status": "nonLOH", "primary_axis": "HP", "hp_axis_valid": True,
                "reason": "both haplotypes present -> HP-axis somatic-controlled"}

    def annotate(self, chrom, pos):
        av = self.axis_validity(chrom, pos)
        return {"chrom": chrom, "pos": pos, "total_cn": self.total_cn(chrom, pos), **av}
