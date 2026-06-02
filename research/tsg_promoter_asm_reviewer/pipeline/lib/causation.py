"""
causation.py — P3 mechanical-cis layer: does the somatic mutation DIRECTLY create or
destroy a CpG dinucleotide? This is the most direct cis-causation mechanism — the
mutation mechanically alters local methylability (a CpG can be methylated; non-CpG C
cannot). Pure deterministic reference-sequence fact (no read data, no probability).

Distinguishes two cis mechanisms for a T3 candidate:
  - mechanical  (DESTROYS_CpG / CREATES_CpG): the mutation itself is the cause — strong, simple.
  - regulatory  (NEUTRAL): the methylation change is NOT a CpG gain/loss; cause is regulatory
    (TFBS/promoter) or the subclone state — needs more evidence (location, subclone-partition).

BRCA2 chr13:32,315,128 G>A = NEUTRAL -> BRCA2's cis is regulatory, not mechanical (harder to prove).

Guardrail: characterization annotation only; emits no TP/FP score.
"""
import os
import pysam

REF_FA = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"


class RefContext:
    """Lazily opens the reference FASTA; annotates the mechanical-cis class of a variant."""

    def __init__(self, ref_fa=REF_FA):
        self._fa = pysam.FastaFile(ref_fa) if os.path.exists(ref_fa) else None

    @property
    def available(self):
        return self._fa is not None

    def mechanical_cis(self, chrom, pos, ref, alt):
        """pos is 1-based. Checks the CpG status of the dinucleotide spanning the variant
        on the reference vs the mutant base. Returns DESTROYS_CpG / CREATES_CpG / NEUTRAL / NA."""
        if self._fa is None:
            return "NA"
        ctx = self._fa.fetch(chrom, pos - 2, pos + 1).upper()  # [pos-1][variant(pos)][pos+1] (0-based pos-1)
        if len(ctx) < 3:
            return "NA"
        left, mid, right = ctx[0], ctx[1], ctx[2]
        ref_cpg = (mid == "C" and right == "G") or (left == "C" and mid == "G")
        alt_mid = alt.upper()
        alt_cpg = (alt_mid == "C" and right == "G") or (left == "C" and alt_mid == "G")
        if ref_cpg and not alt_cpg:
            return "DESTROYS_CpG"
        if alt_cpg and not ref_cpg:
            return "CREATES_CpG"
        return "NEUTRAL"


def causation_summary(mechanical, cis_tier):
    """combine mechanical-cis with the cis-ladder verdict into a causation evidence label.
    Only meaningful at T3 (cis-candidate). Lower tiers -> mechanism question is premature."""
    if cis_tier != "T3":
        return {"causation_mechanism": "n/a (not a cis-candidate)", "mechanical": mechanical}
    if mechanical in ("DESTROYS_CpG", "CREATES_CpG"):
        return {"causation_mechanism": f"mechanical ({mechanical}) — strongest cis evidence",
                "mechanical": mechanical,
                "next": ["confirm focal methylation change at the altered CpG"]}
    return {"causation_mechanism": "regulatory/subclone (not mechanical) — needs more evidence",
            "mechanical": mechanical,
            "next": ["diff-CpG location vs TSS/TFBS", "subclone-partition (mutation-specific vs subclone-general)"]}
