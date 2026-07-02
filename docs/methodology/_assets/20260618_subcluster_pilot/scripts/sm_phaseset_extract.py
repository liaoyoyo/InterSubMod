"""[L3a phase-set 抽取] 對每個 somatic sSNV, tumor pileup → ALT reads 的主要 PS tag + HP。
argv: chroms out_path。輸出 {locus: {ps, hp, n_alt_ps}}。供 sm_phaseset_analysis.py 評估 Tier-PS 延伸。
"""
import sys
import json
from collections import Counter
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20


def main(chroms, out_path):
    cen = json.load(open(f"{A}/sm_linkage_genomewide.json"))["census"]
    tb = pysam.AlignmentFile(TBAM, "rb")
    out = {}
    chset = set(chroms)
    loci = [(c["chrom"], c["pos"], c["ref"], c["alt"]) for c in cen.values()
            if c["chrom"] in chset and c.get("somatic") is True]
    loci.sort()
    for i, (chrom, pos, ref, alt) in enumerate(loci):
        ps_c = Counter()
        hp_c = Counter()
        for col in tb.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0, stepper="samtools"):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                aln = pr.alignment
                if aln.mapping_quality < MAPQ_MIN:
                    continue
                if aln.query_sequence[pr.query_position].upper() != alt:
                    continue
                if aln.has_tag("PS"):
                    ps_c[aln.get_tag("PS")] += 1
                if aln.has_tag("HP"):
                    hp_c[aln.get_tag("HP")] += 1
        ps = ps_c.most_common(1)[0][0] if ps_c else None
        hp = hp_c.most_common(1)[0][0] if hp_c else None
        out[f"{chrom}:{pos}"] = {"ps": ps, "hp": hp, "n_alt_ps": sum(ps_c.values())}
        if i % 2000 == 0:
            print(f"[{chroms[0]}..] {i}/{len(loci)}", flush=True)
    tb.close()
    json.dump(out, open(out_path, "w"), ensure_ascii=False)
    print(f"DONE {chroms} loci={len(loci)} with_ps={sum(1 for v in out.values() if v['ps'])} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
