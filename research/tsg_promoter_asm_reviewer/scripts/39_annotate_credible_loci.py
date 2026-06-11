#!/usr/bin/env python3
"""
39 — Gene + CpG-island annotation for the credible ASM loci (A4).

Annotates:
  - regime-A Tier-A pass loci (23, from regimeA_credible_probe.json)
  - LOH candidate_subclone loci (from loh_diagnostic_classifier.json)
with: nearest RefSeq gene + TSS distance + CpG-island context (island / shore<=2kb / open).

Inputs (local, downloaded 2026-05-27):
  data/hg38.ncbiRefSeq.gtf.gz, data/cpgIslandExt_hg38.txt.gz
No BAM access; deterministic. Output: genome_survey_v2/credible_loci_annotation.tsv
"""
import gzip, json
from pathlib import Path
from collections import defaultdict
import numpy as np

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
GTF = WS / "data/hg38.ncbiRefSeq.gtf.gz"
CPG = WS / "data/cpgIslandExt_hg38.txt.gz"
PROBE = WS / "genome_survey_v2/regimeA_credible_probe.json"
LOHC = WS / "genome_survey_v2/loh_diagnostic_classifier.json"
OUT = WS / "genome_survey_v2/credible_loci_annotation.tsv"
SHORE = 2000


def load_tss():
    """transcript TSS per gene from RefSeq GTF: {chrom: [(tss_pos, strand, gene)]}"""
    tss = defaultdict(list)
    seen = set()
    op = gzip.open(GTF, "rt")
    for line in op:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        chrom, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
        gene = None
        for kv in attr.split(";"):
            kv = kv.strip()
            if kv.startswith("gene_id"):
                gene = kv.split('"')[1] if '"' in kv else kv.split()[-1]
                break
        if gene is None:
            continue
        t = start if strand == "+" else end
        key = (chrom, t, gene)
        if key in seen:
            continue
        seen.add(key)
        tss[chrom].append((t, strand, gene))
    for c in tss:
        tss[c].sort()
    return tss


def load_cpg():
    """CpG islands: {chrom: [(start,end)]} sorted."""
    isl = defaultdict(list)
    op = gzip.open(CPG, "rt")
    for line in op:
        f = line.rstrip("\n").split("\t")
        # UCSC cpgIslandExt: bin chrom start end name ... -> chrom in f[1], start f[2], end f[3]
        if len(f) >= 4 and f[1].startswith("chr"):
            try:
                isl[f[1]].append((int(f[2]), int(f[3])))
            except ValueError:
                continue
        elif len(f) >= 3 and f[0].startswith("chr"):
            try:
                isl[f[0]].append((int(f[1]), int(f[2])))
            except ValueError:
                continue
    for c in isl:
        isl[c].sort()
    return isl


def nearest_tss(tss, chrom, pos):
    arr = tss.get(chrom, [])
    if not arr:
        return None, None, None
    best = None
    for t, strand, gene in arr:
        d = pos - t if strand == "+" else t - pos   # signed: + = downstream of TSS
        if best is None or abs(d) < abs(best[0]):
            best = (d, gene, strand)
    return best[1], best[0], best[2]


def cpg_context(isl, chrom, pos):
    arr = isl.get(chrom, [])
    mind = None
    for s, e in arr:
        if s <= pos <= e:
            return "island", 0
        d = min(abs(pos - s), abs(pos - e))
        if mind is None or d < mind:
            mind = d
    if mind is None:
        return "open", None
    return ("shore" if mind <= SHORE else "open"), mind


def main():
    print("loading RefSeq TSS + CpG islands...")
    tss = load_tss(); isl = load_cpg()
    print(f"  TSS chroms={len(tss)}, CpG-island chroms={len(isl)}")

    loci = []
    probe = json.load(open(PROBE))
    for r in probe['loci']:
        if r['pass_tierA']:
            c, p = r['locus'].split(":")
            loci.append(dict(locus=r['locus'], chrom=c, pos=int(p), set="regimeA_tierA",
                             axis=r['axis'], ari=r['ari'], delta=r['delta'],
                             germ_beta=r['germ_beta'], n_cpg=r['n_cpg_survey']))
    lohc = json.load(open(LOHC))
    for r in lohc['loci']:
        if r.get('klass') == 'candidate_subclone':
            loci.append(dict(locus=f"{r['chrom']}:{r['pos']}", chrom=r['chrom'], pos=int(r['pos']),
                             set="LOH_candidate_subclone", axis=r['axis'], ari=r['ari'],
                             delta=r['delta'], germ_beta=r['germ_beta'], n_cpg=r['n_cpg'],
                             cn_class=r.get('cn_class')))

    rows = []
    for L in loci:
        gene, dist, strand = nearest_tss(tss, L['chrom'], L['pos'])
        ctx, cpgd = cpg_context(isl, L['chrom'], L['pos'])
        L.update(nearest_gene=gene, tss_dist=dist, strand=strand, cpg_context=ctx, cpg_dist=cpgd)
        rows.append(L)

    # write TSV
    cols = ["set", "locus", "axis", "ari", "delta", "germ_beta", "n_cpg",
            "nearest_gene", "tss_dist", "strand", "cpg_context", "cpg_dist", "cn_class"]
    with open(OUT, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # summary
    print(f"\nannotated {len(rows)} credible loci -> {OUT}")
    print(f"{'set':<24}{'locus':<18}{'gene':<14}{'tss_dist':>9}{'cpg_ctx':>9}{'ari':>7}")
    for r in sorted(rows, key=lambda x: (x['set'], -(x['ari'] or 0))):
        print(f"{r['set']:<24}{r['locus']:<18}{str(r['nearest_gene'])[:13]:<14}"
              f"{('' if r['tss_dist'] is None else r['tss_dist']):>9}{r['cpg_context']:>9}{(r['ari'] or 0):>7.3f}")
    # promoter-proximal count
    prox = [r for r in rows if r['tss_dist'] is not None and abs(r['tss_dist']) <= 2000]
    isl_n = [r for r in rows if r['cpg_context'] in ('island', 'shore')]
    print(f"\npromoter-proximal (|TSS|<=2kb): {len(prox)}/{len(rows)}")
    print(f"CpG island/shore: {len(isl_n)}/{len(rows)}")


if __name__ == "__main__":
    main()
