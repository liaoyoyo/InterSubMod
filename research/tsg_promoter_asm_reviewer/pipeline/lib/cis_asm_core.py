"""
cis_asm_core.py — versioned, importable core of the cis-ASM characterization pipeline.

Refactored (behavior-preserved) from exploratory scripts 18 (dual-axis pivot),
34 (normal-anchored cis-test), 35 (per-tag cohesion). Level-1-based functions only
(self-contained, no BAM/annotation needed) so the BRCA2 regression test is hermetic.

Guardrail: cis-ASM CHARACTERIZATION only — computes NO TP/FP discriminator score,
feeds NO F1 filter. HP-axis (HP1 vs HP1-1, same germline haplotype) is the
somatic-controlled axis; ALLELE-axis (ALT vs REF) is baseline-confounded.

Frozen conventions (golden test depends on these — change => [golden-update] commit):
  THR=0.5  MIN_N=3 (reads/CpG/group)  MIN_PAIR=4 (paired CpGs)
"""
import gzip
import numpy as np
from collections import defaultdict
from scipy import stats

THR = 0.5
MIN_N = 3
MIN_PAIR = 4

L1_COLS = ("somatic_pos", "methyl_pos", "haplotype_tag", "meth_call",
           "bam_source_id", "somatic_allele_type", "read_id")


def load_level1(path):
    """Level-1 TSV(.gz) -> D[spos][(bam,hp,allele)][cpg][read_id] = max meth_call.

    5mC/5hmC are emitted as two unlabeled rows per (read,CpG); we max-collapse to
    'any-modification' here (same as script 18). Mod-type separation needs the BAM
    (see PoC 36 — read-level 5hmC clustering is not viable; 5hmC too sparse).
    """
    D = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for col in L1_COLS:
            if col not in ix:
                raise KeyError(f"Level-1 missing column {col!r} in {path}")
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                mc = float(p[ix["meth_call"]])
            except (ValueError, IndexError):
                continue
            cell = D[p[ix["somatic_pos"]]][(p[ix["bam_source_id"]],
                                            p[ix["haplotype_tag"]],
                                            p[ix["somatic_allele_type"]])][p[ix["methyl_pos"]]]
            rid = p[ix["read_id"]]
            if rid not in cell or mc > cell[rid]:
                cell[rid] = mc
    return D


def _beta_cpg(D, spos, bam, hps, alleles=None):
    """per-CpG (frac reads meth>=THR, n_reads) for reads matching bam/hps/alleles."""
    bycpg = defaultdict(dict)
    for (b, hp, al), cpgd in D[spos].items():
        if b != bam or hp not in hps or (alleles is not None and al not in alleles):
            continue
        for cpg, rd in cpgd.items():
            for rid, v in rd.items():
                bycpg[cpg][rid] = v
    return {c: (float(np.mean([1 if x >= THR else 0 for x in d.values()])), len(d))
            for c, d in bycpg.items()}


def dbeta_axis(D, spos, germ_hps, som_hps, alleles_g=None, alleles_s=None, bam="tumor"):
    """paired per-CpG mean delta-beta (som - germ). HP-axis: germ={'1'}, som={'1-1'}.
    ALLELE-axis: hps={'1','2'}, alleles_g={'ref'}, alleles_s={'alt'}. None if untestable."""
    G = _beta_cpg(D, spos, bam, germ_hps, alleles_g)
    S = _beta_cpg(D, spos, bam, som_hps, alleles_s)
    sh = [c for c in set(G) & set(S) if G[c][1] >= MIN_N and S[c][1] >= MIN_N]
    if len(sh) < MIN_PAIR:
        return None
    return round(float(np.mean([S[c][0] - G[c][0] for c in sh])), 3)


def cis_test(D, spos):
    """normal-anchored cis-test (HP1 axis): A=normal HP1, B=tumor HP1, C=tumor HP1-1.
    d_cis = C-A (somatic vs pre-mutation baseline); d_drift = B-A (germline drift);
    d_somatic = C-B (within-tumor ASM). cis-candidate <=> |d_cis| sig AND >> |d_drift|."""
    A = _beta_cpg(D, spos, "normal", {"1"})
    B = _beta_cpg(D, spos, "tumor", {"1"})
    C = _beta_cpg(D, spos, "tumor", {"1-1"})
    sh = [c for c in set(A) & set(B) & set(C)
          if A[c][1] >= MIN_N and B[c][1] >= MIN_N and C[c][1] >= MIN_N]
    if len(sh) < MIN_PAIR:
        return {"testable": False, "n_shared_cpg": len(sh)}
    a = np.array([A[c][0] for c in sh])
    b = np.array([B[c][0] for c in sh])
    c_ = np.array([C[c][0] for c in sh])

    def wp(x, y):
        try:
            return round(float(stats.wilcoxon(x, y).pvalue), 4)
        except Exception:
            return None
    return {"testable": True, "n_shared_cpg": len(sh),
            "d_cis": round(float(np.mean(c_ - a)), 4),
            "d_drift": round(float(np.mean(b - a)), 4),
            "d_somatic": round(float(np.mean(c_ - b)), 4),
            "p_cis": wp(c_, a)}


def cohesion_per_tag(D, spos, min_shared=5):
    """per-HP-tag silhouette (within vs between-group) on observed-only read distance.
    high silhouette = somatic reads form a cohesive cluster. Returns {hp: silhouette}."""
    mat = defaultdict(dict)
    info = {}
    for (b, hp, al), cpgd in D[spos].items():
        if b != "tumor":
            continue
        for cpg, rd in cpgd.items():
            for rid, v in rd.items():
                if cpg not in mat[rid] or v > mat[rid][cpg]:
                    mat[rid][cpg] = v
                info[rid] = hp
    reads = [r for r in mat if len(mat[r]) >= min_shared]
    n = len(reads)
    if n < 6:
        return {}
    dist = np.full((n, n), np.nan)
    for i in range(n):
        dist[i, i] = 0.0
        for j in range(i + 1, n):
            sh = set(mat[reads[i]]) & set(mat[reads[j]])
            if len(sh) >= min_shared:
                dist[i, j] = dist[j, i] = float(np.mean([abs(mat[reads[i]][c] - mat[reads[j]][c]) for c in sh]))
    tags = {r: info[r] for r in reads}
    labs = sorted(set(tags.values()))
    res = {}
    for li in labs:
        idin = [k for k in range(n) if tags[reads[k]] == li]
        if len(idin) < 3:
            continue
        sils = []
        for k in idin:
            a = np.nanmean([dist[k, m] for m in idin if m != k])
            bs = [np.nanmean([dist[k, m] for m in range(n) if tags[reads[m]] == lj])
                  for lj in labs if lj != li]
            bs = [x for x in bs if not np.isnan(x)]
            if bs and not np.isnan(a):
                sils.append((min(bs) - a) / max(a, min(bs)))
        if sils:
            res[li] = round(float(np.mean(sils)), 3)
    return res


def score_locus(D, spos):
    """unified per-locus core record (Level-1-based features only)."""
    coh = cohesion_per_tag(D, spos)
    ct = cis_test(D, spos)
    return {
        "somatic_pos": spos,
        "dbeta_HP": dbeta_axis(D, spos, {"1"}, {"1-1"}),
        "dbeta_ALLELE": dbeta_axis(D, spos, {"1", "2"}, {"1", "2"},
                                   alleles_g={"ref"}, alleles_s={"alt"}),
        "sil_HP1": coh.get("1"), "sil_HP11": coh.get("1-1"), "sil_HP2": coh.get("2"),
        "delta_cohesion": (round(coh["1-1"] - coh["1"], 3)
                           if coh.get("1-1") is not None and coh.get("1") is not None else None),
        "cis": ct,
    }
