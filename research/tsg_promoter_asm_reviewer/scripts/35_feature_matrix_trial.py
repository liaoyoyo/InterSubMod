#!/usr/bin/env python3
"""
35_feature_matrix_trial.py — TRIAL of a rich per-locus feature matrix for cis-ASM
characterization, on BRCA2 + 9 control loci. Purpose: see if multi-angle features
cross-verify the cis verdict, and RECORD which features are appropriate/informative.

Feature groups (per locus):
 (1) methylation signal  : dbeta_HP, dbeta_ALLELE, dbeta_5mC, dbeta_5hmC, percpg_maxabs,
                           maxabs_dist_to_var, entropy_HP1, entropy_HP11, dmr_span
 (2) cis evidence        : d_cis, d_drift, d_somatic, p_cis (from script 34 json)
 (3) subclone structure  : sil_HP1, sil_HP11, sil_HP2, delta_cohesion, kstar, crosstag_chi2_p
 (4) read/haplotype vars : n_HP1, n_HP11, n_HP2, vaf_obs, frac_phased, n_paired_cpg
 (5) genomic context     : mechanical_cis_cpg, total_cn, dist_to_tss, in_cpg_island, gene
 (6) quality/confound    : mapq_mean, nm_rate, softclip_rate, power_class, strand_balance

Read-only inputs: control_msa_out Level-1 + BAM + ref + CN beds + cpgIsland + GTF
                  + control_cohesion_cistest.json + control_loci_manifest.json
Output: genome_survey_v2/feature_matrix_trial.json (+ tsv)
"""
import gzip, json, os, glob, subprocess
import numpy as np
from collections import defaultdict, Counter
from scipy import stats
import pysam

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
GS = f"{ROOT}/genome_survey_v2"
L1 = glob.glob(f"{GS}/control_msa_out/*/level1_raw_methylation_details.tsv.gz")[0]
MAN = {m["locus"]: m for m in json.load(open(f"{GS}/control_loci_manifest.json"))}
CIS = json.load(open(f"{GS}/control_cohesion_cistest.json"))  # keyed by somatic_pos
cis_by_locus = {v["locus"]: v for v in CIS.values()}
REF = pysam.FastaFile("/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta")
TUMOR = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GAIN = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_gain_cn.bed"
LOSS = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_loss_cn.bed"
CPGI = f"{ROOT}/data/cpgIslandExt_hg38.txt.gz"
GTF = f"{ROOT}/data/hg38.ncbiRefSeq.gtf.gz"
THR = 0.5; MIN_N = 3

# ---------- load Level-1 ----------
# D[spos][(bam,hp,allele)][cpg][read] = max meth_call ; strand tracked
D = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
STRAND = defaultdict(Counter)
f = gzip.open(L1, "rt"); h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
for line in f:
    p = line.rstrip("\n").split("\t")
    try: mc = float(p[ix["meth_call"]])
    except ValueError: continue
    spos = p[ix["somatic_pos"]]; cpg = p[ix["methyl_pos"]]; rid = p[ix["read_id"]]
    cell = D[spos][(p[ix["bam_source_id"]], p[ix["haplotype_tag"]], p[ix["somatic_allele_type"]])][cpg]
    if rid not in cell or mc > cell[rid]: cell[rid] = mc
    STRAND[spos][p[ix["strand"]]] += 1

def beta_cpg(spos, bam, hps, alleles=None):
    bycpg = defaultdict(dict)
    for (b, hp, al), cpgd in D[spos].items():
        if b != bam or hp not in hps or (alleles and al not in alleles): continue
        for cpg, rd in cpgd.items():
            for rid, v in rd.items(): bycpg[cpg][rid] = v
    return {c: (np.mean([1 if x >= THR else 0 for x in d.values()]), len(d)) for c, d in bycpg.items()}

def read_matrix(spos, bam):
    mat = defaultdict(dict); info = {}
    for (b, hp, al), cpgd in D[spos].items():
        if b != bam: continue
        for cpg, rd in cpgd.items():
            for rid, v in rd.items():
                if cpg not in mat[rid] or v > mat[rid][cpg]: mat[rid][cpg] = v
                info[rid] = (hp, al)
    return mat, info

def dbeta(spos, germ_hps, som_hps, alleles_g=None, alleles_s=None):
    G = beta_cpg(spos, "tumor", germ_hps, alleles_g); S = beta_cpg(spos, "tumor", som_hps, alleles_s)
    sh = [c for c in set(G) & set(S) if G[c][1] >= MIN_N and S[c][1] >= MIN_N]
    if len(sh) < 4: return None, 0, None, None
    g = np.array([G[c][0] for c in sh]); s = np.array([S[c][0] for c in sh])
    d = s - g
    maxi = int(np.argmax(np.abs(d)))
    return float(np.mean(d)), len(sh), float(np.max(np.abs(d))), int(sh[maxi])

def entropy(spos, hps):
    B = beta_cpg(spos, "tumor", hps)
    vals = [b for b, n in B.values() if n >= MIN_N]
    if not vals: return None
    # mean per-CpG binary entropy
    es = []
    for v in vals:
        v = min(max(v, 1e-6), 1 - 1e-6)
        es.append(-(v*np.log2(v) + (1-v)*np.log2(1-v)))
    return round(float(np.mean(es)), 3)

def silhouette_tags(spos):
    mat, info = read_matrix(spos, "tumor")
    reads = [r for r in mat if len(mat[r]) >= 5]; n = len(reads)
    if n < 6: return {}, None, None
    dist = np.full((n, n), np.nan)
    for i in range(n):
        dist[i, i] = 0
        for j in range(i+1, n):
            sh = set(mat[reads[i]]) & set(mat[reads[j]])
            if len(sh) >= 5:
                dist[i, j] = dist[j, i] = np.mean([abs(mat[reads[i]][c]-mat[reads[j]][c]) for c in sh])
    tags = {r: info[r][0] for r in reads}
    res = {}
    labs = sorted(set(tags.values()))
    for li in labs:
        idin = [k for k in range(n) if tags[reads[k]] == li]
        if len(idin) < 3: continue
        sils = []
        for k in idin:
            a = np.nanmean([dist[k, m] for m in idin if m != k])
            bs = [np.nanmean([dist[k, m] for m in range(n) if tags[reads[m]] == lj]) for lj in labs if lj != li]
            bs = [b for b in bs if not np.isnan(b)]
            if bs and not np.isnan(a): sils.append((min(bs)-a)/max(a, min(bs)))
        if sils: res[li] = round(float(np.mean(sils)), 3)
    # kstar: global silhouette over k=1..4 (unsupervised) — FIX: min-n gate (small-n unstable)
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score
    MIN_SOM = 10
    n_som = sum(1 for r in reads if tags[r] in ("1-1", "2-1"))
    if n_som < MIN_SOM:
        return res, None, None  # too few SOMATIC reads -> k-selection/crosstag for cis-ASM unstable
    Dfill = np.where(np.isnan(dist), np.nanmax(dist), dist)
    best_k, best_s = 1, -1
    for k in (2, 3, 4):
        if n <= k: break
        try:
            lab = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage="average").fit_predict(Dfill)
            s = silhouette_score(Dfill, lab, metric="precomputed")
            if s > best_s: best_s, best_k = s, k
        except Exception: pass
    if best_s < 0.10: best_k = 1  # weak -> single cluster
    # cross-tag membership: tag(germline/somatic) x best_k cluster chi2
    cm_p = None
    if best_k >= 2:
        try:
            lab = AgglomerativeClustering(n_clusters=best_k, metric="precomputed", linkage="average").fit_predict(Dfill)
            grp = ["som" if tags[reads[k]] in ("1-1", "2-1") else "germ" for k in range(n)]
            tab = np.zeros((2, best_k))
            for k in range(n): tab[0 if grp[k] == "germ" else 1, lab[k]] += 1
            tab = tab[:, tab.sum(0) > 0]
            if tab.shape[1] >= 2: cm_p = round(float(stats.chi2_contingency(tab)[1]), 4)
        except Exception: pass
    return res, best_k, cm_p

# ---------- annotations ----------
def cn_of(chrom, pos):
    for src, sign in ((GAIN, ">2"), (LOSS, "<2")):
        try:
            out = subprocess.run(f"awk -F'\\t' '$1==\"{chrom}\" && $2<={pos} && $3>={pos} {{print $4; exit}}' {src}",
                                 shell=True, capture_output=True, text=True, timeout=30).stdout.strip()
            if out:
                try: return float(out)
                except ValueError: return out
        except Exception: pass
    return 2.0  # neutral default

def mechanical_cis(chrom, pos, ref, alt):
    """does somatic mutation create/destroy a CpG dinucleotide?"""
    ctx = REF.fetch(chrom, pos-2, pos+1).upper()  # [pos-1][pos][pos+1] (0-based pos-1=ref)
    if len(ctx) < 3: return "NA"
    left, mid, right = ctx[0], ctx[1], ctx[2]
    ref_cpg = (mid == "C" and right == "G") or (left == "C" and mid == "G")
    alt_mid = alt.upper()
    alt_cpg = (alt_mid == "C" and right == "G") or (left == "C" and alt_mid == "G")
    if ref_cpg and not alt_cpg: return "DESTROYS_CpG"
    if alt_cpg and not ref_cpg: return "CREATES_CpG"
    return "NEUTRAL"

# CpG island load (small)
CPG_ISL = []
with gzip.open(CPGI, "rt") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        # UCSC cpgIslandExt: bin chrom start end name ...
        try: CPG_ISL.append((p[1], int(p[2]), int(p[3])))
        except (IndexError, ValueError): continue
def cpg_island_class(chrom, pos):
    for c, s, e in CPG_ISL:
        if c == chrom:
            if s <= pos <= e: return "island"
            if s-2000 <= pos <= e+2000: return "shore"
            if s-4000 <= pos <= e+4000: return "shelf"
    return "open_sea"

# TSS / gene from GTF (filter per chrom for speed)
def tss_gene(chrom, pos):
    try:
        out = subprocess.run(
            f"zcat {GTF} | awk -F'\\t' '$1==\"{chrom}\" && $3==\"transcript\"' | head -200000 | "
            f"awk -F'\\t' '{{tss=($7==\"+\")?$4:$5; d=(tss>{pos})?tss-{pos}:{pos}-tss; print d\"\\t\"$9}}' | sort -n | head -1",
            shell=True, capture_output=True, text=True, timeout=90).stdout.strip()
        if out:
            dist = int(out.split("\t")[0]); attr = out.split("\t")[-1]
            gname = "NA"
            for field in attr.split(";"):
                if "gene_name" in field: gname = field.split('"')[1] if '"' in field else "NA"
            return dist, gname
    except Exception: pass
    return None, "NA"

def power_class(n_shared, p):
    floor = min(1.0, 2.0/(2**n_shared)) if n_shared <= 25 else 0.0
    BONF = 0.05/51091
    if p is not None and p < BONF: return "DETECTED"
    if floor >= BONF: return "UNDERPOWERED"
    return "NULL"

# BAM features per locus: 5mC/5hmC dbeta + MAPQ/NM/softclip
bam = pysam.AlignmentFile(TUMOR)
def bam_features(chrom, spos, ref, alt):
    by = defaultdict(lambda: defaultdict(dict))  # (hp,mod) -> cpg -> {read:prob}
    mapq = defaultdict(list); nm = defaultdict(list); sc = defaultdict(list)
    for r in bam.fetch(chrom, spos-1000, spos+1000):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        try: hp = str(r.get_tag('HP'))
        except KeyError: hp = 'untag'
        # allele at spos
        al = 'na'
        q2r = dict((q, rp) for q, rp in r.get_aligned_pairs(matches_only=True))
        r2q = {v: k for k, v in q2r.items()}
        if (spos-1) in r2q:
            b = r.query_sequence[r2q[spos-1]].upper(); al = 'alt' if b == alt.upper() else ('ref' if b == ref.upper() else 'oth')
        grp = hp
        mapq[grp].append(r.mapping_quality)
        try: nm[grp].append(r.get_tag('NM')/max(1, r.query_length))
        except KeyError: pass
        cig = r.cigartuples or []
        sclen = sum(l for op, l in cig if op == 4)
        sc[grp].append(sclen/max(1, r.query_length))
        mb = r.modified_bases
        if not mb: continue
        for (canon, strand, mod), calls in mb.items():
            mc = 'm' if str(mod) in ('m', '0') else ('h' if str(mod) in ('h', '17') else str(mod))
            if mc not in ('m', 'h'): continue
            for qpos, qual in calls:
                rp = q2r.get(qpos)
                if rp is None: continue
                if abs(rp - (spos-1)) > 1000: continue  # FIX: window-restrict (ONT reads span 10-30kb)
                by[hp][rp].setdefault(r.query_name, {})[mc] = qual/255.0
    def beta_paired(hp, mode):
        # per-CpG frac methylated, >=MIN_N reads; mode 'm'/'h'/'any'(max over mods per read)
        out = {}
        for cpg, rd in by.get(hp, {}).items():
            vals = []
            for rec in rd.values():
                if mode == 'any':
                    v = max(rec.get('m', 0.0), rec.get('h', 0.0))
                elif mode in rec:
                    v = rec[mode]
                else:
                    continue
                vals.append(1 if v >= 0.5 else 0)
            if len(vals) >= MIN_N:
                out[cpg] = (np.mean(vals), len(vals))
        return out
    def dbeta_mod(mode):
        # FIX: PAIRED per-CpG (same method as dbeta_HP) so 5mC/5hmC/any are comparable
        g = beta_paired("1", mode); s = beta_paired("1-1", mode)
        sh = [c for c in set(g) & set(s) if g[c][1] >= MIN_N and s[c][1] >= MIN_N]
        return round(float(np.mean([s[c][0]-g[c][0] for c in sh])), 3) if len(sh) >= 4 else None
    def m(d, grps):
        vals = [x for gg in grps for x in d.get(gg, [])]
        return round(float(np.mean(vals)), 3) if vals else None
    return {"dbeta_5mC": dbeta_mod("m"), "dbeta_5hmC": dbeta_mod("h"), "dbeta_any_bam": dbeta_mod("any"),
            "mapq_mean": m(mapq, ["1", "1-1", "2"]), "nm_rate": m(nm, ["1", "1-1", "2"]),
            "softclip_rate": m(sc, ["1", "1-1", "2"])}

# ---------- assemble per locus ----------
rows = []
for spos in sorted(D, key=lambda s: int(s)):
    locus = None
    for k, m in MAN.items():
        if m["locus"].split(":")[1] == spos: locus = m; break
    if locus is None: continue
    chrom = locus["locus"].split(":")[0]; ref = locus["ref"]; alt = locus["alt"]
    sil, kstar, cm_p = silhouette_tags(spos)
    dhp, nhp, mx_hp, mxpos = dbeta(spos, {"1"}, {"1-1"})
    dal, nal, _, _ = dbeta(spos, {"1", "2"}, {"1", "2"}, alleles_g={"ref"}, alleles_s={"alt"})
    cis = cis_by_locus.get(locus["locus"], {}).get("cis", {})
    # read counts
    cnt = Counter()
    for (b, hp, al), cpgd in D[spos].items():
        if b != "tumor": continue
        for cpg, rd in cpgd.items(): cnt[hp] += 0
    # count unique reads per tag
    tagreads = defaultdict(set)
    for (b, hp, al), cpgd in D[spos].items():
        if b != "tumor": continue
        for cpg, rd in cpgd.items():
            for rid in rd: tagreads[hp].add(rid)
    n1 = len(tagreads.get("1", set())); n11 = len(tagreads.get("1-1", set())); n2 = len(tagreads.get("2", set()))
    alt_reads = set(); ref_reads = set()
    for (b, hp, al), cpgd in D[spos].items():
        if b != "tumor": continue
        for cpg, rd in cpgd.items():
            for rid in rd: (alt_reads if al == "alt" else ref_reads if al == "ref" else set()).add(rid)
    vaf = round(len(alt_reads)/max(1, len(alt_reads)+len(ref_reads)), 3)
    untag = len(tagreads.get("untag", set()))
    tot = sum(len(v) for v in tagreads.values())
    bf = bam_features(chrom, int(spos), ref, alt)
    cn = cn_of(chrom, int(spos)); mech = mechanical_cis(chrom, int(spos), ref, alt)
    isl = cpg_island_class(chrom, int(spos)); dtss, gene = tss_gene(chrom, int(spos))
    n_sh = cis.get("n_shared_cpg", 0)
    row = {
        "tag": locus["tag"], "locus": locus["locus"], "cls": locus["cls"], "loh": locus["loh"],
        # (1) methylation
        "dbeta_HP": round(dhp, 3) if dhp is not None else None, "dbeta_ALLELE": round(dal, 3) if dal is not None else None,
        "dbeta_5mC": bf["dbeta_5mC"], "dbeta_5hmC": bf["dbeta_5hmC"], "dbeta_any_bam": bf["dbeta_any_bam"],
        "percpg_maxabs": round(mx_hp, 3) if mx_hp is not None else None,
        "maxabs_dist_to_var": (mxpos-int(spos)) if mxpos is not None else None,
        "entropy_HP1": entropy(spos, {"1"}), "entropy_HP11": entropy(spos, {"1-1"}),
        # (2) cis
        "d_cis": cis.get("d_cis_C_minus_normalA"), "d_drift": cis.get("d_drift_B_minus_normalA"),
        "d_somatic": cis.get("d_somatic_C_minus_B"), "p_cis": cis.get("p_cis"),
        # (3) subclone
        "sil_HP1": sil.get("1"), "sil_HP11": sil.get("1-1"), "sil_HP2": sil.get("2"),
        "delta_cohesion": round(sil.get("1-1", 0)-sil.get("1", 0), 3) if sil.get("1-1") is not None and sil.get("1") is not None else None,
        "kstar": kstar, "crosstag_chi2_p": cm_p,
        # (4) read/hap
        "n_HP1": n1, "n_HP11": n11, "n_HP2": n2,
        "n_som": len(tagreads.get("1-1", set())) + len(tagreads.get("2-1", set())), "vaf_obs": vaf,
        "frac_untag": round(untag/max(1, tot), 3), "n_paired_cpg": n_sh,
        # (5) genomic context
        "mechanical_cis": mech, "total_cn": cn, "dist_to_tss": dtss, "cpg_island": isl, "gene": gene,
        # (6) quality
        "mapq_mean": bf["mapq_mean"], "nm_rate": bf["nm_rate"], "softclip_rate": bf["softclip_rate"],
        "power_class": power_class(n_sh, cis.get("p_cis")),
        "strand_balance": round(min(STRAND[spos].values())/max(1, sum(STRAND[spos].values())), 3) if STRAND[spos] else None,
    }
    rows.append(row)

json.dump(rows, open(f"{GS}/feature_matrix_trial.json", "w"), indent=1, ensure_ascii=False)
# tsv
cols = list(rows[0].keys())
with open(f"{GS}/feature_matrix_trial.tsv", "w") as fo:
    fo.write("\t".join(cols)+"\n")
    for r in rows: fo.write("\t".join(str(r[c]) for c in cols)+"\n")
print("WROTE feature_matrix_trial.json/tsv —", len(rows), "loci x", len(cols), "features")
# print compact
key = ["tag", "cls", "loh", "dbeta_HP", "dbeta_5mC", "dbeta_5hmC", "d_cis", "d_drift", "p_cis",
       "sil_HP11", "delta_cohesion", "kstar", "crosstag_chi2_p", "vaf_obs", "n_HP11",
       "mechanical_cis", "total_cn", "dist_to_tss", "cpg_island", "gene", "mapq_mean", "nm_rate", "power_class"]
print("\t".join(key))
for r in rows:
    print("\t".join(str(r.get(k)) for k in key))
