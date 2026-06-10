#!/usr/bin/env python3
"""Large-scale methodology verification (816 loci) of 3 questions:
   Q1 rate(Fisher per-CpG) vs structure(read-read d_between-d_within) — how many loci are
      'structure-visible but rate-weak' (eye-sees-structure-Fisher-misses, generalised);
   Q2 normal-HP-tag reference value — does adding normal HP1-vs-HP2 germline-allelic baseline
      separate tumor-specific from germline-baseline ASM, and validate difference regions;
   Q3 Fisher sensitivity — fraction of large-structure loci the per-CpG rate test misses.
   Proxy for ISM PERMANOVA = d_between-d_within on NHD read-read distance (LabelTest Δ).
   Runs on tsg Level-1 caches (per-read per-CpG meth_call) + cis_scan_full.json tiers."""
import gzip, csv, glob, json, os
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact
rng = np.random.default_rng(20260611)
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod"
TSG = f"{ROOT}/research/tsg_promoter_asm_reviewer"
OUT = f"{ROOT}/docs/method_comparison/20260609_ism_vs_external_methylation_tools/benchmark/runs"
THR = 0.5; MIN_N = 5; MIN_CPG = 4; MAXREAD = 80; NPERM = 99

# tier lookup
cis = json.load(open(f"{TSG}/genome_survey_v2/cis_scan_full.json"))
cis = cis if isinstance(cis, list) else list(cis.values())
tier = {r["locus"].split(":")[1]: r.get("cis_tier", "?") for r in cis}

def load(level1):
    # reads[(bam,hp)][read_id][cpg]=0/1 ; only 5mC ('m')
    R = defaultdict(lambda: defaultdict(dict))
    with gzip.open(level1, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["mod_code"] != "m":
                continue
            R[(row["bam_source"], row["haplotype_tag"])][row["read_id"]][int(row["methyl_pos"])] = \
                1 if float(row["meth_call"]) >= THR else 0
    return R

def beta_cpg(group):  # group = {read:{cpg:0/1}} -> {cpg:(meth,unmeth)}
    bc = defaultdict(lambda: [0, 0])
    for rd in group.values():
        for cpg, b in rd.items():
            bc[cpg][0 if b == 1 else 1] += 1
    return bc

def paired_absdbeta(gA, gB):
    a = beta_cpg(gA); b = beta_cpg(gB)
    shared = [c for c in set(a) & set(b) if sum(a[c]) >= MIN_N and sum(b[c]) >= MIN_N]
    if len(shared) < MIN_CPG:
        return None, 0, None, None
    dabs = []; fsig = 0; ntum_spec = []
    for c in shared:
        ba = a[c][0] / sum(a[c]); bb = b[c][0] / sum(b[c]); d = ba - bb
        dabs.append(abs(d))
        _, p = fisher_exact([[a[c][0], a[c][1]], [b[c][0], b[c][1]]])
        if p < 0.05: fsig += 1
        ntum_spec.append((c, abs(d)))
    return float(np.mean(dabs)), len(shared), fsig / len(shared), {c: d for c, d in ntum_spec}

def struct_delta(gA, gB):
    # read-read NHD; d_between - d_within ; light 99-perm p
    reads = [(r, 0) for r in gA] + [(r, 1) for r in gB]
    prof = {**{r: gA[r] for r in gA}, **{r: gB[r] for r in gB}}
    if len(gA) < MIN_N or len(gB) < MIN_N:
        return None, None
    # subsample
    if len(reads) > 2 * MAXREAD:
        idxA = [r for r in gA]; idxB = [r for r in gB]
        rng.shuffle(idxA); rng.shuffle(idxB)
        reads = [(r, 0) for r in idxA[:MAXREAD]] + [(r, 1) for r in idxB[:MAXREAD]]
    ids = [r for r, _ in reads]; lab = np.array([l for _, l in reads]); n = len(ids)
    D = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            sh = set(prof[ids[i]]) & set(prof[ids[j]])
            if len(sh) >= MIN_N:
                d = np.mean([abs(prof[ids[i]][c] - prof[ids[j]][c]) for c in sh])
                D[i, j] = D[j, i] = d
    def delta(lb):
        w = []; bw = []
        for i in range(n):
            for j in range(i + 1, n):
                if np.isnan(D[i, j]): continue
                (w if lb[i] == lb[j] else bw).append(D[i, j])
        if not w or not bw: return None
        return np.mean(bw) - np.mean(w)
    obs = delta(lab)
    if obs is None: return None, None
    ge = 1
    for _ in range(NPERM):
        sh = lab.copy(); rng.shuffle(sh)
        d = delta(sh)
        if d is not None and d >= obs: ge += 1
    return float(obs), float(ge / (NPERM + 1))

rows = []
caches = sorted(glob.glob(f"{TSG}/pipeline/cache/level1/*.tsv.gz"))
for k, cache in enumerate(caches):
    base = os.path.basename(cache)
    pos = base.split("_")[1]
    try:
        R = load(cache)
    except Exception:
        continue
    tumHP1 = R.get(("tumor", "1"), {}); tumHP11 = R.get(("tumor", "1-1"), {})
    norHP1 = R.get(("normal", "1"), {}); norHP2 = R.get(("normal", "2"), {})
    # tumor somatic axis
    t_abs, t_ncpg, t_fsig, t_cpgd = paired_absdbeta(tumHP1, tumHP11)
    # normal germline-allelic baseline (HP1 vs HP2)
    n_abs, n_ncpg, n_fsig, n_cpgd = paired_absdbeta(norHP1, norHP2)
    # structure (tumor HP1 vs HP1-1)
    sdelta, sp = struct_delta(tumHP1, tumHP11)
    # normal-referenced: tumor-specific CpGs = tumor diff large but normal flat
    ntum_spec = None
    if t_cpgd and n_cpgd is not None:
        shared = set(t_cpgd) & set(n_cpgd)
        if shared:
            ntum_spec = sum(1 for c in shared if t_cpgd[c] >= 0.3 and n_cpgd.get(c, 0) < 0.2)
    rows.append(dict(pos=pos, tier=tier.get(pos, "?"),
                     t_ncpg=t_ncpg, t_absdbeta=t_abs, t_frac_fisher_sig=t_fsig,
                     struct_delta=sdelta, struct_p=sp,
                     n_normal_cpg=n_ncpg, normal_absdbeta=n_abs,
                     normal_ref_net=(None if (t_abs is None or n_abs is None) else round(t_abs - n_abs, 4)),
                     n_tumor_specific_cpg=ntum_spec,
                     has_normal=(n_ncpg or 0) >= MIN_CPG))
    if (k + 1) % 100 == 0:
        print(f"...{k+1}/{len(caches)} processed", flush=True)

# write per-locus
with open(f"{OUT}/largescale_perlocus.tsv", "w") as f:
    cols = ["pos", "tier", "t_ncpg", "t_absdbeta", "t_frac_fisher_sig", "struct_delta", "struct_p",
            "n_normal_cpg", "normal_absdbeta", "normal_ref_net", "n_tumor_specific_cpg", "has_normal"]
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c)) for c in cols) + "\n")

# aggregate
def num(x): return x is not None
testable = [r for r in rows if num(r["t_absdbeta"]) and num(r["struct_delta"])]
struct_sig = [r for r in testable if r["struct_delta"] > 0 and r["struct_p"] is not None and r["struct_p"] <= 0.05]
rate_weak = lambda r: (r["t_frac_fisher_sig"] is not None and r["t_frac_fisher_sig"] < 0.05)
with_normal = [r for r in rows if r["has_normal"] and num(r["t_absdbeta"]) and num(r["normal_absdbeta"])]
summary = {
    "n_loci_total": len(rows),
    "n_testable_tumor_axis": len(testable),
    "Q1_structure_sig (delta>0,p<=0.05)": len(struct_sig),
    "Q1_structure_sig_BUT_rate_weak (Fisher frac<5%)": sum(1 for r in struct_sig if rate_weak(r)),
    "Q1_pct_structure_sig_rate_weak": round(100 * sum(1 for r in struct_sig if rate_weak(r)) / max(1, len(struct_sig)), 1),
    "Q2_n_with_normal_HP1HP2": len(with_normal),
    "Q2_tumor_absdbeta_median": round(float(np.median([r["t_absdbeta"] for r in with_normal])), 4) if with_normal else None,
    "Q2_normal_absdbeta_median": round(float(np.median([r["normal_absdbeta"] for r in with_normal])), 4) if with_normal else None,
    "Q2_loci_tumor>normal (tumor-specific by net)": sum(1 for r in with_normal if r["normal_ref_net"] and r["normal_ref_net"] > 0.05),
    "Q2_loci_tumor<=normal (germline-baseline-driven)": sum(1 for r in with_normal if r["normal_ref_net"] is not None and r["normal_ref_net"] <= 0.05),
    "Q3_struct_sig_loci": len(struct_sig),
    "Q3_of_those_rate_strong (Fisher frac>=10%)": sum(1 for r in struct_sig if r["t_frac_fisher_sig"] is not None and r["t_frac_fisher_sig"] >= 0.10),
}
json.dump({"summary": summary}, open(f"{OUT}/largescale_summary.json", "w"), indent=2)
print(json.dumps(summary, indent=2))
print("DONE", len(rows), "loci ->", f"{OUT}/largescale_perlocus.tsv")
