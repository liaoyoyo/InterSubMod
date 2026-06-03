#!/usr/bin/env python3
"""
57 - Cross-sample verification of HCC1395 key ASM positions (A pilot, targeted).

QUESTION (user, 2026-06-03)
  Do the OTHER 5 samples (COLO829, H1437, H2009, HCC1937, HCC1954) show the same
  allele-specific-methylation (ASM) pattern at HCC1395's 38 characterized key
  positions, and to what degree? Note: cancers carry PRIVATE somatic mutations, so
  same-position somatic recurrence is biologically expected to be RARE. The honest
  decomposition of "similar pattern" is therefore:
    (a) Same somatic event at the same position  -> almost never (private).
    (b) Germline / imprinting ASM (HP1 vs HP2) at the same position -> CAN recur
        across samples (shared germline / imprinted loci, e.g. IGF2).
    (c) Somatic sub-haplotype ASM (HPx vs HPx-1) -> only where THAT sample independently
        carries a somatic event at the locus.

  This script extracts, per (sample, key-position): somatic status, HP-group coverage,
  and per-HP-group 5mC beta, so the synthesis can classify each position as
  HCC1395-private-somatic-ASM vs recurrent-germline/imprinting-ASM vs no-signal/LOH.

METHOD (reuses VALIDATED pysam logic from 54_modkit_extract_crossval.py, which PASSED
  the modkit v0.6.3 cross-val gate at Pearson=1.0)
  - per read: r.modified_bases, code in ('m','h'), get_aligned_pairs query->ref, ML/255.
  - 5mC ('m') is the canonical primary signal (modkit lesson: MAX-collapse loses 5hmC;
    here we keep 5mC and 5hmC SEPARATE). per-read beta = mean 5mC prob over CpGs in
    window. group reads by HP:Z tag string ('1','2','1-1','2-1',...).
  - WINDOW = +/-600 bp (matches 54).
  - MIN_N = 3 reads per HP group to report a beta.
  - somatic_in_sample: this sample's own somatic_pass.vcf.gz has a PASS record within
    +/-SOM_TOL bp of the position.
  - has_somatic_subhap: any HP:Z tag containing '-' (e.g. '1-1') at the locus.

AXES
  - hp_axis_delta (somatic-CONTROLLED, comparable to HCC1395's dbeta_HP):
      for credible axis 'HPk_vs_HPk-1', delta = beta(HPk-1) - beta(HPk), if both present.
  - germline_delta (HP1 vs HP2; imprinting / baseline allelic):
      delta = beta(HP2) - beta(HP1), if both present.

USAGE
  python3 57_cross_sample_keypos.py <SAMPLE>      # one of the 6 sample names
  (writes genome_survey_v2/cn_confound/cross_sample/<SAMPLE>_keypos.json)
"""
import os, sys, json
os.environ.setdefault("TMPDIR", "/big7_disk")
import numpy as np
import pysam

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
GS = f"{ROOT}/genome_survey_v2"
CN = f"{GS}/cn_confound"
OUTDIR = f"{CN}/cross_sample"
CRED = f"{GS}/credible_loci_annotation.tsv"
REF = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"

WINDOW = 600
MIN_N = 3
SOM_TOL = 5          # bp tolerance for "somatic at this position"
MOD_THR = 0.5        # P>=0.5 => methylated call (for frac-methylated)

# sample -> complete_matrix run dir name (under canonical/<sample>/paired_full/)
RUNDIR = {
    "HCC1395": "20260314_HCC1395_paired_full_full_complete_matrix",
    "COLO829": "20260315_COLO829_paired_full_full_complete_matrix",
    "H1437":   "20260315_H1437_paired_full_full_complete_matrix",
    "H2009":   "20260315_H2009_paired_full_full_complete_matrix",
    "HCC1937": "20260315_HCC1937_paired_full_full_complete_matrix",
    "HCC1954": "20260315_HCC1954_paired_full_full_complete_matrix",
}
CANCER = {
    "HCC1395": "breast (ductal carcinoma; SEQC2 ref)",
    "COLO829": "melanoma",
    "H1437":   "lung adenocarcinoma",
    "H2009":   "lung adenocarcinoma",
    "HCC1937": "breast (BRCA1-mutant)",
    "HCC1954": "breast",
}


def paths_for(sample):
    base = (f"/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{sample}/"
            f"paired_full/{RUNDIR[sample]}/longphase_s")
    return (f"{base}/{sample}_tagged.bam",
            f"{base}/somatic_pass.vcf.gz")


def load_positions():
    """37 credible loci + BRCA2(chr13:32315128). Returns list of dicts."""
    out = []
    with open(CRED) as f:
        header = f.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            chrom, pos = p[ci["locus"]].split(":")
            out.append(dict(
                chrom=chrom, pos=int(pos), axis=p[ci["axis"]],
                gene=p[ci["nearest_gene"]], set=p[ci["set"]],
                hcc1395_ari=_f(p[ci["ari"]]), hcc1395_delta=_f(p[ci["delta"]]),
                hcc1395_germ_beta=_f(p[ci["germ_beta"]]), n_cpg=_i(p[ci["n_cpg"]]),
                cn_class=p[ci["cn_class"]] if ci.get("cn_class") is not None and len(p) > ci["cn_class"] else "",
            ))
    # BRCA2 flagship (nearest gene ZAR1L; same position) -- not in credible list
    if not any(l["chrom"] == "chr13" and l["pos"] == 32315128 for l in out):
        out.append(dict(chrom="chr13", pos=32315128, axis="HP1_vs_HP1-1",
                        gene="BRCA2/ZAR1L", set="flagship_brca2",
                        hcc1395_ari=0.790, hcc1395_delta=-0.122,
                        hcc1395_germ_beta=0.211, n_cpg=None, cn_class=""))
    return out


def _f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def _i(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return None


def read_hp(r):
    """HP:Z string tag value or None."""
    try:
        return str(r.get_tag("HP"))
    except KeyError:
        return None


def read_5mc_5hmc_beta(r, s, e):
    """mean 5mC prob and mean 5hmC prob over CpGs in [s,e] for this read.
    Returns (beta_5mC, beta_5hmC, n_cpg_m) or (None,None,0)."""
    try:
        mod = r.modified_bases
    except Exception:
        return None, None, 0
    if not mod:
        return None, None, 0
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False)
          if a is not None and b is not None}
    m_vals, h_vals = [], []
    for k, calls in mod.items():
        code = k[2]
        if code not in ('m', 'h'):
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is None or not (s <= rf <= e):
                continue
            (m_vals if code == 'm' else h_vals).append(ml / 255.0)
    bm = float(np.mean(m_vals)) if m_vals else None
    bh = float(np.mean(h_vals)) if h_vals else None
    return bm, bh, len(m_vals)


def somatic_at(vcf, chrom, pos, tol=SOM_TOL):
    """Is there a PASS somatic record within +/-tol of pos? returns (bool, info_str)."""
    try:
        vf = pysam.VariantFile(vcf)
    except Exception:
        return False, "vcf_open_fail"
    try:
        for rec in vf.fetch(chrom, max(0, pos - 1 - tol), pos + tol):
            filt = list(rec.filter.keys())
            if (not filt) or ("PASS" in filt):
                af = None
                try:
                    s = rec.samples[0]
                    af = s.get("AF")
                    if isinstance(af, (list, tuple)):
                        af = af[0]
                except Exception:
                    pass
                return True, f"{rec.chrom}:{rec.pos}{rec.ref}>{','.join(rec.alts or [])} AF={af}"
    except Exception as ex:
        return False, f"fetch_err:{ex}"
    return False, ""


def grp_stats(betas):
    """list of per-read 5mC betas -> dict."""
    a = np.array([b for b in betas if b is not None])
    if len(a) == 0:
        return dict(n=0, beta=None, frac_meth=None)
    return dict(n=int(len(a)), beta=round(float(a.mean()), 4),
                frac_meth=round(float((a >= MOD_THR).mean()), 4))


def analyze(sample):
    bam_path, vcf_path = paths_for(sample)
    bam = pysam.AlignmentFile(bam_path, "rb")
    positions = load_positions()
    results = []
    for L in positions:
        chrom, pos = L["chrom"], L["pos"]
        var0 = pos - 1
        s, e = max(0, var0 - WINDOW), var0 + WINDOW
        # collect per-read betas grouped by HP tag
        hp_m = {}   # hp_tag -> list of per-read 5mC beta
        hp_h = {}   # hp_tag -> list of per-read 5hmC beta
        n_reads = 0
        n_untagged = 0
        subhap_tags = set()
        for r in bam.fetch(chrom, s, e):
            if r.flag & 0x900 or r.flag & 0x400:   # secondary/suppl/dup
                continue
            n_reads += 1
            hp = read_hp(r)
            bm, bh, ncm = read_5mc_5hmc_beta(r, s, e)
            if hp is None:
                n_untagged += 1
                hp = "untagged"
            if "-" in hp and hp != "untagged":
                subhap_tags.add(hp)
            hp_m.setdefault(hp, []).append(bm)
            hp_h.setdefault(hp, []).append(bh)

        groups = {}
        for tag in sorted(set(list(hp_m) + list(hp_h))):
            gm = grp_stats(hp_m.get(tag, []))
            # 5hmC beta (only over reads with h calls)
            hh = np.array([b for b in hp_h.get(tag, []) if b is not None])
            gm["beta_5hmC"] = round(float(hh.mean()), 4) if len(hh) else None
            groups[tag] = gm

        # somatic status in THIS sample
        som, som_info = somatic_at(vcf_path, chrom, pos)
        has_subhap = len(subhap_tags) > 0

        # axis delta (somatic-controlled): for axis 'HPk_vs_HPk-1'
        axis = L["axis"]
        hp_axis_delta = None
        hp_axis_info = None
        if "_vs_" in axis:
            main_tag, sub_tag = axis.split("_vs_")        # 'HP1','HP1-1'
            main = main_tag.replace("HP", "")              # '1'
            sub = sub_tag.replace("HP", "")                # '1-1'
            gm = groups.get(main)
            gs_ = groups.get(sub)
            if gm and gs_ and gm["n"] >= MIN_N and gs_["n"] >= MIN_N \
                    and gm["beta"] is not None and gs_["beta"] is not None:
                hp_axis_delta = round(gs_["beta"] - gm["beta"], 4)
                hp_axis_info = f"beta({sub})={gs_['beta']} - beta({main})={gm['beta']}"
            else:
                hp_axis_info = "subhap_absent_or_lowN"

        # germline delta HP1 vs HP2 (imprinting / baseline allelic)
        germ_delta = None
        germ_info = None
        g1, g2 = groups.get("1"), groups.get("2")
        if g1 and g2 and g1["n"] >= MIN_N and g2["n"] >= MIN_N \
                and g1["beta"] is not None and g2["beta"] is not None:
            germ_delta = round(g2["beta"] - g1["beta"], 4)
            germ_info = f"beta(HP2)={g2['beta']} - beta(HP1)={g1['beta']}"

        results.append(dict(
            chrom=chrom, pos=pos, gene=L["gene"], set=L["set"], axis=axis,
            cn_class=L.get("cn_class"),
            hcc1395_delta=L["hcc1395_delta"], hcc1395_ari=L["hcc1395_ari"],
            n_reads_window=n_reads, n_untagged=n_untagged,
            somatic_in_sample=bool(som), somatic_info=som_info,
            has_somatic_subhap=has_subhap, subhap_tags=sorted(subhap_tags),
            hp_groups=groups,
            hp_axis_delta=hp_axis_delta, hp_axis_info=hp_axis_info,
            germline_delta=germ_delta, germline_info=germ_info,
        ))
    bam.close()
    return positions, results


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in RUNDIR:
        print(f"usage: {sys.argv[0]} <{'|'.join(RUNDIR)}>")
        sys.exit(2)
    sample = sys.argv[1]
    os.makedirs(OUTDIR, exist_ok=True)
    bam_path, vcf_path = paths_for(sample)
    print(f"[57] sample={sample} cancer={CANCER[sample]}")
    print(f"     bam={bam_path}")
    print(f"     vcf={vcf_path}")
    positions, results = analyze(sample)

    # quick aggregates
    n_pos = len(results)
    n_som = sum(1 for r in results if r["somatic_in_sample"])
    n_subhap = sum(1 for r in results if r["has_somatic_subhap"])
    n_hp_axis = sum(1 for r in results if r["hp_axis_delta"] is not None)
    n_germ = sum(1 for r in results if r["germline_delta"] is not None)
    n_covered = sum(1 for r in results if r["n_reads_window"] >= 10)

    out = dict(
        meta=dict(
            script="57_cross_sample_keypos.py", sample=sample,
            cancer_type=CANCER[sample], task_type="A pilot (targeted cross-sample)",
            bam=bam_path, somatic_vcf=vcf_path, ref=REF,
            window=WINDOW, min_n=MIN_N, som_tol=SOM_TOL, mod_thr=MOD_THR,
            n_key_positions=n_pos,
            method="validated pysam (54 crossval Pearson=1.0); 5mC primary, 5hmC separate; HP:Z split",
        ),
        aggregates=dict(
            n_positions=n_pos, n_somatic_in_sample=n_som,
            n_with_somatic_subhap=n_subhap, n_hp_axis_computable=n_hp_axis,
            n_germline_delta_computable=n_germ, n_covered_ge10reads=n_covered,
        ),
        per_position=results,
    )
    outp = f"{OUTDIR}/{sample}_keypos.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
    print(f"[57] wrote {outp}")
    print(f"[57] n_pos={n_pos} somatic_in_sample={n_som} with_subhap={n_subhap} "
          f"hp_axis_computable={n_hp_axis} germline_computable={n_germ} covered>=10={n_covered}")


if __name__ == "__main__":
    main()
