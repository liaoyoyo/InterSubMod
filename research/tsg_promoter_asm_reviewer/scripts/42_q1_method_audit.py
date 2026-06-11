#!/usr/bin/env python3
"""
42_q1_method_audit.py — Q1 method re-audit (gate) for CN-confound pilot.

PART A (O1 beta recompute):
  The original genome-wide dual-axis ASM survey (scripts/18_dual_axis_pivot.py,
  driven by 19_full_tp_asm_batch.sh) DELETES the per-chr Level1 raw methylation
  TSV right after pivoting (19_full_tp_asm_batch.sh:63 `rm -f "$L1"`,
  comment line 7 "Disk-safe: per-chr Level1 deleted after pivot"). So the raw
  Level1 files under genome_survey_v2/msa_tmp/<chr>_tp/ no longer exist
  (only level2/level3 summaries remain).

  Because the pipeline is fully deterministic (MSA binary + REF + tumor/normal
  BAMs + LOH bed all still present, and --window extracts each variant in an
  independent +/-1000bp window so a subset VCF reproduces identical Level1 rows
  per variant), we REGENERATE Level1 for ~20 picked somatic_pos by re-running
  MSA on a small VCF (exact same flags as 19_full_tp_asm_batch.sh:52-55), then
  re-implement the beta computation from 18_dual_axis_pivot.py and compare the
  recomputed mean_delta (HP1_vs_HP1-1 axis) against the value stored in
  asm_dualaxis_tp.tsv. Tolerance 0.02.

PART B (O2 ruler):
  Read genome_survey_v2/cn_confound/ruler_validation.json (Stage B) and confirm
  ruler_ok == (pc1_ari > 0.5 AND nc1_ari < 0.15).

Output:
  genome_survey_v2/cn_confound/q1_audit.json  (full detail)
  Returns a compact summary (printed JSON).
"""
import sys, os, gzip, json, subprocess, shutil
from collections import defaultdict
import numpy as np
from scipy import stats

# ---------------------------------------------------------------- paths
PROJ = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
WORK = f"{PROJ}/genome_survey_v2"
CN_DIR = f"{WORK}/cn_confound"
TP_TSV = f"{WORK}/asm_dualaxis_tp.tsv"
RULER_JSON = f"{CN_DIR}/ruler_validation.json"
OUT_JSON = f"{CN_DIR}/q1_audit.json"

MSA = "/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/build/bin/msa"
REF = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
TUMOR = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
NORMAL = "/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam"
LOH_BED = f"{PROJ}/output/seqc2_loh_only.bed"
CLAIRS_HDR = f"{WORK}/clairs_header.txt"
SEQC2_TP = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"

AUDIT_DIR = f"{CN_DIR}/q1_audit_msa"   # scratch for regenerated MSA Level1
AUDIT_CHR = "chr19"                    # single chrom -> fast MSA
N_TARGET = 20
TOL = 0.02

# replicate 18_dual_axis_pivot.py constants
METH_THR = 0.5
MIN_N = 3
MIN_PAIRED = 5

os.environ["TMPDIR"] = "/big7_disk"


# ---------------------------------------------------------------- pivot logic (verbatim from 18_dual_axis_pivot.py)
def load_loh(loh_bed):
    loh = defaultdict(list)
    if os.path.exists(loh_bed):
        with open(loh_bed) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    try:
                        loh[p[0]].append((int(p[1]), int(p[2])))
                    except ValueError:
                        continue
    return loh


def in_loh(loh, chrom, pos):
    pos = int(pos)
    for s, e in loh.get(chrom, []):
        if s <= pos <= e:
            return "LOH"
    return "nonLOH"


def beta_of(cell):
    n = len(cell)
    if n == 0:
        return None, 0
    m = sum(1 for v in cell.values() if v >= METH_THR)
    return m / n, n


def paired_test(cpgs, g_germ, g_som):
    pg, ps = [], []
    for mpos, hpd in cpgs.items():
        bg, ng = beta_of(hpd.get(g_germ, {}))
        bs, ns = beta_of(hpd.get(g_som, {}))
        if bg is not None and bs is not None and ng >= MIN_N and ns >= MIN_N:
            pg.append(bg)
            ps.append(bs)
    n = len(pg)
    if n < MIN_PAIRED:
        return None
    deltas = [s - gv for gv, s in zip(pg, ps)]
    try:
        w_p = float(stats.wilcoxon(pg, ps).pvalue)
    except Exception:
        w_p = float("nan")
    return {
        "n_paired_cpg": n,
        "mean_germ_beta": round(float(np.mean(pg)), 4),
        "mean_som_beta": round(float(np.mean(ps)), 4),
        "mean_delta": round(float(np.mean(deltas)), 4),
        "max_abs_delta": round(float(np.max(np.abs(deltas))), 4),
        "wilcoxon_p": w_p,
    }


def pivot_level1(level1, loh):
    """Re-implement 18_dual_axis_pivot.py data load + pivot. Returns
    dict[(chrom,spos)] -> {axis: result}."""
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    opener = gzip.open if level1.endswith(".gz") else open
    with opener(level1, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        i_spos = idx["somatic_pos"]; i_mpos = idx["methyl_pos"]
        i_hp = idx["haplotype_tag"]; i_meth = idx["meth_call"]
        i_bam = idx["bam_source_id"]; i_chrom = idx["chrom"]
        i_allele = idx["somatic_allele_type"]; i_read = idx["read_id"]
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(i_meth, i_allele, i_read):
                continue
            if p[i_bam] != "tumor":
                continue
            try:
                mc = float(p[i_meth])
            except ValueError:
                continue
            key = (p[i_chrom], p[i_spos]); mpos = p[i_mpos]; rid = p[i_read]
            cell = data[key][mpos][p[i_hp]]
            if rid not in cell or mc > cell[rid]:
                cell[rid] = mc
            at = p[i_allele].strip().lower()
            if at in ("alt", "ref"):
                c2 = data[key][mpos]["ALT" if at == "alt" else "REF"]
                if rid not in c2 or mc > c2[rid]:
                    c2[rid] = mc

    out = {}
    for (chrom, spos), cpgs in data.items():
        per_axis = {}
        for (g_germ, g_som, axis) in [("1", "1-1", "HP1_vs_HP1-1"),
                                      ("2", "2-1", "HP2_vs_HP2-1")]:
            r = paired_test(cpgs, g_germ, g_som)
            if r:
                per_axis[axis] = r
        r = paired_test(cpgs, "REF", "ALT")
        if r:
            per_axis["ALT_vs_REF"] = r
        out[(chrom, spos)] = per_axis
    return out


# ---------------------------------------------------------------- PART A
def part_a():
    res = {"status": None, "reason": None, "n_checked": 0,
           "max_abs_diff": None, "reproduces": None, "per_pos": []}

    # 1) confirm raw Level1 really gone (provenance note)
    sample_dir = f"{WORK}/msa_tmp/{AUDIT_CHR}_tp"
    l1_existing = []
    if os.path.isdir(sample_dir):
        for root, _, files in os.walk(sample_dir):
            for fn in files:
                if "level1" in fn.lower():
                    l1_existing.append(os.path.join(root, fn))
    res["raw_level1_present"] = l1_existing  # expected [] -> regenerate

    # 2) pick N_TARGET HP1_vs_HP1-1 positions on AUDIT_CHR from the stored TP TSV
    stored = {}
    with open(TP_TSV) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ci["chrom"]] != AUDIT_CHR:
                continue
            if p[ci["axis"]] != "HP1_vs_HP1-1":
                continue
            spos = p[ci["somatic_pos"]]
            stored[spos] = {
                "mean_delta": float(p[ci["mean_delta"]]),
                "n_paired_cpg": int(p[ci["n_paired_cpg"]]),
                "mean_germ_beta": float(p[ci["mean_germ_beta"]]),
                "mean_som_beta": float(p[ci["mean_som_beta"]]),
                "loh_status": p[ci["loh_status"]],
            }
    # deterministic pick: sort by pos, take first N_TARGET with n_paired_cpg>=MIN_PAIRED
    picks = sorted((int(s) for s, v in stored.items()
                    if v["n_paired_cpg"] >= MIN_PAIRED))[:N_TARGET]
    if not picks:
        res["status"] = "no_positions"
        res["reason"] = f"no HP1_vs_HP1-1 positions on {AUDIT_CHR} in {TP_TSV}"
        return res
    res["picked_positions"] = picks

    # 3a) fetch REAL REF/ALT alleles for picked positions from the SEQC2 TP VCF.
    #     MSA classifies reads into haplotype groups using somatic_allele_type
    #     (ALT vs REF base at POS) AND the HP BAM tag, so the somatic-reconstructed
    #     group ("1-1") read membership depends on the correct ALT allele. Using
    #     placeholder alleles mis-partitions reads -> HP1_vs_HP1-1 axis drops below
    #     MIN_PAIRED. We therefore replicate the exact alleles 19_full_tp_asm_batch.sh
    #     pulled from SEQC2 (single-base SNV REF/ALT, same awk filter).
    pick_set = set(picks)
    alleles = {}
    q = subprocess.run(
        ["bash", "-c",
         f"bcftools view -f PASS -r {AUDIT_CHR} {SEQC2_TP} 2>/dev/null | "
         f"bcftools query -f '%POS\\t%REF\\t%ALT\\n' 2>/dev/null | "
         f"awk 'length($2)==1 && length($3)==1'"],
        capture_output=True, text=True)
    for ln in q.stdout.splitlines():
        f3 = ln.split("\t")
        if len(f3) == 3 and int(f3[0]) in pick_set:
            alleles[int(f3[0])] = (f3[1], f3[2])
    res["alleles_resolved"] = len(alleles)
    missing = [p for p in picks if p not in alleles]
    if missing:
        # keep going but record which positions lack resolved alleles
        res["alleles_missing"] = missing

    # 3b) build subset VCF exactly like 19_full_tp_asm_batch.sh:37-48 (real alleles)
    if os.path.isdir(AUDIT_DIR):
        shutil.rmtree(AUDIT_DIR)
    os.makedirs(AUDIT_DIR, exist_ok=True)
    vcf = f"{AUDIT_DIR}/{AUDIT_CHR}_tp_audit.vcf"
    with open(CLAIRS_HDR) as fh:
        hdr_lines = fh.read()
    with open(vcf, "w") as out:
        out.write(hdr_lines)
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for pos in picks:
            ref, alt = alleles.get(pos, ("A", "T"))  # fallback only if unresolved
            out.write(f"{AUDIT_CHR}\t{pos}\t.\t{ref}\t{alt}\t30\tPASS\tSOMATIC\t"
                      f"GT:GQ:DP:AF:AD\t0/1:30:50:0.2:25,10\n")
    subprocess.run(["bgzip", "-f", vcf], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", vcf + ".gz"], check=True)
    res["audit_vcf"] = vcf + ".gz"

    # 4) run MSA exactly like 19_full_tp_asm_batch.sh:52-55
    out_msa = f"{AUDIT_DIR}/{AUDIT_CHR}_tp"
    os.makedirs(out_msa, exist_ok=True)
    log = f"{AUDIT_DIR}/{AUDIT_CHR}_tp_msa.log"
    cmd = [MSA, "-v", vcf + ".gz", "-r", REF, "-t", TUMOR, "-n", NORMAL,
           "--window", "1000", "--meth-high", "0.8", "--meth-low", "0.2",
           "--min-strand-reads", "1", "--max-read-depth", "500",
           "--threads", "16", "-o", out_msa]
    with open(log, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    res["msa_rc"] = rc
    res["msa_cmd"] = " ".join(cmd)
    if rc != 0:
        res["status"] = "msa_failed"
        res["reason"] = f"MSA exit {rc}; see {log}"
        return res

    # locate regenerated Level1
    l1 = None
    for root, _, files in os.walk(out_msa):
        for fn in files:
            if fn == "level1_raw_methylation_details.tsv.gz":
                l1 = os.path.join(root, fn)
                break
    if not l1:
        res["status"] = "no_level1_regenerated"
        res["reason"] = f"MSA produced no level1 under {out_msa}"
        return res
    res["regenerated_level1"] = l1

    # 5) re-pivot + compare HP1_vs_HP1-1 mean_delta
    loh = load_loh(LOH_BED)
    pv = pivot_level1(l1, loh)

    diffs = []
    per_pos = []
    for pos in picks:
        spos = str(pos)
        st = stored.get(spos)
        recomp = pv.get((AUDIT_CHR, spos), {}).get("HP1_vs_HP1-1")
        rec = {"somatic_pos": pos,
               "stored_mean_delta": None if st is None else st["mean_delta"],
               "recomputed_mean_delta": None,
               "abs_diff": None,
               "stored_n_paired_cpg": None if st is None else st["n_paired_cpg"],
               "recomputed_n_paired_cpg": None}
        if recomp is not None:
            rec["recomputed_mean_delta"] = recomp["mean_delta"]
            rec["recomputed_n_paired_cpg"] = recomp["n_paired_cpg"]
            if st is not None:
                d = abs(recomp["mean_delta"] - st["mean_delta"])
                rec["abs_diff"] = round(d, 4)
                diffs.append(d)
        per_pos.append(rec)

    res["per_pos"] = per_pos
    res["n_checked"] = len(diffs)
    if diffs:
        res["max_abs_diff"] = round(float(max(diffs)), 4)
        res["reproduces"] = bool(max(diffs) <= TOL)
        res["status"] = "ok"
    else:
        res["status"] = "no_overlap"
        res["reason"] = ("MSA regenerated Level1 but no HP1_vs_HP1-1 axis passed "
                         "MIN_PAIRED for any picked pos (cannot compare)")
    return res


# ---------------------------------------------------------------- PART B
def part_b():
    if not os.path.exists(RULER_JSON):
        return {"status": "missing", "ruler_ok": False,
                "reason": f"{RULER_JSON} not found"}
    with open(RULER_JSON) as f:
        rj = json.load(f)
    pc1 = rj.get("pc1_ari")
    nc1 = rj.get("nc1_ari")
    pc1_min = rj.get("thresholds", {}).get("pc1_min", 0.5)
    nc1_max = rj.get("thresholds", {}).get("nc1_max", 0.15)
    recomputed_ok = bool(pc1 is not None and nc1 is not None
                         and pc1 > pc1_min and nc1 < nc1_max)
    return {
        "status": "ok",
        "pc1_ari": pc1,
        "nc1_ari": nc1,
        "pc1_min": pc1_min,
        "nc1_max": nc1_max,
        "ruler_ok_file": bool(rj.get("ruler_ok")),
        "ruler_ok_recomputed": recomputed_ok,
        "ruler_ok": recomputed_ok and bool(rj.get("ruler_ok")) == recomputed_ok,
    }


def main():
    a = part_a()
    b = part_b()

    o1_repro = a.get("reproduces")  # bool or None
    o1_maxdiff = a.get("max_abs_diff")
    n_checked = a.get("n_checked", 0)
    o2_ok = bool(b.get("ruler_ok"))

    # overall pass: ruler must pass AND (beta reproduces OR genuinely unrecoverable
    # is NOT a pass). For a gate, we require o2 ok AND o1 reproduces True.
    overall = bool(o2_ok and (o1_repro is True))

    notes_parts = []
    notes_parts.append(
        "PART A: raw Level1 under msa_tmp/<chr>_tp/ was deleted by "
        "19_full_tp_asm_batch.sh:63 (disk-safe rm). Regenerated via MSA on a "
        f"{N_TARGET}-pos subset VCF ({AUDIT_CHR}, HP1_vs_HP1-1 axis), re-pivoted "
        "with 18_dual_axis_pivot.py logic.")
    if a["status"] == "ok":
        notes_parts.append(
            f"checked {n_checked}/{N_TARGET} HP1 positions; "
            f"max|Δ(mean_delta)|={o1_maxdiff} (tol {TOL}) -> "
            f"reproduces={o1_repro}.")
    else:
        notes_parts.append(f"PART A status={a['status']}: {a.get('reason')}")
    notes_parts.append(
        f"PART B: PC1 ARI={b.get('pc1_ari')} (>{b.get('pc1_min')}), "
        f"NC1 ARI={b.get('nc1_ari')} (<{b.get('nc1_max')}) -> ruler_ok={o2_ok}.")

    full = {
        "part_a_o1_beta": a,
        "part_b_o2_ruler": b,
        "summary": {
            "o1_beta_reproduces": o1_repro,
            "o1_max_abs_diff": o1_maxdiff,
            "n_checked": n_checked,
            "o2_ruler_ok": o2_ok,
            "overall_pass": overall,
            "notes": " ".join(notes_parts),
        },
    }
    os.makedirs(CN_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(full, f, indent=2, default=str)

    print("=== Q1 AUDIT SUMMARY ===")
    print(json.dumps(full["summary"], indent=2, default=str))
    print(f"[written] {OUT_JSON}")


if __name__ == "__main__":
    main()
