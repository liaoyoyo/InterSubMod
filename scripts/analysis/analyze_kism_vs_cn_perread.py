#!/usr/bin/env python3
"""
analyze_kism_vs_cn_perread.py
ISM 甲基分群 k_ISM × CN/SV per-read 整合分析 (HCC1395)

回答 4 問 (plan: ism-k-ism-zany-tiger):
  Q1 k_ISM vs k_CN 相關
  Q2 SV 斷點當 read 標籤軸是否有用
  Q3 量化多少甲基分群受 CN 影響 (cluster split 對齊 HP/allele = CN-explained)
  Q4 per-read cluster 指派

資料來源 (皆現有, 零重跑 ISM):
  - per-region: <ISM>/filtered_snv_tp/<chr>/<chr>_<pos>/<window>/{clustering,reads}/
      clustering/significance.json  -> optimal_k, passed_gating, cluster_structure.{permanova_p,dispersion_p}
      clustering/linkage_matrix.csv -> scipy linkage (cluster_i,cluster_j,distance,new_cluster_id,size)
      reads/reads.tsv               -> read_id(0..n-1), read_name(UUID), hp, alt_support, is_tumor, strand
  - CN ref (已驗證 LOH Jaccard 0.962): SAVANA cna_normalhet segmented_absolute_copy_number.tsv
  - SV: SAVANA run sv_breakpoints_read_support.tsv (UUID), sv_breakpoints.vcf (coords)

k_eff (gate 後真實群數): optimal_k floor=2 永遠>=2 -> 用 PERMANOVA gate:
  permanova_valid & permanova_p<0.05 & dispersion_p>=0.05  -> k_eff=optimal_k (結構存在)
  否則 k_eff=1 (無結構)

k_CN: minorAlleleCopyNumber>=0.5(het-retained)->2 ; <0.5(LOH)->1 ; |CN-round|>0.3 -> fractional flag

用法:
  python3 analyze_kism_vs_cn_perread.py --label tp [--sample N] [--procs K]
"""
import argparse, json, sys, os
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.cluster.hierarchy import linkage, fcluster  # fcluster on provided Z
from multiprocessing import Pool

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    UNKNOWN_CURRENT_CLASS,
    VERIFICATION_PROVENANCE_COLUMNS,
    SchemaContractError,
    extract_provenance_frame,
    read_evidence,
    select_current_view,
    select_legacy_view,
    select_loh_legacy,
)

PROVENANCE_METADATA_FIELDS = (
    "VerificationProvenanceStatus",
    "VerificationProvenanceSourceField",
    "VerificationProvenanceWarnings",
)
PROVENANCE_EXPORT_FIELDS = VERIFICATION_PROVENANCE_COLUMNS + PROVENANCE_METADATA_FIELDS

# ---------- config / paths ----------
CANON = Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2")
SAVANA_CN = Path("/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv")
SV_READSUP = Path("/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/run/HCC1395.sv_breakpoints_read_support.tsv")
SV_VCF = Path("/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/run/HCC1395.sv_breakpoints.vcf")
OUTDIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/_assets")

CRAMERSV_ALIGNED = 0.7   # >=0.7 + Fisher sig -> aligned (per project threshold)
CRAMERSV_AMBIG = 0.5
MIN_READS_TEST = 10      # min reads for a cluster x label contingency test

# ---------- globals filled per-process ----------
_CN_SEG = None       # dict chrom -> list[(start,end,copyNumber,minorCN)]
_SV_READS = None     # set of UUIDs that are tumour-supporting for >=1 SV
_ISM_TP_DIR = None

# ---------- helpers ----------
def hp_family(hp):
    s = str(hp).strip()
    if s in ("1", "1-1"): return "fam1"
    if s in ("2", "2-1"): return "fam2"
    if s == "3": return "somatic"
    return "unphased"   # 0, nan, ""

def load_cn_segments():
    seg = defaultdict(list)
    df = pd.read_csv(SAVANA_CN, sep="\t")
    for _, r in df.iterrows():
        try:
            cn = float(r["copyNumber"]); mn = float(r["minorAlleleCopyNumber"])
        except (ValueError, TypeError):
            cn, mn = np.nan, np.nan
        seg[r["chromosome"]].append((int(r["start"]), int(r["end"]), cn, mn))
    for c in seg:
        seg[c].sort()
    return dict(seg)

def cn_at(chrom, pos):
    """Return (copyNumber, minorCN, k_CN, loh, fractional) for a locus."""
    for s, e, cn, mn in _CN_SEG.get(chrom, []):
        if s <= pos <= e:
            if np.isnan(mn):
                return cn, mn, np.nan, None, None
            loh = mn < 0.5
            k_cn = 1 if loh else 2
            frac = (not np.isnan(cn)) and abs(cn - round(cn)) > 0.3
            return cn, mn, k_cn, loh, frac
    return np.nan, np.nan, np.nan, None, None

def load_sv_reads():
    s = set()
    with open(SV_READSUP) as fh:
        next(fh, None)  # header
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                for u in parts[1].split(","):
                    u = u.strip().strip('"')
                    if u: s.add(u)
    return s

def find_region_dir(chrom, pos):
    parent = _ISM_TP_DIR / f"{chrom}" / f"{chrom}_{pos}"
    if not parent.exists():
        return None
    for d in parent.iterdir():
        if d.is_dir():
            return d
    return None

def cut_tree(window, optimal_k):
    """Recompute UPGMA from the C++ BERNOULLI read x read distance matrix (matrix.csv,
    index/cols = valid read_id), cut at C++ silhouette optimal_k.
    NOTE: scipy fcluster on the C++ linkage_matrix.csv fails (incompatible cluster-id
    scheme: 'uses the same cluster more than once'); recomputing UPGMA from the SAME
    distances with scipy is the robust equivalent (same algo, same distances).
    Returns dict {read_id(int) -> cluster_label} over valid reads."""
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage as sp_linkage, fcluster as sp_fcluster
    dm = window / "distance" / "BERNOULLI" / "matrix.csv"
    if not dm.exists():
        return None
    try:
        D = pd.read_csv(dm, index_col=0)
    except Exception:
        return None
    nval = D.shape[0]
    if nval < 4 or D.shape[0] != D.shape[1]:
        return None
    arr = D.to_numpy(dtype=float)
    arr = (arr + arr.T) / 2.0
    np.fill_diagonal(arr, 0.0)
    try:
        cond = squareform(arr, checks=False)
        Z = sp_linkage(cond, method="average")  # UPGMA
        lab = sp_fcluster(Z, t=int(optimal_k), criterion="maxclust")
    except Exception:
        return None
    return {int(rid): int(c) for rid, c in zip(D.index, lab)}

def cramers_v_fisher(cluster_labels, group_labels):
    """cluster x group contingency -> (cramers_v, fisher_or_chi2_p, n, ngroups, table_shape)."""
    df = pd.DataFrame({"c": cluster_labels, "g": group_labels})
    df = df[df["g"].notna()]
    if len(df) < MIN_READS_TEST:
        return None
    ct = pd.crosstab(df["c"], df["g"])
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return None
    n = ct.values.sum()
    chi2, _, _, _ = stats.chi2_contingency(ct.values, correction=False)
    v = np.sqrt(chi2 / (n * (min(ct.shape) - 1))) if n > 0 else 0.0
    # p: fisher if 2x2 else chi2
    if ct.shape == (2, 2):
        _, p = stats.fisher_exact(ct.values)
    else:
        try:
            p = stats.chi2_contingency(ct.values)[1]
        except Exception:
            p = np.nan
    return dict(cramers_v=float(v), p=float(p), n=int(n), table_shape=ct.shape)

def adjusted_rand(a, b):
    from sklearn.metrics import adjusted_rand_score
    m = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(m) < MIN_READS_TEST or m["b"].nunique() < 2:
        return np.nan
    return float(adjusted_rand_score(m["a"], m["b"]))


def validate_summary_provenance(df, allow_unversioned_v1=False):
    """Fail closed on schema identity and retain every available source field."""
    selected = df.copy()
    if "VerificationSchemaVersion" in selected.columns:
        current = select_current_view(selected)
        select_legacy_view(selected)
        read_evidence(selected)
        select_loh_legacy(selected)
        missing = [field for field in VERIFICATION_PROVENANCE_COLUMNS if field not in selected.columns]
        if missing:
            raise SchemaContractError("kISM provenance: missing required fields: " + ", ".join(missing))
        known = current.values != UNKNOWN_CURRENT_CLASS
        if known.any():
            extract_provenance_frame(selected.loc[known])
        status = "V2"
        source = current.field
        warning_payload = {"unknown_current_counts": current.unknown_counts}
    else:
        if not allow_unversioned_v1:
            raise SchemaContractError(
                "kISM provenance: VerificationSchemaVersion missing; explicit "
                "--allow-unversioned-v1 authorization is required for historical H1 input"
            )
        provenance = extract_provenance_frame(selected, allow_unversioned_v1=True)
        status = "UNVERSIONED_V1"
        source = str(provenance.attrs.get("selection_field", "VerificationClass"))
        warning_payload = {
            "authorization": "--allow-unversioned-v1",
            "messages": list(provenance.attrs.get("warnings", [])),
        }
        for field in VERIFICATION_PROVENANCE_COLUMNS:
            if field not in selected.columns:
                selected[field] = ""
    selected["VerificationProvenanceStatus"] = status
    selected["VerificationProvenanceSourceField"] = source
    selected["VerificationProvenanceWarnings"] = json.dumps(
        warning_payload, ensure_ascii=False, sort_keys=True
    )
    return selected

# ---------- per-region worker ----------
def analyze_region(row):
    chrom, pos = row["Chr"], int(row["Pos"])
    window = find_region_dir(chrom, pos)
    if window is None:
        return None
    sig_p = window / "clustering" / "significance.json"
    reads_p = window / "reads" / "reads.tsv"
    if not sig_p.exists() or not reads_p.exists():
        return None
    try:
        sig = json.loads(sig_p.read_text())
    except Exception:
        return None
    reads = pd.read_csv(reads_p, sep="\t")
    n = len(reads)
    optimal_k = int(sig.get("optimal_k", 0) or 0)
    cs = sig.get("cluster_structure", {})
    perm_valid = bool(cs.get("permanova_valid", False))
    perm_p = float(cs.get("permanova_p", 1.0))
    disp_p = float(cs.get("dispersion_p", 1.0))
    disp_warn = bool(cs.get("dispersion_warning", False))
    structure_exists = perm_valid and (perm_p < 0.05) and (disp_p >= 0.05)
    k_eff = optimal_k if structure_exists else 1
    ls = sig.get("label_structure", {})
    allele_perm_p = float(ls.get("allele_permanova_p", 1.0))
    allele_perm_valid = bool(ls.get("allele_permanova_valid", False))
    hp_perm_p = float(ls.get("hp_permanova_p", 1.0))

    cn, mn, k_cn, loh, frac = cn_at(chrom, pos)

    reads["hp_family"] = reads["hp"].apply(hp_family)
    reads["sv_support"] = reads["read_name"].apply(lambda u: u in _SV_READS).map({True: "SV", False: "noSV"})
    n_tumor = int((reads["is_tumor"] == 1).sum())
    n_normal = int((reads["is_tumor"] == 0).sum())

    rec = dict(
        region_id=row.get("RegionID"), chr=chrom, pos=pos,
        verification_class=row.get("VerificationClass", ""),
        n_reads=n, num_cpgs=row.get("NumCpGs", np.nan),
        optimal_k=optimal_k, k_eff=k_eff, structure_exists=structure_exists,
        permanova_p=perm_p, dispersion_p=disp_p, dispersion_warning=disp_warn,
        copyNumber=cn, minorCN=mn, k_CN=k_cn, loh=loh, fractional=frac,
        n_tumor=n_tumor, n_normal=n_normal,
        allele_permanova_p=allele_perm_p, allele_permanova_valid=allele_perm_valid, hp_permanova_p=hp_perm_p,
        hp_usable=False, alt_usable=False,
    )
    rec.update({field: row.get(field, "") for field in PROVENANCE_EXPORT_FIELDS})

    # cluster labels (per-read) only if we can cut
    lab_map = None
    labels = None
    if k_eff >= 2:
        lab_map = cut_tree(window, optimal_k)
    if lab_map is not None:
        reads["cluster"] = reads["read_id"].map(lab_map)
        reads = reads[reads["cluster"].notna()].copy()  # valid (clustered) reads only
        labels = reads["cluster"].to_numpy()
        # Q3: cluster x hp_family (germline allele) ; x alt_support (somatic allele)
        hp_axis = reads["hp_family"].where(reads["hp_family"].isin(["fam1", "fam2"]))
        hp_res = cramers_v_fisher(reads["cluster"], hp_axis)
        alt_axis = reads["alt_support"].where(reads["alt_support"].isin(["ALT", "REF"]))
        alt_res = cramers_v_fisher(reads["cluster"], alt_axis)
        sv_res = cramers_v_fisher(reads["cluster"], reads["sv_support"])
        n_sv = int((reads["sv_support"] == "SV").sum())
        ari_hp = adjusted_rand(reads["cluster"], hp_axis)
        # cluster purity wrt hp_family (max fraction of one fam per cluster, weighted)
        pur = np.nan
        try:
            sub = reads[hp_axis.notna()]
            if len(sub) >= MIN_READS_TEST:
                pur = float(sub.groupby("cluster")["hp_family"].apply(
                    lambda x: x.value_counts().iloc[0] / len(x)).mean())
        except Exception:
            pass
        ari_alt = adjusted_rand(reads["cluster"], alt_axis)
        rec.update(
            n_clusters=int(pd.Series(labels).nunique()),
            cluster_hp_cramersv=(hp_res or {}).get("cramers_v", np.nan),
            cluster_hp_p=(hp_res or {}).get("p", np.nan),
            cluster_hp_n=(hp_res or {}).get("n", np.nan),
            cluster_alt_cramersv=(alt_res or {}).get("cramers_v", np.nan),
            cluster_alt_p=(alt_res or {}).get("p", np.nan),
            cluster_alt_n=(alt_res or {}).get("n", np.nan),
            cluster_sv_cramersv=(sv_res or {}).get("cramers_v", np.nan),
            cluster_sv_p=(sv_res or {}).get("p", np.nan),
            n_sv_reads=n_sv, ari_cluster_hp=ari_hp, ari_cluster_alt=ari_alt, cluster_purity_hp=pur,
            hp_usable=hp_res is not None, alt_usable=alt_res is not None,
        )
    else:
        rec.update(n_clusters=1, cluster_hp_cramersv=np.nan, cluster_hp_p=np.nan,
                   cluster_hp_n=np.nan, cluster_alt_cramersv=np.nan, cluster_alt_p=np.nan,
                   cluster_sv_cramersv=np.nan, cluster_sv_p=np.nan, n_sv_reads=int((reads["sv_support"]=="SV").sum()),
                   ari_cluster_hp=np.nan, cluster_purity_hp=np.nan)

    # alignment_class
    rec["alignment_class"] = classify_alignment(rec)

    # per-read rows (Q4) only for k_eff>=2 regions to bound size
    perread = None
    if labels is not None:
        pr = reads[["read_name", "cluster", "hp", "hp_family", "alt_support", "is_tumor", "sv_support"]].copy()
        pr.insert(0, "region_id", row.get("RegionID"))
        pr["k_eff"] = k_eff; pr["k_CN"] = k_cn; pr["alignment_class"] = rec["alignment_class"]
        for field in PROVENANCE_EXPORT_FIELDS:
            pr[field] = row.get(field, "")
        perread = pr
    return rec, perread

def _axis_aligned(v, p):
    return (not pd.isna(v)) and (v >= CRAMERSV_ALIGNED) and (not pd.isna(p)) and (p < 0.05)

def classify_alignment(rec):
    """Use best available allele axis (germline hp_family OR somatic ALT/REF).
    germline-HP is the CN/LOH-relevant allele axis; ALT/REF is the somatic allele axis.
    In somatic-SNV regions germline-HP is often unusable (reads hp=3/0), so fall back to ALT/REF."""
    if rec["k_eff"] < 2:
        return "no_structure"
    hp_v, hp_p = rec.get("cluster_hp_cramersv", np.nan), rec.get("cluster_hp_p", np.nan)
    al_v, al_p = rec.get("cluster_alt_cramersv", np.nan), rec.get("cluster_alt_p", np.nan)
    hp_usable, al_usable = not pd.isna(hp_v), not pd.isna(al_v)
    loh = rec.get("loh")
    if not hp_usable and not al_usable:
        return "no_usable_allele_label"
    hp_aligned = _axis_aligned(hp_v, hp_p)
    al_aligned = _axis_aligned(al_v, al_p)
    best_v = np.nanmax([v for v in [hp_v, al_v] if not pd.isna(v)] or [np.nan])
    if loh is True:
        # germline allele lost; a real cluster split here is candidate beyond simple CN dosage
        # (still confounded by cis-ASM / somatic ALT — flagged, not a subclone claim)
        if al_aligned:
            return "aligned_somatic_allele_in_LOH"
        return "candidate_beyond_CN"
    # het-retained or unknown CN
    if hp_aligned:
        return "CN_explained_HP"      # aligns to germline haplotype = CN/allele-dosage explainable
    if al_aligned:
        return "aligned_somatic_allele"
    if not pd.isna(best_v) and best_v >= CRAMERSV_AMBIG:
        return "ambiguous"
    return "unaligned"

# ---------- driver ----------
def init_worker(ism_dir):
    global _CN_SEG, _SV_READS, _ISM_TP_DIR
    _ISM_TP_DIR = ism_dir
    _CN_SEG = load_cn_segments()
    _SV_READS = load_sv_reads()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", choices=["tp", "fp"], default="tp")
    ap.add_argument("--sample", type=int, default=0, help="0=all; else random N regions")
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize historical H1 summaries without VerificationSchemaVersion.",
    )
    args = ap.parse_args()

    ism_dir = CANON / f"intersubmod_{args.label}" / f"filtered_snv_{args.label}"
    sig_summary = CANON / f"intersubmod_{args.label}" / "significance_summary.csv"
    sig = validate_summary_provenance(
        pd.read_csv(sig_summary, keep_default_na=False, low_memory=False),
        allow_unversioned_v1=args.allow_unversioned_v1,
    )
    print(f"[{args.label}] significance_summary: {len(sig)} regions", flush=True)
    if args.sample and len(sig) > args.sample:
        sig = sig.sample(args.sample, random_state=args.seed)
        print(f"  sampled {len(sig)} regions", flush=True)

    rows = [r for _, r in sig.iterrows()]
    init_worker(ism_dir)  # also load globals in main for single-proc fallback
    results, perreads = [], []
    if args.procs > 1:
        with Pool(args.procs, initializer=init_worker, initargs=(ism_dir,)) as pool:
            for out in pool.imap_unordered(analyze_region, rows, chunksize=50):
                if out:
                    rec, pr = out
                    results.append(rec)
                    if pr is not None: perreads.append(pr)
    else:
        for r in rows:
            out = analyze_region(r)
            if out:
                rec, pr = out
                results.append(rec);
                if pr is not None: perreads.append(pr)

    res = pd.DataFrame(results)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tag = args.label + (f"_sample{args.sample}" if args.sample else "")
    res.to_csv(OUTDIR / f"region_table_{tag}.tsv", sep="\t", index=False)
    print(f"\n[{args.label}] analyzed regions: {len(res)}", flush=True)

    if perreads:
        pr_all = pd.concat(perreads, ignore_index=True)
        pr_all.to_csv(OUTDIR / f"perread_table_{tag}.tsv.gz", sep="\t", index=False, compression="gzip")
        print(f"  per-read rows (k_eff>=2 regions): {len(pr_all)}", flush=True)

    # ---- aggregate ----
    agg = build_aggregate(res)
    (OUTDIR / f"aggregate_{tag}.json").write_text(json.dumps(agg, indent=2, default=str))
    receipt = {
        "input": str(sig_summary),
        "schema_status": str(sig["VerificationProvenanceStatus"].iloc[0]),
        "selection_field": str(sig["VerificationProvenanceSourceField"].iloc[0]),
        "allow_unversioned_v1": args.allow_unversioned_v1,
        "provenance_fields": list(PROVENANCE_EXPORT_FIELDS),
        "warnings": sorted(sig["VerificationProvenanceWarnings"].astype(str).unique().tolist()),
    }
    (OUTDIR / f"verification_schema_contract_{tag}.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print("\n===== AGGREGATE =====", flush=True)
    print(json.dumps(agg, indent=2, default=str), flush=True)

def build_aggregate(res):
    a = {}
    a["n_regions"] = int(len(res))
    a["n_with_CN"] = int(res["k_CN"].notna().sum())
    a["n_structure_exists"] = int(res["structure_exists"].sum())
    a["pct_structure_exists"] = float(res["structure_exists"].mean()) if len(res) else np.nan
    # Q1: k_eff vs k_CN
    q1 = res.dropna(subset=["k_CN"]).copy()
    q1["k_CN_i"] = q1["k_CN"].astype(int)
    a["Q1_kEff_vs_kCN_crosstab"] = pd.crosstab(q1["k_eff"], q1["k_CN_i"]).to_dict()
    try:
        rho, prho = stats.spearmanr(q1["k_eff"], q1["k_CN_i"])
        a["Q1_spearman_rho"] = float(rho); a["Q1_spearman_p"] = float(prho)
    except Exception:
        a["Q1_spearman_rho"] = None
    a["Q1_keff_median_by_LOH"] = q1.groupby("loh")["k_eff"].median().to_dict()
    # label availability among structured regions
    s = res[res["structure_exists"]]
    n_struct = len(s)
    a["Q3_structured_n"] = int(n_struct)
    if n_struct:
        a["Q3_hp_usable_pct"] = float(s["hp_usable"].mean() * 100)
        a["Q3_alt_usable_pct"] = float(s["alt_usable"].mean() * 100)
        a["Q3_alignment_class_counts"] = s["alignment_class"].value_counts().to_dict()
        a["Q3_alignment_class_pct"] = (s["alignment_class"].value_counts() / n_struct * 100).round(2).to_dict()
    # CN-influence headline: among het-retained structured regions
    het = s[s["loh"] == False]
    a["Q3_het_structured_n"] = int(len(het))
    if len(het):
        a["Q3_pct_CN_explained_HP_among_het"] = float((het["alignment_class"] == "CN_explained_HP").mean() * 100)
        a["Q3_pct_aligned_somatic_allele_among_het"] = float((het["alignment_class"] == "aligned_somatic_allele").mean() * 100)
        a["Q3_pct_unaligned_among_het"] = float((het["alignment_class"] == "unaligned").mean() * 100)
    loh = s[s["loh"] == True]
    a["Q3_loh_structured_n"] = int(len(loh))
    if len(loh):
        a["Q3_pct_candidate_beyond_CN_among_loh"] = float((loh["alignment_class"] == "candidate_beyond_CN").mean() * 100)
    # context: somatic allele-PERMANOVA significance rate (ASM signal, label-level not cluster-level)
    av = res.dropna(subset=["allele_permanova_p"])
    a["context_allele_permanova_sig_pct"] = float((av["allele_permanova_p"] < 0.05).mean() * 100) if len(av) else None
    # Q2: SV axis on SV-informative regions
    svr = s.dropna(subset=["cluster_sv_cramersv"])
    svr = svr[svr["n_sv_reads"] >= 3]
    a["Q2_sv_informative_regions"] = int(len(svr))
    if len(svr):
        a["Q2_sv_aligned_pct"] = float(((svr["cluster_sv_cramersv"] >= CRAMERSV_ALIGNED) & (svr["cluster_sv_p"] < 0.05)).mean() * 100)
        a["Q2_hp_aligned_pct_same_subset"] = float(((svr["cluster_hp_cramersv"] >= CRAMERSV_ALIGNED) & (svr["cluster_hp_p"] < 0.05)).mean() * 100)
    # Q4: per-read coverage
    a["Q4_regions_with_perread"] = int((res["k_eff"] >= 2).sum())
    a["Q4_mean_cluster_purity_hp"] = float(res["cluster_purity_hp"].dropna().mean()) if res["cluster_purity_hp"].notna().any() else None
    a["Q4_median_ari_cluster_hp"] = float(res["ari_cluster_hp"].dropna().median()) if res["ari_cluster_hp"].notna().any() else None
    return a

if __name__ == "__main__":
    main()
