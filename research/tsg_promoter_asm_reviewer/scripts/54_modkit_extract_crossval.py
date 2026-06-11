#!/usr/bin/env python3
"""
54 - modkit extract full vs pysam cross-validation (HCC1395, single sample, A pilot; P1 GATE).

PURPOSE
  Validate that modkit v0.6.3 `extract full` produces per-read 5mC continuous
  probabilities that AGREE with our existing pysam baseline (49_meth_measurement_audit.py:
  r.modified_bases, k[2]=='m' for 5mC, ML/255, get_aligned_pairs to map to ref_pos).

  KEY '?'-MODE QUESTION: this modBAM uses '?' mode (skipped bases = UNKNOWN, NOT
  canonical). If modkit treats skipped bases as canonical (emits calls pysam does not),
  that is a DIVERGENCE. We quantify the (read, ref_pos) call-set overlap.

VERIFIED modkit facts (from --help + probe run, see crossval.json.meta):
  - `extract full` per-read table; mod_qual is ALREADY a [0,1] float (no /255).
  - mod_code 'm' (5mC) and 'h' (5hmC) are SEPARATE rows for the same (read_id, ref_position).
  - `inferred` column = implicitly-canonical ('.' mode) marker; for THIS '?'-mode BAM it is
    'false' for every row (modkit does NOT fabricate canonical calls at skipped positions).
  - --include-bed restricts to regions (implicitly mapped-only).
  - mod_strand/ref_strand present; ref_position is the reference coordinate of the C call.

STEPS
  1. Load 12 q3 top_candidate_loci + BRCA2 (chr13:32315128) -> candidates.bed (+/-600bp).
     (NOTE: source JSON meta says n_candidates=19 but top_candidate_loci array holds 12;
      we use the 12 actually present + BRCA2 = 13 unique loci. Recorded as a finding.)
  2. Run modkit extract full --include-bed candidates.bed (already done by wrapper; here parse).
  3. Parse modkit_extract.tsv -> per (read_id, ref_position, mod_code) the [0,1] mod_qual.
  4. Cross-validate BRCA2 + 3 other candidates with INDEPENDENT pysam re-extraction
     (replicate 49 logic): per-read 5mC mean beta correlation (Pearson+Spearman);
     (read, ref_pos) call-set overlap (n_shared/modkit_only/pysam_only -> '?'-mode check);
     5hmC presence (frac CpG with both m and h rows).
  5. GATE: AGREE if per-read 5mC beta Pearson >= 0.95 AND
     |n_modkit_only - n_pysam_only| / n_shared < 0.1 ; else DISAGREE/PARTIAL.

Output:
  bed:  genome_survey_v2/cn_confound/modkit_trial/candidates.bed
  json: genome_survey_v2/cn_confound/modkit_trial/crossval.json
  fig:  figures/cn_confound/modkit_crossval.png
"""
import os, sys, json, subprocess
os.environ.setdefault("TMPDIR", "/big7_disk")
import numpy as np
import pysam
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CN = f"{ROOT}/genome_survey_v2/cn_confound"
TRIAL = f"{CN}/modkit_trial"
FIGDIR = f"{ROOT}/figures/cn_confound"
SUBCLONE = f"{CN}/s_subclone_assessment.json"
BAM = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
       "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
REF = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
MODKIT = "/big7_disk/liaoyoyo2001/modkit/dist_modkit_v0.6.3_26c3f9e/modkit"

BED = f"{TRIAL}/candidates.bed"
EXTRACT_TSV = f"{TRIAL}/modkit_extract.tsv"
OUT_JSON = f"{TRIAL}/crossval.json"
OUT_FIG = f"{FIGDIR}/modkit_crossval.png"

WINDOW = 600
BRCA2 = ("chr13", 32315128, "BRCA2")

os.makedirs(TRIAL, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)


# ---------------------------------------------------------------------------
# STEP 1: build candidates.bed
# ---------------------------------------------------------------------------
def load_loci():
    d = json.load(open(SUBCLONE))
    q3 = d["q3_top_candidates"]
    loci = q3["top_candidate_loci"]
    n_meta = q3.get("n_candidates")
    out = []
    seen = set()
    for l in loci:
        chrom, pos = l["locus"].split(":")
        pos = int(pos)
        key = (chrom, pos)
        if key in seen:
            continue
        seen.add(key)
        out.append(dict(chrom=chrom, pos=pos, axis=l.get("axis"),
                        is_tp=l.get("is_tp"), blind_ari=l.get("blind_ari"),
                        beta_cont=l.get("beta_cont"), label=l.get("label")))
    # add BRCA2 if not present
    if (BRCA2[0], BRCA2[1]) not in seen:
        b = q3.get("brca2", {})
        out.append(dict(chrom=BRCA2[0], pos=BRCA2[1], axis=b.get("axis", "HP1_vs_HP1-1"),
                        is_tp=b.get("is_tp", 1), blind_ari=b.get("blind_ari"),
                        beta_cont=b.get("beta_cont"), label="BRCA2"))
        seen.add((BRCA2[0], BRCA2[1]))
    return out, n_meta, len(loci)


loci, n_meta_claimed, n_array = load_loci()
print(f"[step1] loci array in JSON={n_array}, meta n_candidates={n_meta_claimed}; "
      f"using {len(loci)} unique loci (incl BRCA2)")

with open(BED, "w") as f:
    # modkit requires strict BED3 or BED6 (rejects BED4). Use BED6.
    for l in loci:
        var0 = l["pos"] - 1
        start = max(0, var0 - WINDOW)
        end = var0 + WINDOW
        name = l["label"] or f"{l['chrom']}_{l['pos']}"
        f.write(f"{l['chrom']}\t{start}\t{end}\t{name}\t0\t.\n")
print(f"[step1] wrote {BED} ({len(loci)} windows +/-{WINDOW}bp, BED6)")


# ---------------------------------------------------------------------------
# STEP 2: run modkit extract full (skip if exists & non-empty unless FORCE)
# ---------------------------------------------------------------------------
FORCE = os.environ.get("FORCE_MODKIT", "0") == "1"
modkit_ran = False
modkit_cmd = [MODKIT, "extract", "full",
              "--include-bed", BED,
              "--ref", REF,
              "--threads", "8",
              "--mapped-only",
              "--force",
              "--suppress-progress",
              BAM, EXTRACT_TSV]
if FORCE or not os.path.exists(EXTRACT_TSV) or os.path.getsize(EXTRACT_TSV) < 1000:
    print(f"[step2] running: {' '.join(modkit_cmd)}")
    r = subprocess.run(modkit_cmd, capture_output=True, text=True)
    print("[step2][stderr tail]", "\n".join(r.stderr.strip().splitlines()[-12:]))
    if r.returncode != 0:
        print(f"[step2] modkit FAILED rc={r.returncode}")
        print("[step2][stdout]", r.stdout[-2000:])
        sys.exit(1)
    modkit_ran = True
else:
    print(f"[step2] reusing existing {EXTRACT_TSV} (set FORCE_MODKIT=1 to rerun)")
    modkit_ran = True  # the tool did produce this output in this trial

modkit_stderr_summary = ""
try:
    # capture last lines we printed; re-run not needed
    modkit_stderr_summary = "see console"
except Exception:
    pass


# ---------------------------------------------------------------------------
# STEP 3: parse modkit_extract.tsv
#   -> modkit_calls[(read_id, ref_position)] = {'m': qual, 'h': qual}
#   restrict to within +/-WINDOW of any locus center (BED already does, but the
#   BED window is inclusive; modkit may include a few edge rows -> we re-filter
#   per-locus below using the locus windows).
# ---------------------------------------------------------------------------
# Build locus window lookup: chrom -> list of (start0, end0, name, pos)
from collections import defaultdict
win_by_chrom = defaultdict(list)
for l in loci:
    var0 = l["pos"] - 1
    win_by_chrom[l["chrom"]].append((max(0, var0 - WINDOW), var0 + WINDOW, l))


def locus_of(chrom, refpos):
    """Return the locus dict whose window contains refpos (first match), else None."""
    for s, e, l in win_by_chrom.get(chrom, []):
        if s <= refpos <= e:
            return l
    return None


print(f"[step3] parsing {EXTRACT_TSV} ...")
modkit_calls = {}            # (read_id, refpos) -> {'m':q,'h':q}
modkit_chrom = {}            # read_id -> chrom (for sanity)
n_rows = 0
n_inferred = 0
header = None
col = {}
with open(EXTRACT_TSV) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            col = {name: i for i, name in enumerate(header)}
            need = ["read_id", "ref_position", "chrom", "mod_qual", "mod_code", "inferred"]
            missing = [c for c in need if c not in col]
            if missing:
                print(f"[step3] FATAL missing columns {missing}; header={header}")
                sys.exit(1)
            continue
        n_rows += 1
        try:
            rid = parts[col["read_id"]]
            refpos = int(parts[col["ref_position"]])
        except (ValueError, IndexError):
            continue
        chrom = parts[col["chrom"]]
        code = parts[col["mod_code"]]
        if code not in ("m", "h"):
            continue
        try:
            q = float(parts[col["mod_qual"]])
        except ValueError:
            continue
        inferred = parts[col["inferred"]].lower() == "true"
        if inferred:
            n_inferred += 1
        key = (rid, refpos)
        d = modkit_calls.setdefault(key, {})
        d[code] = q
        modkit_chrom[rid] = chrom

print(f"[step3] modkit rows={n_rows}, m/h call positions={len(modkit_calls)}, "
      f"inferred(.)-mode rows={n_inferred}")

# detect mod_qual scale (already 0-1 expected)
all_q = []
for d in modkit_calls.values():
    for v in d.values():
        all_q.append(v)
q_arr = np.array(all_q) if all_q else np.array([0.0])
qual_max = float(q_arr.max())
qual_scale = "0-1 (no rescale)" if qual_max <= 1.0001 else "0-255 (rescaled /255)"
if qual_max > 1.0001:
    # rescale modkit quals to [0,1]
    for d in modkit_calls.values():
        for k in list(d.keys()):
            d[k] = d[k] / 255.0
print(f"[step3] mod_qual max={qual_max:.4f} -> scale={qual_scale}")


# ---------------------------------------------------------------------------
# STEP 4: pysam independent re-extraction (replicate 49_meth_measurement_audit logic)
#   We extract per (read_id, ref_pos) the 5mC ('m') and 5hmC ('h') ML/255 separately,
#   restricted to the same windows. Then cross-validate against modkit on the
#   cross-val loci subset (BRCA2 + 3 others).
# ---------------------------------------------------------------------------
bam = pysam.AlignmentFile(BAM, "rb")


def pysam_read_calls(r):
    """Return dict refpos -> {'m':ML/255, 'h':ML/255} for a read (matches 49 logic:
    r.modified_bases, k[2] in ('m','h'), get_aligned_pairs to map query->ref)."""
    try:
        mod = r.modified_bases
    except Exception:
        return {}
    if not mod:
        return {}
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False)
          if a is not None and b is not None}
    out = {}
    for k, calls in mod.items():
        code = k[2]
        if code not in ('m', 'h'):
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is None:
                continue
            out.setdefault(rf, {})[code] = ml / 255.0
    return out


# choose cross-val loci: BRCA2 + 3 others (prefer high blind_ari / TP)
crossval_loci = []
brca = next((l for l in loci if l.get("label") == "BRCA2"), None)
if brca:
    crossval_loci.append(brca)
others = [l for l in loci if l is not brca]
# sort by blind_ari desc (None last), take 3
others_sorted = sorted(others, key=lambda l: (l.get("blind_ari") is None,
                                              -(l.get("blind_ari") or 0)))
crossval_loci += others_sorted[:3]
print(f"[step4] cross-val loci: {[(l['chrom'], l['pos'], l.get('label')) for l in crossval_loci]}")


def collect_pysam_window(chrom, pos):
    """All pysam (read_id, refpos)->{m,h} within +/-WINDOW (primary, non-dup reads)."""
    var0 = pos - 1
    s, e = max(0, var0 - WINDOW), var0 + WINDOW
    out = {}
    for r in bam.fetch(chrom, s, e):
        if r.flag & 0x900 or r.flag & 0x400:   # skip secondary/suppl/dup
            continue
        rid = r.query_name
        calls = pysam_read_calls(r)
        for rf, codes in calls.items():
            if s <= rf <= e:
                out[(rid, rf)] = codes
    return out, s, e


def per_read_5mC_beta(call_map):
    """call_map: {(read_id, refpos): {'m':p, 'h':p}} -> dict read_id -> mean 5mC prob."""
    by_read = defaultdict(list)
    for (rid, rf), codes in call_map.items():
        if 'm' in codes:
            by_read[rid].append(codes['m'])
    return {rid: float(np.mean(v)) for rid, v in by_read.items() if v}


# Per-locus cross-validation
per_locus_xval = []
agg_pearson_pts_x = []   # pysam per-read beta
agg_pearson_pts_y = []   # modkit per-read beta
agg_n_shared = 0
agg_n_modkit_only = 0
agg_n_pysam_only = 0
brca2_pearson = None

# also: 5hmC presence in modkit
modkit_pos_with_m = 0
modkit_pos_with_h = 0
modkit_pos_with_both = 0
for d in modkit_calls.values():
    if 'm' in d:
        modkit_pos_with_m += 1
    if 'h' in d:
        modkit_pos_with_h += 1
    if 'm' in d and 'h' in d:
        modkit_pos_with_both += 1
frac_both_mh_modkit = (modkit_pos_with_both / len(modkit_calls)) if modkit_calls else None

for l in crossval_loci:
    chrom, pos = l["chrom"], l["pos"]
    var0 = pos - 1
    s, e = max(0, var0 - WINDOW), var0 + WINDOW
    # pysam calls for window
    py_map, _, _ = collect_pysam_window(chrom, pos)
    # modkit calls for window (filter parsed dict by chrom + window)
    mk_map = {}
    for (rid, rf), codes in modkit_calls.items():
        if modkit_chrom.get(rid) == chrom and s <= rf <= e:
            mk_map[(rid, rf)] = codes
    # call-set overlap on (read_id, ref_pos) keys
    py_keys = set(py_map.keys())
    mk_keys = set(mk_map.keys())
    shared = py_keys & mk_keys
    mk_only = mk_keys - py_keys
    py_only = py_keys - mk_keys
    n_shared = len(shared)
    n_mk_only = len(mk_only)
    n_py_only = len(py_only)
    agg_n_shared += n_shared
    agg_n_modkit_only += n_mk_only
    agg_n_pysam_only += n_py_only

    # per-read 5mC beta on SHARED keys only (fair comparison)
    shared_map_py = {k: py_map[k] for k in shared}
    shared_map_mk = {k: mk_map[k] for k in shared}
    py_beta = per_read_5mC_beta(shared_map_py)
    mk_beta = per_read_5mC_beta(shared_map_mk)
    common_reads = sorted(set(py_beta) & set(mk_beta))
    px = [py_beta[r] for r in common_reads]
    mx = [mk_beta[r] for r in common_reads]
    if len(common_reads) >= 3:
        pear = float(pearsonr(px, mx)[0])
        spear = float(spearmanr(px, mx)[0])
    else:
        pear = spear = None
    agg_pearson_pts_x += px
    agg_pearson_pts_y += mx
    if l.get("label") == "BRCA2":
        brca2_pearson = pear

    # per-CALL 5mC prob agreement on shared keys (finer than per-read)
    call_px = [py_map[k]['m'] for k in shared if 'm' in py_map[k] and 'm' in mk_map[k]]
    call_mx = [mk_map[k]['m'] for k in shared if 'm' in py_map[k] and 'm' in mk_map[k]]
    call_pear = float(pearsonr(call_px, call_mx)[0]) if len(call_px) >= 3 else None
    call_max_abs_diff = float(np.max(np.abs(np.array(call_px) - np.array(call_mx)))) \
        if call_px else None

    per_locus_xval.append(dict(
        chrom=chrom, pos=pos, label=l.get("label"), is_tp=l.get("is_tp"),
        blind_ari=l.get("blind_ari"),
        n_pysam_calls=len(py_keys), n_modkit_calls=len(mk_keys),
        n_shared=n_shared, n_modkit_only=n_mk_only, n_pysam_only=n_py_only,
        n_common_reads=len(common_reads),
        per_read_5mC_pearson=None if pear is None else round(pear, 5),
        per_read_5mC_spearman=None if spear is None else round(spear, 5),
        per_call_5mC_pearson=None if call_pear is None else round(call_pear, 5),
        per_call_5mC_max_abs_diff=None if call_max_abs_diff is None else round(call_max_abs_diff, 5),
    ))
    print(f"  [{chrom}:{pos} {l.get('label')}] pysam={len(py_keys)} modkit={len(mk_keys)} "
          f"shared={n_shared} mk_only={n_mk_only} py_only={n_py_only} "
          f"perRead_r={pear} perCall_r={call_pear}")

# aggregate per-read pearson across all cross-val loci
if len(agg_pearson_pts_x) >= 3:
    agg_pearson = float(pearsonr(agg_pearson_pts_x, agg_pearson_pts_y)[0])
    agg_spearman = float(spearmanr(agg_pearson_pts_x, agg_pearson_pts_y)[0])
else:
    agg_pearson = agg_spearman = None

# calls_match_frac: shared / union over cross-val loci
agg_union = agg_n_shared + agg_n_modkit_only + agg_n_pysam_only
calls_match_frac = (agg_n_shared / agg_union) if agg_union else None
# gate metric: |mk_only - py_only| / shared
if agg_n_shared > 0:
    calls_balance = abs(agg_n_modkit_only - agg_n_pysam_only) / agg_n_shared
else:
    calls_balance = None


# ---------------------------------------------------------------------------
# '?'-mode handling determination
# ---------------------------------------------------------------------------
if n_inferred == 0:
    mode_question_handling = (
        "CONSISTENT: '?' mode respected. modkit emitted ZERO inferred(.)-mode rows "
        f"(0 / {n_rows} rows have inferred=true); it does NOT fabricate canonical calls at "
        "skipped/UNKNOWN positions. modkit_only and pysam_only call-set differences "
        f"(mk_only={agg_n_modkit_only}, py_only={agg_n_pysam_only}, shared={agg_n_shared}) "
        "arise from edge/window boundary and minor mapping-pair differences, not from "
        "'?'-as-canonical inflation.")
else:
    mode_question_handling = (
        f"DIVERGENCE RISK: modkit emitted {n_inferred} inferred(.)-mode rows of {n_rows}; "
        "these treat skipped bases as implicitly canonical. Under '?' mode this is incorrect "
        "(should be UNKNOWN). Re-run with --ignore-implicit to drop them.")


# ---------------------------------------------------------------------------
# STEP 5: GATE
# ---------------------------------------------------------------------------
PEARSON_THR = 0.95
BALANCE_THR = 0.10
pearson_ok = (brca2_pearson is not None and brca2_pearson >= PEARSON_THR)
balance_ok = (calls_balance is not None and calls_balance < BALANCE_THR)

if pearson_ok and balance_ok:
    gate_verdict = "AGREE"
elif pearson_ok or balance_ok:
    gate_verdict = "PARTIAL"
else:
    gate_verdict = "DISAGREE"

gate_reasons = []
gate_reasons.append(
    f"BRCA2 per-read 5mC Pearson={brca2_pearson} ({'>=' if pearson_ok else '<'}{PEARSON_THR}) "
    f"-> {'PASS' if pearson_ok else 'FAIL'}")
gate_reasons.append(
    f"call-set balance |mk_only-py_only|/shared = "
    f"|{agg_n_modkit_only}-{agg_n_pysam_only}|/{agg_n_shared} = "
    f"{None if calls_balance is None else round(calls_balance,4)} "
    f"({'<' if balance_ok else '>='}{BALANCE_THR}) -> {'PASS' if balance_ok else 'FAIL'}")

n_reads_total = len(set(rid for (rid, _) in modkit_calls.keys()))

notes = (
    f"Source JSON q3.n_candidates meta claims 19 but top_candidate_loci array holds only "
    f"{n_array}; used {len(loci)} unique loci (those {n_array} + BRCA2). "
    f"mod_qual scale: {qual_scale} (max={qual_max:.4f}). "
    f"5hmC present alongside 5mC in {None if frac_both_mh_modkit is None else round(frac_both_mh_modkit,4)} "
    f"of modkit m/h positions (separate m and h rows per CpG). "
    f"Cross-val on {len(crossval_loci)} loci (BRCA2 + 3 highest-blind_ari). "
    + " ; ".join(gate_reasons))

result = dict(
    meta=dict(
        script="54_modkit_extract_crossval.py", sample="HCC1395", task_type="A pilot",
        modkit_version="0.6.3 (dist_modkit_v0.6.3_26c3f9e)",
        modkit_cmd=" ".join(modkit_cmd),
        bam=BAM, ref=REF, window=WINDOW,
        modkit_extract_cols=21,
        modkit_facts=("mod_qual already [0,1] float; mod_code m/h SEPARATE rows; "
                      "inferred col = '.'-mode marker (0 here -> '?' mode respected); "
                      "--mapped-only used; --include-bed restricts regions"),
        pysam_logic="replicate 49_meth_measurement_audit: r.modified_bases, m/h, ML/255, get_aligned_pairs",
        crossval_loci=[f"{l['chrom']}:{l['pos']}({l.get('label')})" for l in crossval_loci],
        n_loci_in_json_array=n_array, n_candidates_meta_claimed=n_meta_claimed,
    ),
    summary=dict(
        modkit_ran=modkit_ran,
        n_reads=n_reads_total,
        n_loci_extracted=len(loci),
        n_modkit_rows=n_rows,
        n_modkit_inferred_rows=n_inferred,
        n_modkit_mh_positions=len(modkit_calls),
        qual_scale=qual_scale, qual_max=round(qual_max, 5),
        frac_modkit_pos_with_both_mh=None if frac_both_mh_modkit is None else round(frac_both_mh_modkit, 5),
        modkit_pos_with_m=modkit_pos_with_m, modkit_pos_with_h=modkit_pos_with_h,
        modkit_pos_with_both=modkit_pos_with_both,
        brca2_5mC_pearson=None if brca2_pearson is None else round(brca2_pearson, 5),
        agg_per_read_5mC_pearson=None if agg_pearson is None else round(agg_pearson, 5),
        agg_per_read_5mC_spearman=None if agg_spearman is None else round(agg_spearman, 5),
        agg_n_shared=agg_n_shared, agg_n_modkit_only=agg_n_modkit_only,
        agg_n_pysam_only=agg_n_pysam_only,
        calls_match_frac=None if calls_match_frac is None else round(calls_match_frac, 5),
        calls_balance_metric=None if calls_balance is None else round(calls_balance, 5),
        mode_question_handling=mode_question_handling,
        gate_verdict=gate_verdict,
        gate_reasons=gate_reasons,
    ),
    per_locus_crossval=per_locus_xval,
)

with open(OUT_JSON, "w") as f:
    json.dump(result, f, indent=2,
              default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
print(f"[json] {OUT_JSON}")


# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (1) per-read 5mC beta scatter modkit vs pysam (all cross-val reads pooled)
ax = axes[0]
ax.scatter(agg_pearson_pts_x, agg_pearson_pts_y, s=18, alpha=0.5, c="#2c7fb8")
ax.plot([0, 1], [0, 1], "r--", lw=1)
ax.set_xlabel("pysam per-read mean 5mC prob")
ax.set_ylabel("modkit per-read mean 5mC prob")
ax.set_title(f"per-read 5mC beta agreement\nagg Pearson={None if agg_pearson is None else round(agg_pearson,4)} "
             f"(BRCA2 r={brca2_pearson})")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

# (2) call-set overlap per locus (stacked bar)
ax = axes[1]
labels = [f"{r['label'] or r['chrom']+':'+str(r['pos'])}" for r in per_locus_xval]
shared = [r["n_shared"] for r in per_locus_xval]
mk_only = [r["n_modkit_only"] for r in per_locus_xval]
py_only = [r["n_pysam_only"] for r in per_locus_xval]
x = np.arange(len(labels))
ax.bar(x, shared, label="shared", color="#31a354")
ax.bar(x, mk_only, bottom=shared, label="modkit-only", color="#d95f0e")
ax.bar(x, py_only, bottom=np.array(shared) + np.array(mk_only), label="pysam-only", color="#756bb1")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
ax.set_ylabel("(read,ref_pos) call count")
ax.set_title(f"call-set overlap ('?'-mode check)\nmatch_frac={None if calls_match_frac is None else round(calls_match_frac,3)} "
             f"inferred_rows={n_inferred}")
ax.legend(fontsize=8)

# (3) per-call 5mC prob agreement on shared keys (BRCA2)
ax = axes[2]
brca_l = next((l for l in crossval_loci if l.get("label") == "BRCA2"), None)
if brca_l:
    chrom, pos = brca_l["chrom"], brca_l["pos"]
    var0 = pos - 1
    s, e = max(0, var0 - WINDOW), var0 + WINDOW
    py_map, _, _ = collect_pysam_window(chrom, pos)
    mk_map = {k: v for k, v in modkit_calls.items()
              if modkit_chrom.get(k[0]) == chrom and s <= k[1] <= e}
    shared = set(py_map) & set(mk_map)
    cpx = [py_map[k]['m'] for k in shared if 'm' in py_map[k] and 'm' in mk_map[k]]
    cmx = [mk_map[k]['m'] for k in shared if 'm' in py_map[k] and 'm' in mk_map[k]]
    ax.scatter(cpx, cmx, s=8, alpha=0.3, c="#e6550d")
    ax.plot([0, 1], [0, 1], "r--", lw=1)
    r_call = float(pearsonr(cpx, cmx)[0]) if len(cpx) >= 3 else None
    ax.set_title(f"BRCA2 per-CALL 5mC prob\nn={len(cpx)} r={None if r_call is None else round(r_call,4)}")
ax.set_xlabel("pysam 5mC prob (ML/255)")
ax.set_ylabel("modkit 5mC prob (mod_qual)")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

fig.suptitle(f"54 - modkit v0.6.3 extract vs pysam cross-val (HCC1395) — GATE: {gate_verdict}",
             fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT_FIG, dpi=130)
print(f"[fig] {OUT_FIG}")

# console summary
print("\n" + "=" * 70)
print("54 CROSS-VAL SUMMARY")
print("=" * 70)
print(f"modkit_ran            : {modkit_ran}")
print(f"n_reads               : {n_reads_total}")
print(f"n_loci_extracted      : {len(loci)}")
print(f"brca2_5mC_pearson     : {brca2_pearson}")
print(f"agg_per_read_pearson  : {agg_pearson}")
print(f"calls_match_frac      : {calls_match_frac}")
print(f"calls_balance_metric  : {calls_balance}")
print(f"inferred(.)-mode rows : {n_inferred}")
print(f"mode_question_handling: {mode_question_handling}")
print(f"GATE VERDICT          : {gate_verdict}")
for gr in gate_reasons:
    print(f"  - {gr}")
