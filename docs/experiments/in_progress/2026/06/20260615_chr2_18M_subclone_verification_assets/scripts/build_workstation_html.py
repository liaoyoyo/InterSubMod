#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the chr2:18M HCC1395 molecular-state candidate JUDGMENT WORKSTATION.

Design (required-field omission guard):
  - EVERY metric is pulled from independent_audit.json / concordance.tsv by key path.
  - req() raises if a REQUIRED metric is missing -> the build REFUSES rather than
    silently rendering a dash (this is exactly the bug that froze in page-04 Fig4).
  - Legitimately-absent values (e.g. DORADO n=1 -> FDR NA) render as explicit "NA".
  - The legacy internal CLEAN/CONFOUNDED/DORADO-replicate classification is
    COMPUTED from the data (normal HP1-vs-HP2 FDR<0.05 = confounded), never
    hand-assigned. Public output spells CLEAN as normal-ASM-screen-negative;
    it does not mean acquired or clone-specific methylation.

Sources (machine-deterministic; re-run independent_subclone_audit.py = byte-identical):
  independent_audit.json   -> all read-level metrics
  chr2_18M_seqc2_concordance.tsv -> SEQC2 TVAF + historical approximate transform
Narrative rulings are bounded to the current claim ceiling and retain verdict_02
file:line citations; historical lineage-shaped labels are not emitted verbatim.

Run:
  python3 build_workstation_html.py
  python3 build_workstation_html.py --html-only

Revision provenance:
  - A dirty generator or data source is always rendered as UNCOMMITTED_CANDIDATE.
  - --source-commit is accepted only for a clean checkout whose HEAD is that exact
    40-hex commit. It cannot be used to relabel dirty content as committed.
  - published_at remains null until a separate publication workflow records it.
Outputs:
  <assets>/display/workstation_data.json
  InterSubMod/docs/explain/07_subclone-judgment-workstation-chr2-18M.standalone.html
"""
import argparse, datetime, hashlib, html, json, os, re, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.dirname(HERE)                       # ...verification_assets
DEFAULT_DATA = os.path.join(ASSETS, "data")
DISPLAY = os.path.join(ASSETS, "display")
REPO = os.path.abspath(os.path.join(ASSETS, "..", "..", "..", "..", "..", ".."))  # InterSubMod
EXPLAIN = os.path.join(REPO, "docs", "explain")
HTML_OUT = os.path.join(EXPLAIN, "07_subclone-judgment-workstation-chr2-18M.standalone.html")
JSON_OUT = os.path.join(DISPLAY, "workstation_data.json")

AUDIT_REL = "docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/independent_audit.json"
VERDICT_REL = "docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_independent_verdict_02.md"
CONC_REL = "docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/chr2_18M_seqc2_concordance.tsv"
GENERATOR_REL = "docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/build_workstation_html.py"

# Preserve the identity of the original producer separately from later revisions.
ORIGINAL_PRODUCER = {
    "name": "InterSubMod Research",
    "commit": "ea26e06cac860bd66727dec55ad84646c25456c5",
    "source_created_at": "2026-06-16T19:27:44+08:00",
}
REVISION_INPUTS = [GENERATOR_REL, AUDIT_REL, VERDICT_REL, CONC_REL]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-commit",
        help="Exact 40-hex commit containing this clean generator and all source inputs.",
    )
    parser.add_argument(
        "--html-only",
        action="store_true",
        help="Regenerate only page 07; leave workstation_data.json and page-04 Fig4 untouched.",
    )
    parser.add_argument(
        "--data-dir",
        default=DEFAULT_DATA,
        help="Directory containing independent_audit.json and chr2_18M_seqc2_concordance.tsv.",
    )
    return parser.parse_args()


ARGS = parse_args()
DATA = os.path.abspath(ARGS.data_dir)

# ----------------------------------------------------------------------------- helpers
class MissingMetric(RuntimeError):
    pass

def req(value, name):
    """Refuse the build if a required, schema-listed metric is missing."""
    if value is None:
        raise MissingMetric(f"REQUIRED metric missing -> refuse to render: {name}")
    return value

def isnum(x):
    return isinstance(x, (int, float)) and not (isinstance(x, float) and x != x)

def fpct(x, nd=0):
    if not isnum(x):
        return "NA"
    return f"{round(x*100):.0f}%" if nd == 0 else f"{x*100:.{nd}f}%"

def f2(x):
    if not isnum(x):
        return "NA"
    return f"{x:.2f}"

def fq(x):
    """format a p/FDR value compactly."""
    if not isnum(x):
        return "NA"
    if x < 1e-3:
        return f"{x:.1e}"
    return f"{x:.3f}".rstrip("0").rstrip(".")

def esc(s):
    return html.escape(str(s))

def git_full_head():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO).decode().strip()
    except Exception:
        return "UNKNOWN"

def git_branch():
    try:
        return subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=REPO).decode().strip()
    except Exception:
        return "UNKNOWN"


def dirty_revision_inputs():
    """Return tracked/untracked source paths that differ from HEAD."""
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain=v1", "--", *REVISION_INPUTS],
            cwd=REPO,
        ).decode()
    except Exception:
        return ["UNKNOWN_GIT_STATUS"]
    dirty = []
    for line in out.splitlines():
        if line:
            dirty.append(line[3:])
    return dirty


def git_tracked(path):
    try:
        result = subprocess.run(
            ["git", "ls-files", "--error-unmatch", "--", path],
            cwd=REPO,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return result.returncode == 0
    except Exception:
        return False


def resolve_revision_provenance():
    base_commit = git_full_head()
    dirty_inputs = dirty_revision_inputs()
    unversioned_inputs = []
    for logical_name, physical_path in (
        (AUDIT_REL, os.path.join(DATA, "independent_audit.json")),
        (CONC_REL, os.path.join(DATA, "chr2_18M_seqc2_concordance.tsv")),
    ):
        if os.path.abspath(physical_path) != os.path.join(REPO, logical_name) or not git_tracked(logical_name):
            unversioned_inputs.append(os.path.abspath(physical_path))
    requested = ARGS.source_commit
    if requested:
        if not re.fullmatch(r"[0-9a-fA-F]{40}", requested):
            raise ValueError("--source-commit must be an exact 40-hex commit")
        requested = requested.lower()
        if dirty_inputs:
            raise RuntimeError(
                "--source-commit refused: revision inputs are dirty: " + ", ".join(dirty_inputs)
            )
        if unversioned_inputs:
            raise RuntimeError(
                "--source-commit refused: source inputs are not versioned by that commit: "
                + ", ".join(unversioned_inputs)
            )
        if base_commit != requested:
            raise RuntimeError(
                f"--source-commit refused: checkout HEAD is {base_commit}, not {requested}"
            )
        revision_state = "COMMITTED_SOURCE"
        source_commit = requested
    elif dirty_inputs or unversioned_inputs:
        revision_state = "UNCOMMITTED_CANDIDATE"
        source_commit = None
    else:
        revision_state = "COMMITTED_SOURCE"
        source_commit = base_commit
    return {
        "original_producer": ORIGINAL_PRODUCER,
        "revision_state": revision_state,
        "source_commit": source_commit,
        "base_commit": base_commit,
        "branch": git_branch(),
        "dirty_revision_inputs": dirty_inputs,
        "unversioned_source_inputs": unversioned_inputs,
        "data_dir": DATA,
        "generated_at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "last_verified_at": "2026-08-13",
        "verification_commit": None,
        "verification_receipt": "research/20260813_public_docs_p0_correction/html_browser_qa_receipt.json",
        "published_at": None,
        "science_recomputed": False,
        "claim_ceiling_note": "本輪只校正公開敘述與呈現，未重算 science；不得把文件驗證升格為科學結果重驗。",
    }


PROVENANCE = resolve_revision_provenance()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

# ----------------------------------------------------------------------------- load
with open(os.path.join(DATA, "independent_audit.json"), encoding="utf-8") as fh:
    AUD = json.load(fh)

def load_concordance():
    rows = []
    with open(os.path.join(DATA, "chr2_18M_seqc2_concordance.tsv"), encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, vals)))
    return rows
CONC = load_concordance()

INPUT_PATHS = {
    "audit_json": os.path.join(DATA, "independent_audit.json"),
    "verdict_02": os.path.join(REPO, VERDICT_REL),
    "concordance_tsv": os.path.join(DATA, "chr2_18M_seqc2_concordance.tsv"),
}
PROVENANCE["input_files"] = {
    key: {"machine_path": path, "sha256": sha256_file(path)}
    for key, path in INPUT_PATHS.items()
}

SAMPLES = AUD["samples"]
HKU_T = SAMPLES["HKU_T"]
HKU_N = SAMPLES["HKU_N"]
DOR_T = SAMPLES["DORADO_TAGGED_T"]
DOR_N = SAMPLES["DORADO_N"]
BC = AUD["benchmark_context"]
SITE = BC["site_status"]
SNV = AUD["snvs"]
CPGC = AUD["cpg_c_positions"]

# ----------------------------------------------------------------------------- extract sSNV
ROLE = {
    "1": "beta-like local marker set A",
    "2": "beta-like local marker set A",
    "3": "alpha local read-state label（HC gap）",
    "4": "beta-like local marker set B（homopolymer-uncertain）",
    "5": "alpha+(5) nested-read pattern（low VAF；非方向）",
    "6": "beta-like local marker set B",
}
CONC_BY_POS = {r["pos1"]: r for r in CONC}

def allele_str(d):
    # d = {"counts":{base:n}, "strands":..., "hp_by_base":...}
    counts = d.get("counts", d) if isinstance(d, dict) else None
    if not counts:
        return "NA"
    return ", ".join(f"{b}:{counts[b]}" for b in sorted(counts, key=lambda k: -counts[k]))

snv_cards = []
for k in ["1", "2", "3", "4", "5", "6"]:
    s = SNV[k]
    st = SITE[k]
    pos = req(s.get("pos"), f"snv[{k}].pos")
    conc = CONC_BY_POS.get(str(pos), {})
    t_bq = HKU_T["allele_counts_bq20"][k]
    d_bq = DOR_T["allele_counts_bq20"][k]
    n_bq = HKU_N["allele_counts_bq20"][k]
    dn_bq = DOR_N["allele_counts_bq20"][k]
    snv_cards.append({
        "idx": k,
        "pos": pos,
        "ref": req(s.get("ref"), f"snv[{k}].ref"),
        "alt": req(s.get("alt"), f"snv[{k}].alt"),
        "seqc2_status": req(st.get("status"), f"site[{k}].status"),
        "in_hc": st.get("in_hc"),
        "truth_records": st.get("truth_records", []),
        # Do not expose lineage-shaped historical source labels as current findings.
        "role": ROLE[k],
        "tvaf": conc.get("tvaf", "."),
        "ccf": conc.get("ccf_est", "."),
        # The source column stores a historical, case-specific TVAF transform.
        # It is not canonical CN/LOH/purity/multiplicity-corrected CCF.
        "clonal_class": ("historical_tvaf_transform_noncanonical"
                         if conc.get("ccf_est", ".") not in (".", "", "n/a") else "n/a"),
        "hku_t_bq20": allele_str(t_bq),
        "dor_t_bq20": allele_str(d_bq),
        "hku_n_bq20": allele_str(n_bq),
        "dor_n_bq20": allele_str(dn_bq),
    })

# ----------------------------------------------------------------------------- LOH / backbone
def hp(sample):
    c = sample["hp_parent_counts"]
    h2 = c.get("2", 0); h1 = c.get("1", 0); none = c.get("None", 0)
    tagged = h2 + h1
    return {"HP2": h2, "HP1": h1, "None": none, "tagged": tagged,
            "HP2_pct": (h2 / tagged) if tagged else None,
            "HP1_pct": (h1 / tagged) if tagged else None}
loh = {
    "hku": hp(HKU_T), "dorado": hp(DOR_T),
    "loh_interval": "chr2:16,146,119-22,100,000",
    "hku_anchor": HKU_T["lineage_anchor_counts"],
    "dor_anchor": DOR_T["lineage_anchor_counts"],
    "hku_allref": HKU_T["all_ref_coverage_check"],
    "dor_allref": DOR_T["all_ref_coverage_check"],
}

# ----------------------------------------------------------------------------- linkage
LINK_PAIRS = ["E1_vs_E2", "E1_vs_E3", "E2_vs_E3", "E3_vs_E5", "E3_vs_E4ANY", "E3_vs_E6",
              "E4ANY_vs_E6", "E5_vs_E6"]
LINK_LABEL = {
    "E1_vs_E2": "(1)G vs (2)C", "E1_vs_E3": "(1)G vs (3)A", "E2_vs_E3": "(2)C vs (3)A",
    "E3_vs_E5": "(3)A vs (5)C", "E3_vs_E4ANY": "(3)A vs (4)alt", "E3_vs_E6": "(3)A vs (6)G",
    "E4ANY_vs_E6": "(4)alt vs (6)G", "E5_vs_E6": "(5)C vs (6)G",
}
LINK_READ = {  # Bounded read-level interpretation of the observed 00/10/01/11 counts.
    "E1_vs_E2": "同一 read 共現", "E1_vs_E3": "本 read 集合未觀測共同 ALT", "E2_vs_E3": "本 read 集合未觀測共同 ALT",
    "E3_vs_E5": "(5)C 僅見於 (3)A-bearing shared-coverage reads（局部 nesting；非次序）",
    "E3_vs_E4ANY": "主要未觀測共同 ALT", "E3_vs_E6": "本 read 集合未觀測共同 ALT",
    "E4ANY_vs_E6": "同一 read-state 共現", "E5_vs_E6": "本 read 集合未觀測共同 ALT（非 lineage）",
}
def link_rows():
    out = []
    hl = HKU_T["linkage_bq10"]; dl = DOR_T["linkage_bq10"]
    for p in LINK_PAIRS:
        h = hl.get(p); d = dl.get(p)
        out.append({
            "pair": LINK_LABEL[p],
            "hku": (f"{h['00']}/{h['10']}/{h['01']}/{h['11']} (n={h['n']})" if h else "NA"),
            "dor": (f"{d['00']}/{d['10']}/{d['01']}/{d['11']} (n={d['n']})" if d else "NA"),
            "read": LINK_READ[p],
        })
    return out
LINKS = link_rows()

# ----------------------------------------------------------------------------- methylation 3-group classification (DATA-DRIVEN)
CONF_FDR = 0.05
def methyl_map(rows):
    return {r["site"]: r for r in rows}
hku_ab = methyl_map(HKU_T["methyl_alpha_vs_pos3_ref"]["rows"])
dor_ab = methyl_map(DOR_T["methyl_alpha_vs_pos3_ref"]["rows"])
norm_hp = methyl_map(HKU_N["methyl_normal_hp1_vs_hp2"]["rows"])

CPG_ORDER = ["2.1", "2.2", "3.1", "3.2", "3.3", "3.4", "3.5", "4.1", "5.1", "5.2"]
methyl_rows = []
for cg in CPG_ORDER:
    h = hku_ab.get(cg, {}); d = dor_ab.get(cg, {}); nz = norm_hp.get(cg, {})
    h_fdr = h.get("mann_whitney_fdr"); d_fdr = d.get("mann_whitney_fdr"); n_fdr = nz.get("mann_whitney_fdr")
    h_delta = h.get("mean_delta_a_minus_b"); d_delta = d.get("mean_delta_a_minus_b")
    confounded = isnum(n_fdr) and n_fdr < CONF_FDR
    dor_repl = (isnum(d_fdr) and d_fdr < CONF_FDR and isnum(h_delta) and isnum(d_delta)
                and ((h_delta > 0) == (d_delta > 0)))
    cls = "CONFOUNDED" if confounded else "CLEAN"
    methyl_rows.append({
        "cpg": cg, "c_pos": CPGC.get(cg),
        "hku_a": h.get("a_mean"), "hku_b": h.get("b_mean"),
        "hku_an": h.get("a_n"), "hku_bn": h.get("b_n"),
        "hku_delta": h_delta, "hku_fdr": h_fdr,
        "dor_a": d.get("a_mean"), "dor_b": d.get("b_mean"),
        "dor_an": d.get("a_n"), "dor_bn": d.get("b_n"),
        "dor_delta": d_delta, "dor_fdr": d_fdr,
        "norm_hp1": nz.get("a_mean"), "norm_hp2": nz.get("b_mean"), "norm_fdr": n_fdr,
        "class": cls, "dorado_replicates": bool(dor_repl),
        "strongest": (cls == "CLEAN" and dor_repl),
    })
clean_set = [r["cpg"] for r in methyl_rows if r["class"] == "CLEAN"]
conf_set = [r["cpg"] for r in methyl_rows if r["class"] == "CONFOUNDED"]
repl_set = [r["cpg"] for r in methyl_rows if r["dorado_replicates"]]
strong_set = [r["cpg"] for r in methyl_rows if r["strongest"]]

# ----------------------------------------------------------------------------- pos4 homopolymer
def pos4(sample):
    p = sample["pos4_subtype_summary"]
    bq = sample["allele_counts_bq20"]["4"]["counts"]
    return {"REF_C": bq.get("C", 0), "G": p.get("G", {}).get("n", 0),
            "T": p.get("T", {}).get("n", 0), "DEL": p.get("DEL", {}).get("n", 0),
            "G_strand": p.get("G", {}).get("strand_counts", {})}
POS4 = {"hku": pos4(HKU_T), "dor": pos4(DOR_T)}

# ----------------------------------------------------------------------------- narrative (verdict_02, cited)
RULINGS = [
    ("1", "此區與 SEQC2 LOH 標註相符", "支持・有界", "bounded", "SEQC2 LOH BED 涵蓋此區，且 HKU/DORADO tumor tagged reads 呈 HP2 imbalance；本頁未獨立確認 LOH 機制、allele-specific CN 或單親保留", "verdict_02:267"),
    ("2", "存在局部 read-level 分子狀態候選", "支持・有界", "bounded", "可稱 local molecular-state candidates；非 biological clone truth，亦非 lineage", "verdict_02:268"),
    ("3", "HP2 沒抽到／突變到不見", "否定", "negated", "HP2 是主體；all-REF 短 linkage reads 存在；完整 read 缺失是 coverage", "verdict_02:269"),
    ("4", "(1)(2)(3)(6) 先突變", "不成立・需重寫", "rewrite", "read-level 共現／未共現不足以建立 trunk、先後或 ancestry；僅支持互斥的 local read-state candidates", "verdict_02:270"),
    ("5", "(5) 是否呈現 (3)A-conditioned 局部 pattern", "支持 (provisional)", "provisional", "在有共同覆蓋的 reads 中，(5)C 僅與 (3)A 共現；HKU 4 / DORADO 2。這是局部 nesting pattern，不代表 parent/child 或發生次序", "verdict_02:271"),
    ("6", "(4) 多突變造成三群", "否定", "negated", "合併為 pos4-altered beta-like state；G/T/DEL 為 homopolymer-uncertain", "verdict_02:272"),
]
RETRACTIONS = [
    ("「5 群 / 5 subclone」", "(4) 的 G/T/DEL 應合併成一個 homopolymer-uncertain pos4-altered state；本頁僅保留 alpha、alpha+(5) nested-read pattern、beta-like 等操作性 local-state labels，不能換算為 clone 數", "verdict_02:49-59"),
    ("「VAF 排序突變」", "局部 VAF 受 CN/LOH/純度/multiplicity 影響；不能由 VAF 或局部共現 pattern 排 parent/child、先後或 ancestry；(4)(6) 共現亦不可排序", "verdict_02:154"),
    ("「無 ancestral REF read / HP2 突變到不見」", "HKU 有 17 條 reads 在 ≥2 已覆蓋 SNV 全為 REF；這只證明 all-REF local read patterns 存在，不代表已觀測 ancestor；任何 root 仍是未確認的 parsimony 假說", "verdict_02:128-132"),
]
FINAL_DEFENSIBLE = "在 SEQC2 LOH-annotated 區域中，read-level allele linkage 支持局部 molecular-state candidates；在 shared-coverage reads 內，(5)C 只於 (3)A-bearing reads 觀測到，屬 provisional nesting pattern，不代表 parent/child 或發生次序。甲基只支持 mutation-defined 分組後的有界關聯與跨 basecaller 技術重現，不確認 cellular subclone、linear ancestry、完整演化順序或因果。"
FINAL_PAPER = "Local read-level molecular-state candidates in a SEQC2 LOH-annotated interval with bounded, mutation-conditioned methylation association."
FINAL_NOT = "Confirmed five-subclone evolutionary reconstruction."
CLAIM_CEILING = "confirmed cellular subclone = 0；confirmed linear ancestry = 0。alpha / beta 等皆為本區域的操作性 read-state labels。"
CAVEATS = [
    ("整體 tier = L2", "單一區域案例 × 單一 biological sample × 兩個 technical datasets；cross-basecaller (HKU vs DORADO) 是同細胞株的技術/資料重現，非獨立 biological replicate（升 ★4 需 COLO829 等第二樣本）。", "verdict_02:314 / critical_audit_03"),
    ("pos3 (3) 落 SEQC2 HC 空隙", "out-of-HC truth-unevaluable，不能由此頁判為 SEQC2 FP 或 TP；alpha 操作性標籤只靠內部 read linkage，外部無真值可確認。", "verdict_02:81-91"),
    ("甲基 confound 未完全形式排除", "3.3/3.4/3.5/4.1 在 matched-normal 已有顯著 germline ASM；其餘只能稱 tumor-associated、normal-ASM-screen-negative，不能據此證明 acquired 或 cellular-subclone-specific。", "verdict_02:205-225"),
    ("「0 違反」僅 strict 定義", "broadened pos4-beta 在 DORADO 有 1 條 discordant read；多組遠距 pair 無跨距 read，故「0 observed」非強 cross-data 複製。", "critical_audit_03 High-2"),
    ("DORADO normal 未 tag", "germline-ASM confound 檢定僅靠 HKU normal；confound 側未跨 basecaller 複製。", "independent_audit.md:244-257"),
]
# A-vs-audit number conflicts footnote (critical_audit_03)
CONFLICTS = [
    ("all-REF reads (≥2 clean pts)", "報告A: 22/17", "權威(audit): 17/15", "independent_audit.md:63,152"),
    ("pos3 tumor 分子數", "報告A: 30/28 (BQ10)", "權威 BQ20: HKU 29/30・DORADO 21/27", "independent_audit.md:36,125"),
    ("DORADO 複製的乾淨 CpG", "報告A: {2.1,2.2,3.1,3.2,3.3}", "權威: 乾淨∩DORADO複製 = {3.1,3.2}（3.3 屬 confound）", "verdict_02:221 / methyl table"),
    ("「乾淨 somatic 甲基案例」", "報告A: 乾淨案例", "權威: 混合 (normal ASM + tumor-associated)；normal-ASM-screen-negative 子集={2.1,2.2,3.1,3.2,5.1,5.2}，未證明 acquired", "verdict_02:215-225"),
    ("TVAF 口徑", "早期 bwaTVAF: .396/.297/.050/.389", "權威 consensus: .403/.389/.242/.048/.382", "concordance_demo:43"),
]

# ----------------------------------------------------------------------------- data.json
BUILD = {
    # commit means the exact committed source revision, never merely the checkout base.
    "commit": PROVENANCE["source_commit"],
    "base_commit": PROVENANCE["base_commit"],
    "branch": PROVENANCE["branch"],
    "revision_state": PROVENANCE["revision_state"],
    "built_utc": PROVENANCE["generated_at"],
    "published_at": PROVENANCE["published_at"],
    "tier": "L2",
    "scope": "PARTIAL — chr2:18,066,481–18,110,828 單一 region/case（6 sSNV）× 1 biological sample × 2 technical datasets（HKU+DORADO）",
}
DATAJSON = {
    "build": BUILD,
    "provenance": PROVENANCE,
    "region": AUD["region"],
    "ccf_semantics": "historical_case_specific_approximate_tvaf_transform_noncanonical_not_clone_corroboration",
    "methyl_class_semantics": "legacy_CLEAN_means_normal_ASM_screen_negative_not_acquired_or_clone_specific",
    "role_semantics": "historical_operational_labels_not_confirmed_cellular_lineage",
    "claim_ceiling": CLAIM_CEILING,
    "sources": {"audit_json": AUDIT_REL, "verdict_02": VERDICT_REL, "concordance_tsv": CONC_REL},
    "snv_cards": snv_cards,
    "loh": loh,
    "linkage": LINKS,
    "methyl_rows": methyl_rows,
    "methyl_sets": {"clean": clean_set, "confounded": conf_set,
                    "dorado_replicate": repl_set, "strongest_clean_replicate": strong_set},
    "pos4": POS4,
    "rulings": [dict(zip(["id", "claim", "verdict", "kind", "fix", "src"], r)) for r in RULINGS],
    "retractions": [dict(zip(["was", "now", "src"], r)) for r in RETRACTIONS],
    "final_defensible": FINAL_DEFENSIBLE, "final_paper": FINAL_PAPER, "final_not": FINAL_NOT,
    "caveats": [dict(zip(["title", "body", "src"], c)) for c in CAVEATS],
    "conflicts": [dict(zip(["metric", "reportA", "authoritative", "src"], c)) for c in CONFLICTS],
}

# ----------------------------------------------------------------------------- SVG figures (data-bound)
def meth_color(v):
    """blue(0) -> grey(.5) -> red(1) for methylation level; NA -> hatched grey."""
    if not isnum(v):
        return "#3a3a44"
    # interpolate blue(40,90,170) -> light(230,230,225) -> red(200,50,55)
    if v <= 0.5:
        t = v / 0.5
        r = int(40 + (230 - 40) * t); g = int(90 + (230 - 90) * t); b = int(170 + (225 - 170) * t)
    else:
        t = (v - 0.5) / 0.5
        r = int(230 + (200 - 230) * t); g = int(230 + (50 - 230) * t); b = int(225 + (55 - 225) * t)
    return f"rgb({r},{g},{b})"

def methyl_heatmap_svg():
    """rows = 10 CpG; cols = HKU a/b, DOR a/b, normal HP1/HP2; right = class badge.
    Data-bound: every cell value comes from methyl_rows; NA rendered explicitly."""
    cols = [("HKU α", "hku_a", "hku_an"), ("HKU β", "hku_b", "hku_bn"),
            ("DOR α", "dor_a", "dor_an"), ("DOR β", "dor_b", "dor_bn"),
            ("N·HP1", "norm_hp1", None), ("N·HP2", "norm_hp2", None)]
    cw, ch = 86, 34
    x0, y0 = 150, 56
    W = x0 + cw * len(cols) + 230
    H = y0 + ch * len(methyl_rows) + 30
    out = [f'<svg viewBox="0 0 {W} {H}" width="100%" preserveAspectRatio="xMidYMid meet" role="img" aria-labelledby="methyl-heatmap-title methyl-heatmap-desc">',
           '<title id="methyl-heatmap-title">甲基化三組熱圖</title>',
           '<desc id="methyl-heatmap-desc">逐 CpG 比較 HKU 與 DORADO 的 mutation-conditioned 甲基值及 matched-normal ASM screen；screen-negative 不表示 acquired 或 clone-specific。</desc>']
    out.append(f'<rect width="{W}" height="{H}" fill="#16161d"/>')
    # col headers
    for ci, (lab, _, _) in enumerate(cols):
        cx = x0 + ci * cw + cw / 2
        sep = "#cfcfe0" if ci in (0, 2, 4) else "#9aa"
        out.append(f'<text x="{cx:.0f}" y="38" fill="{sep}" font-size="13" font-weight="600" text-anchor="middle">{esc(lab)}</text>')
    out.append(f'<text x="14" y="38" fill="#cfcfe0" font-size="13" font-weight="600">CpG</text>')
    out.append(f'<text x="{x0+cw*4-6:.0f}" y="20" fill="#7fd1c4" font-size="11" text-anchor="middle">tumor α(=pos3 A) vs β(=pos3 G)</text>')
    out.append(f'<text x="{x0+cw*5:.0f}" y="20" fill="#e0a96d" font-size="11" text-anchor="middle">normal ASM</text>')
    for ri, r in enumerate(methyl_rows):
        ry = y0 + ri * ch
        # cpg label + coord
        lblc = "#f2c94c" if r["strongest"] else ("#9aa" if r["class"] == "CONFOUNDED" else "#dfe")
        out.append(f'<text x="14" y="{ry+ch/2+4:.0f}" fill="{lblc}" font-size="12" font-weight="600">{esc(r["cpg"])}</text>')
        out.append(f'<text x="52" y="{ry+ch/2+4:.0f}" fill="#666" font-size="9">{r["c_pos"]}</text>')
        for ci, (lab, key, nkey) in enumerate(cols):
            v = r.get(key)
            cx = x0 + ci * cw
            fill = meth_color(v)
            out.append(f'<rect x="{cx}" y="{ry+2}" width="{cw-3}" height="{ch-4}" fill="{fill}" rx="3"/>')
            tcol = "#111" if (isnum(v) and 0.18 < v < 0.82) else ("#fff" if isnum(v) else "#aaa")
            txt = f2(v) if isnum(v) else "NA"
            ntxt = ""
            if nkey and isnum(r.get(nkey)):
                ntxt = f' n{r.get(nkey)}'
            out.append(f'<text x="{cx+(cw-3)/2:.0f}" y="{ry+ch/2+4:.0f}" fill="{tcol}" font-size="11" text-anchor="middle">{esc(txt)}{esc(ntxt)}</text>')
        # right badges: class + dorado replicate + FDR
        bx = x0 + cw * len(cols) + 8
        if r["class"] == "CONFOUNDED":
            out.append(f'<rect x="{bx}" y="{ry+5}" width="96" height="{ch-12}" fill="#5a2d2d" rx="4"/>')
            out.append(f'<text x="{bx+48}" y="{ry+ch/2+4:.0f}" fill="#ffb3b3" font-size="10" text-anchor="middle" font-weight="600">germline-ASM</text>')
        else:
            col = "#2e5d44" if not r["strongest"] else "#3a6e1f"
            out.append(f'<rect x="{bx}" y="{ry+5}" width="96" height="{ch-12}" fill="{col}" rx="4"/>')
            lab = "SCREEN-NEG★" if r["strongest"] else "SCREEN-NEG"
            out.append(f'<text x="{bx+48}" y="{ry+ch/2+4:.0f}" fill="#cdebcf" font-size="10" text-anchor="middle" font-weight="600">{lab}</text>')
        # dorado replicate dot
        rep = "#4ad" if r["dorado_replicates"] else "#555"
        out.append(f'<circle cx="{bx+112}" cy="{ry+ch/2:.0f}" r="6" fill="{rep}"><title>DORADO 複製 FDR&lt;0.05: {r["dorado_replicates"]}</title></circle>')
        out.append(f'<text x="{bx+124}" y="{ry+ch/2+4:.0f}" fill="#789" font-size="9">{"DOR✓" if r["dorado_replicates"] else "DOR—"}</text>')
    out.append('</svg>')
    return "".join(out)

def ccf_bar_svg():
    """Historical, case-specific TVAF transform; not canonical CCF."""
    rows = [r for r in CONC if r["ccf_est"] not in (".", "", "n/a")]
    W, H = 760, 230; x0, y0 = 60, 28; bw = 84; gap = 26
    out = [f'<svg viewBox="0 0 {W} {H}" width="100%" preserveAspectRatio="xMidYMid meet" role="img" aria-labelledby="tvaf-transform-title tvaf-transform-desc">',
           '<title id="tvaf-transform-title">歷史案例特定 TVAF 近似轉換</title>',
           '<desc id="tvaf-transform-desc">這不是 allele-specific CN、LOH、純度與 multiplicity 校正後的 canonical CCF，不能作 clone 佐證。</desc>']
    out.append(f'<rect width="{W}" height="{H}" fill="#16161d"/>')
    base = y0 + 150
    # Historical reference line retained for auditability; it is not a current gate.
    yth = base - 0.85 * 150
    out.append(f'<line x1="{x0-8}" y1="{yth:.0f}" x2="{W-20}" y2="{yth:.0f}" stroke="#f2994a" stroke-dasharray="5 4" stroke-width="1"/>')
    out.append(f'<text x="{W-22}" y="{yth-4:.0f}" fill="#f2994a" font-size="10" text-anchor="end">歷史 0.85 參考線（非現行 gate）</text>')
    order = sorted(rows, key=lambda r: float(r["ccf_est"]))
    for i, r in enumerate(order):
        ccf = float(r["ccf_est"])
        bx = x0 + i * (bw + gap)
        bh = ccf * 150
        col = "#c0504d" if ccf >= 0.85 else ("#e0a96d" if ccf >= 0.4 else "#6aa6d8")
        out.append(f'<rect x="{bx}" y="{base-bh:.0f}" width="{bw}" height="{bh:.0f}" fill="{col}" rx="3"/>')
        out.append(f'<text x="{bx+bw/2:.0f}" y="{base-bh-6:.0f}" fill="#dfe" font-size="12" font-weight="600" text-anchor="middle">{ccf:.2f}</text>')
        pos = int(r["pos1"])
        out.append(f'<text x="{bx+bw/2:.0f}" y="{base+16:.0f}" fill="#9aa" font-size="10" text-anchor="middle">chr2:{pos:,}</text>')
        out.append(f'<text x="{bx+bw/2:.0f}" y="{base+30:.0f}" fill="#789" font-size="9" text-anchor="middle">VAF {r["tvaf"]}</text>')
    out.append(f'<text x="{x0-8}" y="{base+30:.0f}" fill="#c0504d" font-size="11">↑ 近似轉換值</text>')
    out.append('</svg>')
    return "".join(out)

def tree_svg():
    """Non-temporal local read-state candidate graph; not a lineage tree."""
    return """
<svg viewBox="0 0 760 300" width="100%" preserveAspectRatio="xMidYMid meet" role="img" aria-labelledby="candidate-tree-title candidate-tree-desc">
<title id="candidate-tree-title">有界的局部 read-state 候選圖</title>
<desc id="candidate-tree-desc">無方向、非時間、非唯一的 read-level 共現圖；不代表已確認的 cellular clone、parent-child 或 lineage tree。</desc>
<rect width="760" height="300" fill="#16161d"/>
<text x="380" y="26" fill="#9aa" font-size="11" text-anchor="middle">layout 僅整理局部 read 共現／未共現；線段無方向，不編碼先後或 ancestry</text>
<rect x="270" y="42" width="220" height="36" rx="8" fill="#252530" stroke="#777"/>
<text x="380" y="64" fill="#ddd" font-size="11" text-anchor="middle">all-REF local read pattern（非 ancestor）</text>
<line x1="380" y1="78" x2="200" y2="112" stroke="#666" stroke-width="1.5" stroke-dasharray="5 5"/>
<line x1="380" y1="78" x2="560" y2="112" stroke="#666" stroke-width="1.5" stroke-dasharray="5 5"/>
<circle cx="200" cy="120" r="7" fill="#7fd1c4"/><text x="200" y="106" fill="#7fd1c4" font-size="12" text-anchor="middle" font-weight="600">alpha label: (3) G&gt;A</text>
<line x1="200" y1="127" x2="200" y2="171" stroke="#7fd1c4" stroke-width="2" stroke-dasharray="4 3"/>
<circle cx="200" cy="180" r="7" fill="#7fd1c4"/><text x="200" y="201" fill="#9fe3d6" font-size="11" text-anchor="middle">alpha+(5) nested-read pattern</text>
<text x="200" y="217" fill="#789" font-size="9" text-anchor="middle">subset pattern；非 parent/child，非時間箭頭</text>
<circle cx="560" cy="120" r="7" fill="#e0a96d"/><text x="560" y="106" fill="#e0a96d" font-size="12" text-anchor="middle" font-weight="600">beta-like operational label</text>
<line x1="560" y1="127" x2="480" y2="166" stroke="#e0a96d" stroke-width="2"/>
<line x1="560" y1="127" x2="640" y2="166" stroke="#e0a96d" stroke-width="2"/>
<text x="470" y="184" fill="#f2d6ad" font-size="11" text-anchor="middle">(1)C&gt;G + (2)G&gt;C</text>
<text x="650" y="184" fill="#f2d6ad" font-size="11" text-anchor="middle">(4) C&gt;G* + (6)C&gt;G</text>
<line x1="470" y1="194" x2="650" y2="194" stroke="#b06ad8" stroke-width="2" stroke-dasharray="6 5"/>
<text x="560" y="211" fill="#c79fe0" font-size="9" text-anchor="middle">兩 local patterns 缺直接 spanning-read bridge；關係未定</text>
<text x="650" y="228" fill="#789" font-size="9" text-anchor="middle">(4) G/T/DEL 合併為 homopolymer-uncertain state</text>
<text x="380" y="270" fill="#9aa" font-size="10" text-anchor="middle">confirmed cellular subclone = 0；confirmed linear ancestry = 0</text>
</svg>
"""

# ----------------------------------------------------------------------------- page-04 Fig4 (data-bound replacement for the v1 buggy SVG)
def fig4_page04_svg():
    """Re-render page-04 Fig4 in its 4-column (HKU a/b, DORADO a/b) teaching style,
    but FULLY data-bound from independent_audit.json. No dash-on-missing: a required
    cell that is missing would render a loud marker. Class = computed."""
    rh = 30; y0 = 56
    H = y0 + rh * len(methyl_rows) + 60
    xs = {"ha": (215, 250), "hb": (295, 330), "da": (435, 470), "db": (515, 550)}
    out = [f'<svg viewBox="0 0 880 {H}" role="img" aria-labelledby="fig4-title fig4-desc">',
           '<title id="fig4-title">甲基化翻轉熱圖資料綁定修正版</title>',
           '<desc id="fig4-desc">比較兩個 basecaller 的 mutation-conditioned 甲基值，並標示 matched-normal ASM confound；不作 acquired 或 clone-specific 結論。</desc>']
    # legend
    out.append('<g font-size="11">'
               '<rect x="600" y="14" width="16" height="16" fill="rgb(200,50,55)"/><text x="622" y="27">高甲基 (meanM→1)</text>'
               '<rect x="600" y="34" width="16" height="16" fill="rgb(230,230,225)" stroke="#E3DACC"/><text x="622" y="47">中</text>'
               '<rect x="600" y="54" width="16" height="16" fill="rgb(40,90,170)"/><text x="622" y="67">低甲基 (meanM→0)</text></g>')
    out.append('<g font-size="12" text-anchor="middle" font-weight="700">'
               '<text x="250" y="34">HKU · α</text><text x="330" y="34">HKU · β</text>'
               '<text x="470" y="34">DORADO · α</text><text x="550" y="34">DORADO · β</text></g>')
    out.append('<g font-family="monospace" font-size="11" text-anchor="middle">')
    strong_y = []
    for ri, r in enumerate(methyl_rows):
        ry = y0 + ri * rh
        lblcol = "#2F7D5B" if r["strongest"] else ("#8a8a83" if r["class"] == "CONFOUNDED" else "#1c1b19")
        out.append(f'<text x="150" y="{ry+18}" text-anchor="start" font-weight="700" fill="{lblcol}">CpG {esc(r["cpg"])}</text>')
        for key, val in [("ha", r["hku_a"]), ("hb", r["hku_b"]), ("da", r["dor_a"]), ("db", r["dor_b"])]:
            rx, tx = xs[key]
            if not isnum(val):
                # Never silently dash a required value; render an explicit small-n / NA marker.
                out.append(f'<rect x="{rx}" y="{ry+2}" width="70" height="24" fill="#2b2b30" stroke="#555"/>'
                           f'<text x="{tx}" y="{ry+18}" fill="#caa">NA</text>')
            else:
                fill = meth_color(val)
                tc = "#111" if 0.18 < val < 0.82 else "#fff"
                out.append(f'<rect x="{rx}" y="{ry+2}" width="70" height="24" fill="{fill}"/>'
                           f'<text x="{tx}" y="{ry+18}" fill="{tc}">{f2(val)}</text>')
        # annotation
        if r["class"] == "CONFOUNDED":
            out.append(f'<text x="600" y="{ry+18}" text-anchor="start" font-size="10" fill="#9A3B3B">germline-ASM confound（normal FDR {fq(r["norm_fdr"])}）</text>')
        else:
            tag = "screen-negative+DORADO複製" if r["strongest"] else ("screen-negative" + ("・DORADO複製" if r["dorado_replicates"] else "・DORADO未過FDR/小n"))
            out.append(f'<text x="600" y="{ry+18}" text-anchor="start" font-size="10" fill="#2F7D5B">{tag}</text>')
        if r["strongest"]:
            strong_y.append(ry)
    # box around strongest rows
    if strong_y:
        top = min(strong_y); bot = max(strong_y) + rh
        out.append(f'<rect x="210" y="{top}" width="380" height="{bot-top}" fill="none" stroke="#2F7D5B" stroke-width="2.5" rx="6"/>')
    out.append('</g>')
    out.append(f'<text x="400" y="{y0+rh*len(methyl_rows)+24}" text-anchor="middle" font-size="11" fill="#2F7D5B" font-weight="700">'
               f'normal-ASM-screen-negative 且跨 basecaller 技術重現 = {{{", ".join(strong_set)}}}（不證明 acquired/clone）</text>')
    out.append(f'<text x="400" y="{y0+rh*len(methyl_rows)+42}" text-anchor="middle" font-size="11" fill="#9A3B3B" font-weight="700">'
               f'被既存 germline ASM confound = {{{", ".join(conf_set)}}} → 不可當 acquired 或 clone-specific 甲基</text>')
    out.append('</svg>')
    return "".join(out)

# ----------------------------------------------------------------------------- HTML
CSS = """
:root{--bg:#0f0f14;--panel:#16161d;--panel2:#1c1c25;--ink:#e8e8ef;--mut:#9aa;--line:#2a2a36;
--teal:#7fd1c4;--amber:#e0a96d;--red:#e07a7a;--green:#7fc98a;--blue:#6aa6d8;--gold:#f2c94c;}
*{box-sizing:border-box}
html,body{max-width:100%}
body{margin:0;background:var(--bg);color:var(--ink);font-family:-apple-system,"Segoe UI",
"Noto Sans CJK TC","PingFang TC",Roboto,sans-serif;line-height:1.6;font-size:15px}
.wrap{display:grid;grid-template-columns:210px 1fr;max-width:1240px;margin:0 auto}
nav.toc{position:sticky;top:0;align-self:start;height:100vh;overflow:auto;padding:18px 12px;
border-right:1px solid var(--line);font-size:13px}
nav.toc a{display:block;color:var(--mut);text-decoration:none;padding:5px 8px;border-radius:6px}
nav.toc a:hover{background:var(--panel2);color:var(--ink)}
nav.toc .lv2{padding-left:18px;font-size:12px}
main{padding:24px 30px;min-width:0}
header.hd{border-bottom:1px solid var(--line);padding-bottom:14px;margin-bottom:8px}
h1{font-size:23px;margin:.2em 0}
h2{font-size:19px;margin:1.5em 0 .5em;padding-top:8px;border-top:1px solid var(--line)}
h3{font-size:16px;margin:1.1em 0 .4em;color:var(--teal)}
.sub{color:var(--mut);font-size:13px}
.badge{display:inline-block;padding:2px 9px;border-radius:11px;font-size:11.5px;font-weight:600;
vertical-align:middle;margin:0 3px}
.b-tier{background:#3a2d5a;color:#cdbfff}.b-partial{background:#5a4a1f;color:#f2d98a}
.b-hc{background:#234a32;color:#9fe3b0}.b-med{background:#4a4423;color:#e8d98a}
.b-gap{background:#4a2d2d;color:#f0aaaa}.b-loh{background:#23404a;color:#9fd6e3}
.b-src{background:#22222c;color:#889;font-weight:500;font-size:10.5px;cursor:help;border:1px solid var(--line)}
.verdict{background:linear-gradient(180deg,#1a1a24,#15151c);border:1px solid var(--line);
border-radius:12px;padding:16px 18px;margin:10px 0}
.vquote{border-left:3px solid var(--teal);padding:6px 12px;color:#cfeee7;background:#13201e;
border-radius:0 8px 8px 0;margin:8px 0}
.vnot{border-left:3px solid var(--red);padding:6px 12px;color:#f2cccc;background:#201313;
border-radius:0 8px 8px 0;margin:8px 0}
.grid{display:grid;gap:12px;min-width:0}
.grid>*{min-width:0}
.cards{grid-template-columns:repeat(auto-fill,minmax(min(330px,100%),1fr))}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.card.gap{border-color:#5a3a3a}.card h4{margin:.1em 0 .4em;font-size:15px}
.kv{display:flex;justify-content:space-between;gap:8px;flex-wrap:wrap;min-width:0;font-size:13px;padding:2px 0;border-bottom:1px dotted #24242e}
.kv>*{min-width:0;overflow-wrap:anywhere}.kv .k{color:var(--mut)}.kv .v{text-align:right;font-variant-numeric:tabular-nums}
table{display:block;max-width:100%;overflow-x:auto;border-collapse:collapse;width:100%;font-size:13px;margin:6px 0}
th,td{border:1px solid var(--line);padding:5px 8px;text-align:center;font-variant-numeric:tabular-nums}
th{background:var(--panel2);color:var(--mut);font-weight:600}
td.l,th.l{text-align:left}
.conf{background:#2a1d1d}.clean{background:#1d2a20}.strong{outline:2px solid var(--gold);outline-offset:-2px}
.fig{max-width:100%;overflow-x:auto;background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px;margin:8px 0}
.fig figure{margin:0;min-width:0}.fig svg{display:block;max-width:100%;height:auto}
.fig figcaption{color:var(--mut);font-size:12px;margin-top:6px}
details.cav{background:var(--panel);border:1px solid var(--line);border-radius:9px;padding:8px 12px;margin:6px 0}
details.cav summary{cursor:pointer;font-weight:600;color:var(--amber)}
.judge{display:inline-flex;gap:4px;margin-top:8px}
.judge button{background:var(--panel2);border:1px solid var(--line);color:var(--mut);border-radius:7px;
padding:3px 9px;cursor:pointer;font-size:12px}
.judge button.on-ok{background:#234a32;color:#9fe3b0;border-color:#3a6e4a}
.judge button.on-doubt{background:#4a4423;color:#f2d98a;border-color:#6e6433}
.judge button.on-no{background:#4a2d2d;color:#f0aaaa;border-color:#6e3a3a}
.bar{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:8px 0}
.bar button{background:var(--teal);color:#08110f;border:none;border-radius:7px;padding:6px 12px;
font-weight:600;cursor:pointer}.bar button.ghost{background:var(--panel2);color:var(--ink);border:1px solid var(--line)}
.sub,.note,footer,code{overflow-wrap:anywhere;word-break:break-word}
.note{color:var(--mut);font-size:12.5px}
footer{margin-top:30px;border-top:1px solid var(--line);padding-top:12px;color:var(--mut);font-size:12px}
code{background:#22222c;padding:1px 5px;border-radius:4px;font-size:12.5px}
@media(max-width:820px){.wrap{grid-template-columns:1fr}nav.toc{display:none}main{padding:16px 12px}.grid[style]{grid-template-columns:1fr!important}}
"""

def src_badge(s):
    return f'<span class="b-src" title="provenance: {esc(s)}">⌖ {esc(s)}</span>'

def status_badge(st):
    if "HighConf" in st: return '<span class="badge b-hc">SEQC2 HighConf</span>'
    if "MedConf" in st: return '<span class="badge b-med">SEQC2 MedConf</span>'
    if "out_of_hc" in st or "out-of-HC" in st: return '<span class="badge b-gap">out-of-HC・truth-unevaluable（非 FP）</span>'
    return f'<span class="badge b-med">{esc(st)}</span>'

VK = {"confirm": ("#234a32", "#9fe3b0"), "bounded": ("#234a4a", "#9fe3e0"),
      "negated": ("#4a2d2d", "#f0aaaa"), "rewrite": ("#4a3a23", "#f2d098"),
      "provisional": ("#3a3a23", "#e8e08a")}

def build_html():
    source_commit = PROVENANCE["source_commit"]
    base_short = PROVENANCE["base_commit"][:12] if PROVENANCE["base_commit"] != "UNKNOWN" else "UNKNOWN"
    build_label = source_commit[:12] if source_commit else f"UNCOMMITTED_CANDIDATE atop {base_short}"
    # ----- sSNV cards
    cards = []
    for c in snv_cards:
        gapcls = " gap" if "out_of_hc" in c["seqc2_status"] else ""
        cards.append(f"""
<div class="card{gapcls}" data-snv="{c['idx']}">
  <h4>({c['idx']}) chr2:{c['pos']:,} <span class="note">{esc(c['ref'])}&gt;{esc(c['alt'])}</span></h4>
  <div>{status_badge(c['seqc2_status'])}</div>
  <div class="kv"><span class="k">歷史操作性標籤（非譜系）</span><span class="v">{esc(c['role'])}</span></div>
  <div class="kv"><span class="k">SEQC2 TVAF / 歷史近似轉換</span><span class="v">{esc(c['tvaf'])} / {esc(c['ccf'])} <span class="note">{esc('非 canonical CCF；不作 clone 判定' if c['ccf'] not in ('.', '', 'n/a') else 'n/a')}</span></span></div>
  <div class="kv"><span class="k">HKU tumor BQ20</span><span class="v">{esc(c['hku_t_bq20'])}</span></div>
  <div class="kv"><span class="k">DORADO tumor BQ20</span><span class="v">{esc(c['dor_t_bq20'])}</span></div>
  <div class="kv"><span class="k">HKU normal BQ20</span><span class="v">{esc(c['hku_n_bq20'])}</span></div>
  <div class="judge" data-snv="{c['idx']}">
    <span class="note" style="margin-right:4px">我的判讀:</span>
    <button data-v="ok">同意</button><button data-v="doubt">存疑</button><button data-v="no">否定</button>
  </div>
</div>""")
    # ----- rulings
    rule_rows = []
    for r in RULINGS:
        bg, fg = VK[r[3]]
        rule_rows.append(f'<tr><td>{r[0]}</td><td class="l">{esc(r[1])}</td>'
                         f'<td style="background:{bg};color:{fg};font-weight:600">{esc(r[2])}</td>'
                         f'<td class="l note">{esc(r[4])} {src_badge(r[5])}</td></tr>')
    # ----- retractions
    retr = "".join(f'<li><b>{esc(w)}</b> → {esc(n)} {src_badge(s)}</li>' for w, n, s in RETRACTIONS)
    # ----- linkage
    lrows = "".join(f'<tr><td class="l">{esc(x["pair"])}</td><td>{esc(x["hku"])}</td>'
                    f'<td>{esc(x["dor"])}</td><td class="l">{esc(x["read"])}</td></tr>' for x in LINKS)
    # ----- methyl table
    mrows = []
    for r in methyl_rows:
        cls = "conf" if r["class"] == "CONFOUNDED" else "clean"
        st = " strong" if r["strongest"] else ""
        mrows.append(
            f'<tr class="{cls}{st}"><td class="l">{esc(r["cpg"])}</td>'
            f'<td>{f2(r["hku_a"])}/{f2(r["hku_b"])}</td><td>{fq(r["hku_fdr"])}</td>'
            f'<td>{f2(r["dor_a"])}/{f2(r["dor_b"])}</td><td>{fq(r["dor_fdr"])}</td>'
            f'<td>{f2(r["norm_hp1"])}/{f2(r["norm_hp2"])}</td><td>{fq(r["norm_fdr"])}</td>'
            f'<td>{"germline-ASM" if r["class"]=="CONFOUNDED" else ("SCREEN-NEG★" if r["strongest"] else "SCREEN-NEG")}</td>'
            f'<td>{"✓" if r["dorado_replicates"] else "—"}</td></tr>')
    # ----- pos4
    p4 = POS4
    # ----- caveats
    cavs = "".join(f'<details class="cav"><summary>{esc(t)}</summary><div class="note">{esc(b)} {src_badge(s)}</div></details>'
                   for t, b, s in CAVEATS)
    # ----- conflicts
    crows = "".join(f'<tr><td class="l">{esc(m)}</td><td class="note">{esc(a)}</td>'
                    f'<td style="color:var(--green)">{esc(au)}</td><td class="l">{src_badge(s)}</td></tr>'
                    for m, a, au, s in CONFLICTS)

    loh_hku = loh["hku"]; loh_dor = loh["dorado"]
    hkanc = loh["hku_anchor"]; danc = loh["dor_anchor"]
    input_receipts = " · ".join(
        f"{name}: {meta['machine_path']} sha256={meta['sha256']}"
        for name, meta in PROVENANCE["input_files"].items()
    )

    body = f"""
<header class="hd">
  <h1>🔬 chr2:18M 分子狀態候選 — 判讀工作站 <span class="badge b-tier">tier L2</span> <span class="badge b-partial">PARTIAL</span></h1>
  <div class="sub">HCC1395 · {esc(AUD['region'])} · 6 sSNV × 10 CpG · HKU(5mCG_5hmCG haplotag) + DORADO cross-basecaller ·
  權威 SoT = <code>independent_audit.json</code> + <code>verdict_02</code> · revision <code>{esc(build_label)}</code></div>
  <div class="sub" style="margin-top:6px">
    規格所列 metric 由 generator 從來源檔注入，缺必填值即 refuse；
    hover <span class="b-src">⌖ 來源</span> 看出處。下方 sSNV 卡可記你的判讀，右上匯出 JSON/CSV。
  </div>
</header>

<div class="verdict">
  <h3 style="margin-top:0">最終可防守結論（verdict_02 §最終裁決）</h3>
  <div class="vquote">{esc(FINAL_DEFENSIBLE)}</div>
  <div class="vnot">Claim ceiling：{esc(CLAIM_CEILING)}</div>
  <div class="kv"><span class="k">適合論文用語</span><span class="v" style="color:var(--teal)">{esc(FINAL_PAPER)}</span></div>
  <div class="vnot">✗ 不適合：{esc(FINAL_NOT)}</div>
</div>

<h2 id="snv">① 6 個 sSNV 證據卡</h2>
<p class="note">tumor/normal allele 為 unique primary reads BQ20 count。pos3 落 HC 空隙＝truth 無法評估（<b>非 FP</b>）。
來源 <span class="b-src">⌖ independent_audit.json: snvs/benchmark_context/allele_counts_bq20</span> + <span class="b-src">⌖ concordance.tsv</span>。</p>
<div class="grid cards">{''.join(cards)}</div>

<h2 id="backbone">② α / β-strict 操作性 read-state labels</h2>
<div class="grid cards">
  <div class="card"><h4>operational read-state counts</h4>
    <div class="kv"><span class="k">HKU α(pos3=A) / β-strict</span><span class="v">{hkanc['alpha_pos3_A']} / {hkanc['beta_strict_E1_or_E2_or_E6']}</span></div>
    <div class="kv"><span class="k">HKU α∧β 違反 read</span><span class="v" style="color:var(--green)">{hkanc['alpha_and_beta_strict_violation']}</span></div>
    <div class="kv"><span class="k">DORADO α / β-strict</span><span class="v">{danc['alpha_pos3_A']} / {danc['beta_strict_E1_or_E2_or_E6']}</span></div>
    <div class="kv"><span class="k">DORADO α∧β 違反 read</span><span class="v" style="color:var(--green)">{danc['alpha_and_beta_strict_violation']}</span></div>
    <div class="note">⌖ independent_audit.json: samples.*.lineage_anchor_counts（source key 為歷史命名；此頁不作 lineage claim）</div>
  </div>
  <div class="card"><h4>判讀</h4>
    <div class="note">α = (3)A；β-strict = (1)G∨(2)C∨(6)G。於已評估 reads、此 strict 定義下，兩 basecaller 均觀測到 <b>0 strict 違反</b>。
    ⚠ 這不是 lineage 證明；broadened pos4-β 在 DORADO 有 1 discordant（見限制）。</div></div>
</div>
<h3>read-level linkage（00/10/01/11；11=兩事件同 read）</h3>
<table><thead><tr><th class="l">事件組合</th><th>HKU</th><th>DORADO</th><th class="l">判讀</th></tr></thead>
<tbody>{lrows}</tbody></table>
<p class="note">⌖ independent_audit.json: samples.{{HKU_T,DORADO_TAGGED_T}}.linkage_bq10 · 對照 verdict_02:140-147</p>

<h2 id="loh">③ SEQC2 LOH 區域標註 + tumor HP-tag imbalance</h2>
<div class="grid cards">
  <div class="card"><h4>HKU tumor HP</h4>
    <div class="kv"><span class="k">HP2-tag</span><span class="v">{loh_hku['HP2']}</span></div>
    <div class="kv"><span class="k">HP1-tag</span><span class="v">{loh_hku['HP1']}（{fpct(loh_hku['HP1_pct'],1)}）</span></div>
    <div class="kv"><span class="k">HP2 佔 tagged</span><span class="v" style="color:var(--teal)">{fpct(loh_hku['HP2_pct'],1)}</span></div></div>
  <div class="card"><h4>DORADO tumor HP</h4>
    <div class="kv"><span class="k">HP2-tag</span><span class="v">{loh_dor['HP2']}</span></div>
    <div class="kv"><span class="k">HP1-tag</span><span class="v">{loh_dor['HP1']}</span></div>
    <div class="kv"><span class="k">HP2 佔 tagged</span><span class="v" style="color:var(--teal)">{fpct(loh_dor['HP2_pct'],1)}</span></div></div>
  <div class="card"><h4>外部 / normal</h4>
    <div class="kv"><span class="k">SEQC2 LOH BED</span><span class="v"><span class="badge b-loh">{esc(loh['loh_interval'])}</span></span></div>
    <div class="kv"><span class="k">normal pos3 (A/G)</span><span class="v">0 / 18（100% REF）</span></div>
    <div class="note">all-REF 短 linkage read 存在（HKU ≥2pt={loh['hku_allref']['2']['all_ref']}）→「HP2 突變到不見」否定</div></div>
</div>
<p class="note">有界裁決：SEQC2 BED 標記此區為 LOH，且兩個 tumor tagged-read 集合皆呈 HP2 imbalance；本頁未整合 allele-specific CN，也不由此獨立確認 LOH 機制或「單親保留」。</p>

<h2 id="methyl">④ 甲基化三組（這是最容易誤讀的地方）</h2>
<p class="note">tumor α(pos3=A) vs β(pos3=G) 每 CpG meanM；右側 normal HP1/HP2 = matched-normal germline ASM。
<b>screen 由資料算出</b>：normal MW-FDR&lt;{CONF_FDR} ⇒ <span style="color:var(--red)">germline-ASM confounded</span>；否則 <span style="color:var(--green)">NORMAL-ASM-SCREEN-NEGATIVE</span>（不證明 acquired）。
DORADO✓ = DORADO tumor α/β FDR&lt;{CONF_FDR} 且同方向。</p>
<div class="grid" style="grid-template-columns:1.3fr 1fr">
  <div class="fig"><figure>{methyl_heatmap_svg()}<figcaption>圖4 · 甲基三組熱圖（藍=低甲基→紅=高甲基）。
  SCREEN-NEG★ = normal-ASM-screen-negative 且 DORADO 技術重現（{', '.join(strong_set)}；不證明 acquired/clone）。資料源 ⌖ independent_audit.json methyl_alpha_vs_pos3_ref + methyl_normal_hp1_vs_hp2。
  缺必填值會 refuse 而非顯示「—」（修掉舊 Fig4 的 v1/v2 漂移）。</figcaption></figure></div>
  <div>
    <div class="card"><h4>三組結論</h4>
      <div class="kv"><span class="k clean" style="padding:1px 5px;border-radius:4px">normal-ASM-screen-negative</span><span class="v">{esc(', '.join(clean_set))} <span class="note">（tumor-associated；非 acquired 證明）</span></span></div>
      <div class="kv"><span class="k conf" style="padding:1px 5px;border-radius:4px">CONFOUNDED（germline-ASM）</span><span class="v">{esc(', '.join(conf_set))}</span></div>
      <div class="kv"><span class="k">DORADO 複製</span><span class="v">{esc(', '.join(repl_set))}</span></div>
      <div class="kv"><span class="k">screen-negative ∩ DORADO 技術重現</span><span class="v" style="color:var(--gold)">{esc(', '.join(strong_set))}</span></div>
      <div class="note" style="margin-top:6px">⚠ 報告A 曾寫 DORADO 複製={{2.1,2.2,3.1,3.2,3.3}}，已修正（見限制末表）。</div>
    </div>
  </div>
</div>
<table><thead><tr><th class="l">CpG</th><th>HKU α/β</th><th>HKU FDR</th><th>DOR α/β</th><th>DOR FDR</th>
<th>normal HP1/HP2</th><th>normal FDR</th><th>分類</th><th>DOR✓</th></tr></thead><tbody>{''.join(mrows)}</tbody></table>

<h2 id="pos4">⑤ pos4 G/T/DEL 為何不足以支持三個 cellular subclone</h2>
<p class="note">chr2:18,096,269 鄰接約 20 bp poly-T homopolymer；觀測到的 G/T/DEL 混合 calls 與 homopolymer-associated calling uncertainty 相容，不能當作三個 cellular subclone 的證據。</p>
<table><thead><tr><th>Dataset</th><th>REF C</th><th>G</th><th>T</th><th>DEL</th></tr></thead><tbody>
<tr><td class="l">HKU BQ20</td><td>{p4['hku']['REF_C']}</td><td>{p4['hku']['G']}</td><td>{p4['hku']['T']}</td><td>{p4['hku']['DEL']}</td></tr>
<tr><td class="l">DORADO BQ20</td><td>{p4['dor']['REF_C']}</td><td>{p4['dor']['G']}</td><td>{p4['dor']['T']}</td><td>{p4['dor']['DEL']}</td></tr>
</tbody></table>
<p class="note">HKU 的 G 僅 4 條且全 reverse strand → 不可視為乾淨獨立群。⌖ independent_audit.json: pos4_subtype_summary（[推論6] 否定）。</p>

<h2 id="ccf">⑥ 歷史、案例特定的 TVAF 近似轉換（非 canonical CCF）</h2>
<div class="fig"><figure>{ccf_bar_svg()}<figcaption>圖 5 · 這些值是早期以 SEQC2 consensus TVAF 套入簡化式的案例特定近似轉換。現行 canonical 分析未整合 allele-specific CN、LOH、purity、multiplicity，因此本圖<b>不是 CN/LOH-corrected CCF，不能作 clone 結構佐證或 clonal/minor 分類</b>。in_truth=5/6（pos3 無外部真值）。⌖ chr2_18M_seqc2_concordance.tsv。</figcaption></figure></div>

<h2 id="tree">⑦ 有界 local read-state candidate graph + 6 推論裁決</h2>
<div class="fig"><figure>{tree_svg()}<figcaption>圖6 · 非時間、無方向的 local read-state candidate graph。圖形只整理共現／未共現與資料缺口；不代表 cellular clone、parent-child 或 linear ancestry。</figcaption></figure></div>
<table><thead><tr><th>#</th><th class="l">推論</th><th>裁決</th><th class="l">修正版 / 依據</th></tr></thead><tbody>{''.join(rule_rows)}</tbody></table>
<h3>3 處過度宣稱已撤回</h3><ul>{retr}</ul>

<h2 id="limit">⑧ 限制 + 報告A vs 權威數字校正</h2>
{cavs}
<h3>報告A（草稿）↔ 權威值 校正表</h3>
<table><thead><tr><th class="l">指標</th><th>報告A（已過時）</th><th>權威值</th><th class="l">來源</th></tr></thead><tbody>{crows}</tbody></table>

<footer>
  <b>Provenance</b> · revision state <code>{esc(PROVENANCE['revision_state'])}</code> · source commit
  <code>{esc(PROVENANCE['source_commit']) if PROVENANCE['source_commit'] else 'null'}</code> · base commit
  <code>{esc(PROVENANCE['base_commit'])}</code> · branch <code>{esc(PROVENANCE['branch'])}</code> · 產生器
  <code>build_workstation_html.py</code> · 資料 <code>workstation_data.json</code><br>
  original producer: <code>{esc(ORIGINAL_PRODUCER['name'])}</code> · original commit
  <code>{esc(ORIGINAL_PRODUCER['commit'])}</code> · source_created_at
  <code>{esc(ORIGINAL_PRODUCER['source_created_at'])}</code><br>
  generated_at: <code>{esc(PROVENANCE['generated_at'])}</code> · published_at: <code>null</code> ·
  dirty revision inputs: <code>{esc(', '.join(PROVENANCE['dirty_revision_inputs']) or 'none')}</code><br>
  last_verified_at: <code>{esc(PROVENANCE['last_verified_at'])}</code> · verification_commit: <code>null</code> ·
  verification_receipt: <code>{esc(PROVENANCE['verification_receipt'])}</code><br>
  science_recomputed: <code>false</code> · {esc(PROVENANCE['claim_ceiling_note'])}<br>
  input receipts: <code>{esc(input_receipts)}</code><br>
  來源檔：<code>{esc(AUDIT_REL)}</code> · <code>{esc(VERDICT_REL)}</code> · <code>{esc(CONC_REL)}</code><br>
  scope: {esc(BUILD['scope'])}。生成器規格列出的數值欄位由 JSON/TSV 注入；qualitative 裁決經有界改寫並標 verdict_02 file:line。
  cross-basecaller = 技術重現非獨立 biological replicate（tier 封頂 L2）。
</footer>
"""
    datalit = json.dumps(DATAJSON, ensure_ascii=False)
    page = f"""<!DOCTYPE html>
<!--
  07 chr2:18M molecular-state candidate JUDGMENT WORKSTATION (generated, do not hand-edit)
  revision_state: {esc(PROVENANCE['revision_state'])}
  source_commit: {esc(PROVENANCE['source_commit']) if PROVENANCE['source_commit'] else 'null'}
  base_commit: {esc(PROVENANCE['base_commit'])}
  branch: {esc(PROVENANCE['branch'])}
  original_producer: {esc(ORIGINAL_PRODUCER['name'])} @ {esc(ORIGINAL_PRODUCER['commit'])}
  source_created_at: {esc(ORIGINAL_PRODUCER['source_created_at'])}
  generated_at: {BUILD['built_utc']}
  last_verified_at: {esc(PROVENANCE['last_verified_at'])}
  verification_commit: null
  verification_receipt: {esc(PROVENANCE['verification_receipt'])}
  published_at: null
  science_recomputed: false
  claim_ceiling_note: {esc(PROVENANCE['claim_ceiling_note'])}
  regenerate: python3 .../verification_assets/scripts/build_workstation_html.py --html-only
  data_sources: {AUDIT_REL} ; {VERDICT_REL} ; {CONC_REL}
  machine_data_dir: {esc(PROVENANCE['data_dir'])}
  input_receipts: {esc(input_receipts)}
  required-field-policy: 規格列出的數值欄位由 independent_audit.json/concordance.tsv 注入；缺必填值即 refuse
-->
<html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>chr2:18M 分子狀態候選判讀工作站 · HCC1395 · L2</title>
<style>{CSS}</style></head>
<body>
<div class="wrap">
<nav class="toc">
  <div style="font-weight:700;color:var(--ink);padding:4px 8px">chr2:18M 工作站</div>
  <div class="bar" style="padding:0 6px">
    <button onclick="exportJSON()">⤓ JSON</button>
    <button class="ghost" onclick="exportCSV()">⤓ CSV</button>
  </div>
  <a href="#snv">① 6 sSNV 證據卡</a>
  <a href="#backbone">② α/β read-state labels</a>
  <a href="#loh">③ LOH annotation + HP imbalance</a>
  <a href="#methyl">④ 甲基三組</a>
  <a href="#pos4">⑤ pos4 homopolymer</a>
  <a href="#ccf">⑥ 歷史 TVAF 近似轉換</a>
  <a href="#tree">⑦ candidate graph + 裁決</a>
  <a href="#limit">⑧ 限制 + 校正</a>
  <div class="note" style="padding:10px 8px">tier L2 · PARTIAL<br>{esc(build_label)}</div>
</nav>
<main>{body}</main>
</div>
<script id="wsdata" type="application/json">{datalit}</script>
<script>{JS}</script>
</body></html>"""
    return page

JS = """
const DATA = JSON.parse(document.getElementById('wsdata').textContent);
const LSKEY = 'ism_chr2_18m_ws_v1';
function loadV(){try{return JSON.parse(localStorage.getItem(LSKEY))||{}}catch(e){return {}}}
function saveV(v){localStorage.setItem(LSKEY, JSON.stringify(v))}
function paint(){
  const v = loadV();
  document.querySelectorAll('.judge').forEach(j=>{
    const id = j.dataset.snv;
    j.querySelectorAll('button').forEach(b=>{
      b.classList.remove('on-ok','on-doubt','on-no');
      if(v[id]===b.dataset.v) b.classList.add('on-'+b.dataset.v);
    });
  });
}
document.querySelectorAll('.judge button').forEach(b=>{
  b.addEventListener('click',()=>{
    const id=b.parentElement.dataset.snv, v=loadV();
    v[id] = (v[id]===b.dataset.v)? null : b.dataset.v;
    saveV(v); paint();
  });
});
paint();
function snapshot(){
  const v=loadV();
  return DATA.snv_cards.map(c=>({snv:c.idx, pos:c.pos, ref:c.ref, alt:c.alt,
    seqc2:c.seqc2_status, role:c.role, ccf:c.ccf, my_verdict:v[c.idx]||'(未判)'}));
}
function dl(name, text, mime){
  const a=document.createElement('a');
  a.href=URL.createObjectURL(new Blob([text],{type:mime})); a.download=name; a.click();
}
function exportJSON(){
  dl('chr2_18M_my_verdicts.json', JSON.stringify({build:DATA.build, verdicts:snapshot()},null,2),'application/json');
}
function exportCSV(){
  const rows=snapshot(); const hdr=Object.keys(rows[0]);
  const csv=[hdr.join(',')].concat(rows.map(r=>hdr.map(h=>JSON.stringify(r[h]??'')).join(','))).join('\\n');
  dl('chr2_18M_my_verdicts.csv', csv, 'text/csv');
}
"""

# ----------------------------------------------------------------------------- write
html_doc = build_html()
os.makedirs(EXPLAIN, exist_ok=True)
with open(HTML_OUT, "w", encoding="utf-8") as fh:
    fh.write(html_doc)

if not ARGS.html_only:
    os.makedirs(DISPLAY, exist_ok=True)
    with open(JSON_OUT, "w", encoding="utf-8") as fh:
        json.dump(DATAJSON, fh, ensure_ascii=False, indent=2)

    # data-bound Fig4 fragment for page-04 injection
    FIG4_OUT = os.path.join(DISPLAY, "fig4_corrected.svg")
    with open(FIG4_OUT, "w", encoding="utf-8") as fh:
        fh.write(fig4_page04_svg())

print("[ok] html      ->", os.path.relpath(HTML_OUT, REPO))
if ARGS.html_only:
    print("[skip] --html-only: workstation_data.json and fig4_corrected.svg unchanged")
else:
    print("[ok] data.json ->", os.path.relpath(JSON_OUT, REPO))
    print("[ok] fig4 svg  ->", os.path.relpath(FIG4_OUT, REPO))
print(f"[provenance] state={PROVENANCE['revision_state']} source_commit={PROVENANCE['source_commit']} "
      f"base_commit={PROVENANCE['base_commit']} published_at={PROVENANCE['published_at']}")
print(f"[classify] clean={clean_set}")
print(f"[classify] confounded={conf_set}")
print(f"[classify] dorado_replicate={repl_set}")
print(f"[classify] strongest(clean∩replicate)={strong_set}")
print(f"[loh] HKU HP2/HP1={loh['hku']['HP2']}/{loh['hku']['HP1']} ({loh['hku']['HP2_pct']*100:.1f}%) ; "
      f"DORADO {loh['dorado']['HP2']}/{loh['dorado']['HP1']} ({loh['dorado']['HP2_pct']*100:.1f}%)")
print(f"[backbone] HKU violation={loh['hku_anchor']['alpha_and_beta_strict_violation']} ; "
      f"DORADO violation={loh['dor_anchor']['alpha_and_beta_strict_violation']}")
