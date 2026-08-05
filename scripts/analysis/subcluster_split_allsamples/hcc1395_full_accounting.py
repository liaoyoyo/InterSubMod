#!/usr/bin/env python3
"""HCC1395 per-SNV-locus 全面整理 — 切割(無監督) + 統計(PERMANOVA) + 方法學分類(C++ verdict) + cross-tab。
零捏造: 全讀 significance_summary_HCC1395.csv (180欄, C++ 真值) + records_wg2_HCC1395.json (best_k)。
輸出: hcc1395_full_accounting.json
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    extract_provenance_frame,
    select_current_view,
    select_legacy_view,
    select_loh_legacy,
    validate_region_strata,
)

DEFAULT_BASE = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
    "20260620_allsample_subcluster_split/results"
)
parser = argparse.ArgumentParser(description="Build schema-aware HCC1395 accounting JSON.")
parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
parser.add_argument("--csv", type=Path, default=None)
parser.add_argument("--records", type=Path, default=None)
parser.add_argument("--output", type=Path, default=None)
parser.add_argument("--region-status", type=Path, default=None)
parser.add_argument(
    "--allow-unversioned-v1",
    action="store_true",
    help="Explicitly authorize the frozen historical v1 VerificationClass/LOH fields; region strata remain unavailable.",
)
args = parser.parse_args()

BASE = args.base
CSV = args.csv or (BASE / "significance_summary_HCC1395.csv")
RECS = args.records or (BASE / "records_wg2_HCC1395.json")
OUTPUT = args.output or (BASE / "hcc1395_full_accounting.json")

with CSV.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle))
N_sig = len(rows)
with RECS.open("r", encoding="utf-8") as handle:
    recs = json.load(handle)
N_rec = len(recs)
# best_k 查表 (切割維度) keyed by Chr_Pos
bestk = {f"{r['chrom']}_{r['pos']}": r["all"]["best_k"] for r in recs}
nall = {f"{r['chrom']}_{r['pos']}": r["all"]["n"] for r in recs}

frame = pd.DataFrame(rows)
try:
    if "VerificationSchemaVersion" in frame.columns:
        current_view = select_current_view(frame)
        legacy_view = select_legacy_view(frame)
        loh_view = select_loh_legacy(frame)
        extract_provenance_frame(frame)

        status_path = args.region_status or (CSV.parent / "region_stratification_status.tsv")
        if not status_path.is_file():
            raise SchemaContractError(
                "v2 HCC1395 accounting requires region_stratification_status.tsv; "
                "offline Verification-only migration is not a RegionStratifier run"
            )
        with status_path.open("r", encoding="utf-8", newline="") as handle:
            status_rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(status_rows) != 1:
            raise SchemaContractError("region status must contain exactly one data row")
        region_view = validate_region_strata(frame, status_rows[0])
        region_status = region_view.status
    else:
        if not args.allow_unversioned_v1:
            raise SchemaContractError(
                "historical input has no VerificationSchemaVersion; rerun only with explicit "
                "--allow-unversioned-v1"
            )
        current_view = select_current_view(frame, allow_unversioned_raw=True)
        legacy_view = select_legacy_view(
            frame, allow_unversioned_v1=True, unversioned_unknown_policy="fail"
        )
        loh_view = select_loh_legacy(frame, allow_unversioned_v1=True)
        region_status = {
            "status": "UNVERSIONED_REGION_SCHEMA",
            "reason": "HISTORICAL_SUBCLONE_ID_NOT_INTERPRETED",
        }
except SchemaContractError as error:
    raise SystemExit(f"REFUSE: schema contract failed: {error}") from error

for index, row in enumerate(rows):
    row["_VerificationClass_Current_Selected"] = str(current_view.values.iloc[index])
    row["_VerificationClass_Legacy_Selected"] = str(legacy_view.values.iloc[index])
    row["_LOH_Subtype_Legacy_Selected"] = str(loh_view.values.iloc[index])

def tru(v): return str(v).strip().lower() in ("true", "1", "yes")
def ff(v):
    try: return float(v)
    except: return None

MIN_SZ = 3
def dist(col):
    c = Counter(str(r.get(col, "")).strip() for r in rows)
    return {k: v for k, v in c.most_common()}
def pct(n, tot=N_sig): return round(100 * n / tot, 1) if tot else None

out = {
    "sample": "HCC1395",
    "N_significance": N_sig,
    "N_records": N_rec,
    "note": "切割維度以 records_wg2 best_k; 分類/統計維度以 significance_summary C++ 真值。",
    "schema_contract": {
        "verification_selection_field": current_view.field,
        "verification_schema_status": current_view.schema_status,
        "verification_categories": list(current_view.categories),
        "verification_unknown_counts": current_view.unknown_counts,
        "legacy_selection_field": legacy_view.field,
        "legacy_schema_status": legacy_view.schema_status,
        "legacy_unknown_counts": legacy_view.unknown_counts,
        "loh_selection_field": loh_view.field,
        "loh_schema_status": loh_view.schema_status,
        "region_status": region_status,
        "source_csv": str(CSV.resolve()),
    },
}

# ===== 維度 1: 切割(無監督 UPGMA+silhouette) =====
split = Counter()
for r in rows:
    key = f"{r['Chr']}_{r['Pos']}"
    bk = bestk.get(key); n = nall.get(key)
    if key not in bestk:
        split["no_record"] += 1
    elif n is not None and n < MIN_SZ * 2:
        split["insuff_n<6"] += 1
    elif bk is None or bk < 2:
        split["cant_切不出"] += 1
    else:
        split[f"cansplit_k={min(bk,5)}"] += 1
cansplit_total = sum(v for k, v in split.items() if k.startswith("cansplit"))
out["dim1_split"] = {"detail": dict(split), "cansplit_total": cansplit_total, "cansplit_pct": pct(cansplit_total),
                     "cant_total": split["cant_切不出"], "cant_pct": pct(split["cant_切不出"]),
                     "insuff": split["insuff_n<6"], "insuff_pct": pct(split["insuff_n<6"])}

# ===== 維度 2: 統計 PERMANOVA / GlobalTest =====
def axis_stat(pcol, vcol, wcol):
    sig = clean = disp = 0
    for r in rows:
        v = tru(r.get(vcol)); p = ff(r.get(pcol)); w = tru(r.get(wcol))
        if v and p is not None and p < 0.05:
            sig += 1
            if w: disp += 1
            else: clean += 1
    return {"sig": sig, "sig_pct": pct(sig), "clean": clean, "clean_pct": pct(clean), "disp": disp, "disp_pct": pct(disp)}
out["dim2_permanova"] = {
    "LabelHP": axis_stat("LabelHPPermanovaP", "LabelHPPermanovaValid", "LabelHPDispersionWarn"),
    "LabelAllele": axis_stat("LabelAllelePermanovaP", "LabelAllelePermanovaValid", "LabelAlleleDispersionWarn"),
    "Cluster": axis_stat("ClusterPermanovaP", "ClusterPermanovaValid", "ClusterDispersionWarn"),
}
# C = clean (任一 label 軸)
C = 0
for r in rows:
    hv = tru(r.get("LabelHPPermanovaValid")); hp = ff(r.get("LabelHPPermanovaP")); hw = tru(r.get("LabelHPDispersionWarn"))
    av = tru(r.get("LabelAllelePermanovaValid")); ap = ff(r.get("LabelAllelePermanovaP")); aw = tru(r.get("LabelAlleleDispersionWarn"))
    if (hv and hp is not None and hp < 0.05 and not hw) or (av and ap is not None and ap < 0.05 and not aw):
        C += 1
out["dim2_permanova"]["C_clean_anyaxis"] = {"n": C, "pct": pct(C)}
# GlobalP / CramersV bands
gsig = sum(1 for r in rows if (ff(r.get("GlobalP")) is not None and ff(r["GlobalP"]) < 0.05))
out["dim2_permanova"]["GlobalP_lt0.05"] = {"n": gsig, "pct": pct(gsig)}
cv = Counter()
for r in rows:
    v = ff(r.get("CramersV"))
    if v is None: cv["NA"] += 1
    elif v >= 0.7: cv["≥0.7"] += 1
    elif v >= 0.5: cv["0.5-0.7"] += 1
    elif v >= 0.3: cv["0.3-0.5"] += 1
    else: cv["<0.3"] += 1
out["dim2_permanova"]["CramersV_bands"] = dict(cv)

# ===== 維度 3: 方法學分類(C++ verdict) =====
out["dim3_classification"] = {
    "VerificationClass_Current": dist("_VerificationClass_Current_Selected"),
    "VerificationClass_Legacy_4state": dist("_VerificationClass_Legacy_Selected"),
    "StrengthGrade_AtoE": dist("StrengthGrade"),
    "PassedGating": dist("PassedGating"),
    "Significant": dist("Significant"),
    "SuggestFilter": dist("SuggestFilter"),
    "DominantLabel": dist("DominantLabel"),
    "TumorIntrinsic": dist("TumorIntrinsic"),
}

# ===== 維度 4: CNV / LOH / WithinHP subclone =====
out["dim4_context"] = {
    "Coverage_Category": dist("Coverage_Category"),
    "Potential_LOH": dist("Potential_LOH"),
    "LOH_Subtype_LegacyVC": dist("_LOH_Subtype_Legacy_Selected"),
    "Quality_Tier": dist("Quality_Tier"),
    "WithinHP_CleanMultigroup": dist("WithinHP_CleanMultigroup"),
    "WithinHP_LevelBimodal": dist("WithinHP_LevelBimodal"),
    "WithinHP_SubcloneSig": dist("WithinHP_SubcloneSig"),
    "TumorOnlyPermanova_sig": None,  # 下面填
}
to_sig = sum(1 for r in rows if tru(r.get("TumorOnlyPermanovaValid")) and ff(r.get("TumorOnlyPermanovaP")) is not None and ff(r["TumorOnlyPermanovaP"]) < 0.05)
out["dim4_context"]["TumorOnlyPermanova_sig"] = {"n": to_sig, "pct": pct(to_sig)}

# ===== 維度 5: A/C Venn (一致性) =====
A = set(); Cset = set()
for r in rows:
    key = f"{r['Chr']}_{r['Pos']}"
    bk = bestk.get(key)
    if bk is not None and bk >= 2: A.add(key)
    hv = tru(r.get("LabelHPPermanovaValid")); hp = ff(r.get("LabelHPPermanovaP")); hw = tru(r.get("LabelHPDispersionWarn"))
    av = tru(r.get("LabelAllelePermanovaValid")); ap = ff(r.get("LabelAllelePermanovaP")); aw = tru(r.get("LabelAlleleDispersionWarn"))
    if (hv and hp is not None and hp < 0.05 and not hw) or (av and ap is not None and ap < 0.05 and not aw): Cset.add(key)
inter = len(A & Cset); uni = len(A | Cset)
out["dim5_venn"] = {"A_cansplit": len(A), "C_clean": len(Cset), "AnC": inter, "A_only": len(A - Cset),
                    "C_only": len(Cset - A), "jaccard": round(inter / uni, 3) if uni else None,
                    "A_pct": pct(len(A)), "C_pct": pct(len(Cset))}

# ===== 維度 6: cross-tabs =====
def crosstab(rowkey_fn, colkey_fn):
    t = defaultdict(lambda: Counter())
    for r in rows:
        t[rowkey_fn(r)][colkey_fn(r)] += 1
    return {k: dict(v) for k, v in t.items()}

def split_state(r):
    key = f"{r['Chr']}_{r['Pos']}"; bk = bestk.get(key); n = nall.get(key)
    if key not in bestk: return "no_record"
    if n is not None and n < 6: return "insuff"
    if bk is None or bk < 2: return "cant"
    return "cansplit"

# cansplit × VerificationClass
out["dim6_crosstab"] = {
    "split_x_VerificationClass": crosstab(
        split_state, lambda r: str(r.get("_VerificationClass_Current_Selected", "")).strip()
    ),
    "split_x_Significant": crosstab(split_state, lambda r: "sig" if tru(r.get("Significant")) else "nonsig"),
    "split_x_StrengthGrade": crosstab(split_state, lambda r: str(r.get("StrengthGrade", "")).strip()),
    "VerificationClass_x_StrengthGrade": crosstab(
        lambda r: str(r.get("_VerificationClass_Current_Selected", "")).strip(),
        lambda r: str(r.get("StrengthGrade", "")).strip(),
    ),
    "split_x_CoverageCategory": crosstab(split_state, lambda r: str(r.get("Coverage_Category", "")).strip()),
}

# ===== mode 偵測 (tumor-only vs paired) =====
mode = "unknown"
try:
    for ln in open(f"{BASE}/ism_run_HCC1395.log"):
        if "Normal BAM" in ln:
            mode = "tumor-only" if "None" in ln else "paired"
            break
except FileNotFoundError:
    pass
nn0 = sum(1 for r in rows if str(r.get("NNormalReads", "")).strip() in ("0", "", "0.0"))
out["mode"] = {"run_mode": mode, "NNormalReads_zero": nn0, "NNormalReads_zero_pct": pct(nn0)}

# ===== dim2b: CramersV exactly-0 拆解 =====
cv0 = sum(1 for r in rows if (ff(r.get("CramersV")) is not None and abs(ff(r["CramersV"])) < 1e-9))
out["dim2_permanova"]["CramersV_exactly0_gated"] = {"n": cv0, "pct": pct(cv0),
    "note": "Cochran-gated 或無關聯→V=0；佔 <0.3 帶絕大多數，非連續弱關聯"}

# ===== dim3b: Significant readback + Significant×StrengthGrade =====
# Significant is authoritative. Do not reconstruct truth from class names.
sig_flag = out["dim3_classification"]["Significant"].get("true", 0)
sigE = sum(1 for r in rows if tru(r.get("Significant")) and r["StrengthGrade"].strip() == "E")
out["dim3_classification"]["Significant_reconciliation"] = {
    "significant_flag": sig_flag,
    "source_field": "Significant",
    "class_name_reconstruction_performed": False,
    "of_which_StrengthGrade_E": sigE,
    "E_pct_of_sig": round(100 * sigE / max(1, sig_flag), 1),
    "note": "Significant 直接讀 canonical boolean；不從 VerificationClass 名稱重建 truth。",
}

# ===== dim5b: 三透鏡 (cansplit / C_clean / Significant) 非巢狀重疊 =====
Sset = set(f"{r['Chr']}_{r['Pos']}" for r in rows if tru(r.get("Significant")))
out["dim5_venn"]["three_lens"] = {
    "cansplit": len(A), "C_clean": len(Cset), "Significant": len(Sset),
    "cansplit_and_Sig": len(A & Sset), "Sig_only": len(Sset - A), "cansplit_only": len(A - Sset),
    "cansplit_and_Cclean": len(A & Cset),
    "note": "三者部分重疊非漏斗；cansplit=無監督k≥2(可為cis-ASM/germline-HP/CN非subclone)",
}

# ===== dim7: 補軸 (HP_AUC / Fisher per-CpG / NHP_Somatic / TumorOnlyK / StrengthScore 分量) =====
def auc_bands(col):
    b = Counter()
    for r in rows:
        v = ff(r.get(col))
        if v is None: b["NA"] += 1
        elif v >= 0.7: b["≥0.7"] += 1
        elif v >= 0.6: b["0.6-0.7"] += 1
        elif v >= 0.5: b["0.5-0.6"] += 1
        else: b["<0.5"] += 1
    return dict(b)
fisher_nonzero = sum(1 for r in rows if (ff(r.get("Fisher_N_Sig")) or 0) > 0)
def nz(col): return sum(1 for r in rows if (ff(r.get(col)) or 0) > 0)
def mean_of(col):
    vs = [ff(r.get(col)) for r in rows if ff(r.get(col)) is not None]
    return round(sum(vs) / len(vs), 3) if vs else None
out["dim7_extra_axes"] = {
    "HP_AUC_Tumor_bands": auc_bands("HP_AUC_Tumor"),
    "HP_AUC_All_bands": auc_bands("HP_AUC_All"),
    "Fisher_perCpG_ASM_nonzero": {"n": fisher_nonzero, "pct": pct(fisher_nonzero), "note": "至少一個 per-CpG Fisher 顯著位點"},
    "NHP_Somatic11_nonzero": {"n": nz("NHP_Somatic11"), "pct": pct(nz("NHP_Somatic11"))},
    "NHP_Somatic21_nonzero": {"n": nz("NHP_Somatic21"), "pct": pct(nz("NHP_Somatic21"))},
    "NHP_Somatic33_nonzero": {"n": nz("NHP_Somatic33"), "pct": pct(nz("NHP_Somatic33"))},
    "TumorOnlyClusterK": dist("TumorOnlyClusterK"),
    "StrengthScore_component_means": {c: mean_of(c) for c in
        ["StrengthScore", "StrengthStruct", "StrengthTumor", "StrengthSomatic", "StrengthAssoc", "StrengthGermline"]},
}

# ===== dim8: tumor-only 空欄稽核 (normal-anchored 欄位「不可計算」非「陰性」) =====
def blank_or_zero(col):
    c = 0
    for r in rows:
        v = str(r.get(col, "")).strip()
        if v in ("", "NA", "nan", "0", "0.0", "false", "-1"):
            c += 1
    return c
normal_anchored = ["SampleASM_NNormal", "SomaticResidualDbeta_Sig", "GermlineAsmDbeta_Sig",
                   "HP_Residual_Sig", "SubcloneDbeta_HP1_Sig", "SubcloneDbeta_HP2_Sig",
                   "SampleASM_Sig", "PerCpgASM_Valid", "NNormalReads"]
out["dim8_tumor_only_empty_audit"] = {
    "mode": mode,
    "columns_not_computable": {c: blank_or_zero(c) for c in normal_anchored},
    "region_status": region_status,
    "note": "tumor-only run：normal-anchored 欄位為『不可計算』非『陰性發現』；region assignment 依 canonical status/reason 解讀。",
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT.open("w", encoding="utf-8") as handle:
    json.dump(out, handle, indent=1, ensure_ascii=False)
# 摘要列印
print(f"N_sig={N_sig} N_rec={N_rec}")
print(f"切割: cansplit={out['dim1_split']['cansplit_total']}({out['dim1_split']['cansplit_pct']}%) cant={out['dim1_split']['cant_total']}({out['dim1_split']['cant_pct']}%) insuff={out['dim1_split']['insuff']}")
print(f"Venn: A={out['dim5_venn']['A_cansplit']} C={out['dim5_venn']['C_clean']} AnC={out['dim5_venn']['AnC']} Jaccard={out['dim5_venn']['jaccard']}")
print(f"VerificationClass: {out['dim3_classification']['VerificationClass_Current']}")
print(f"Significant true={out['dim3_classification']['Significant'].get('true')} PassedGating true={out['dim3_classification']['PassedGating'].get('true')}")
print(f"[json -> {OUTPUT}]")
