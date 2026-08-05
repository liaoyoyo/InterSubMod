#!/usr/bin/env python3
"""
TO TP 全量特徵畫像 + FP Germline/LOH/Artifact 分類分析
2026-03-22

分析目標：
1. TO 28,505 個 TP 的 AF 分佈 (Subclonal/Somatic/High-AF 亞型)
2. TO 11,598 個殘餘 PASS FP 按 normal AF 分類：
   - Germline-escaped (naf >= 0.2)
   - LOH-region FP (naf < 0.05, to_af >= 0.7)
   - Tumor-specific FP (naf < 0.05, to_af < 0.4)
   - Uncertain (paired_raw_ndp = 0)
3. ISM 候選池 FP (298) 的甲基化特徵差異
4. 高 AF TP 的 LOH 特徵探索
"""

import argparse
import gzip
import csv
import os
import sys
import math
from collections import defaultdict, Counter
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import select_legacy_view

# ──────────────────────────────────────────────────────────────────
# 路徑設定
# ──────────────────────────────────────────────────────────────────
BASE_DIR = "/big7_disk/liaoyoyo2001/InterSubMod"
PROVENANCE_GZ = os.path.join(
    BASE_DIR,
    "research/fp_provenance/20260322_hcc1395_5khz_to/hcc1395_to_fp_provenance_master.tsv.gz"
)
RESCUE_FEATURES = (
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/"
    "big8_disk_output/research_rounds/"
    "20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv"
)
OUTPUT_DIR = os.path.join(
    BASE_DIR,
    "research/fp_provenance/20260322_hcc1395_5khz_to"
)

# F1 baseline
BASELINE_TP = 28505
BASELINE_FP = 11598
BASELINE_FN = 10942
BASELINE_F1 = 2 * BASELINE_TP / (2 * BASELINE_TP + BASELINE_FP + BASELINE_FN)


def attach_historical_legacy_records(rows, columns=None, allow_unversioned_v1=False):
    """Attach an explicit legacy-four-state view to historical rescue rows."""
    frame = pd.DataFrame(rows, columns=columns)
    view = select_legacy_view(
        frame,
        allow_unversioned_v1=allow_unversioned_v1,
        unversioned_unknown_policy="fail",
    )
    selected_rows = []
    for row, selected_class in zip(rows, view.values.tolist()):
        selected = dict(row)
        selected["VerificationClass_Legacy_Selected"] = selected_class
        selected["verification_selection_field"] = view.field
        selected["verification_schema_status"] = view.schema_status
        selected_rows.append(selected)
    return selected_rows, {
        "verification_view": "legacy4_historical",
        "verification_selection_field": view.field,
        "verification_schema_status": view.schema_status,
    }

# ──────────────────────────────────────────────────────────────────
# 分類函數
# ──────────────────────────────────────────────────────────────────

def classify_fp(to_af, naf, ndp):
    """FP 類型分類"""
    try:
        naf_val = float(naf) if naf and naf not in ('', 'nan', 'NA', 'None') else None
        ndp_val = int(float(ndp)) if ndp and ndp not in ('', 'nan', 'NA', 'None') else 0
        af_val = float(to_af) if to_af and to_af not in ('', 'nan', 'NA', 'None') else None
    except (ValueError, TypeError):
        return "Uncertain"

    # 無 normal coverage 記錄
    if ndp_val == 0:
        return "Uncertain"

    # 有 normal coverage
    if naf_val is not None:
        if naf_val >= 0.2:
            return "Germline-escaped"
        elif naf_val < 0.05:
            if af_val is not None and af_val >= 0.7:
                return "LOH-region"
            elif af_val is not None and af_val < 0.4:
                return "Tumor-specific"
            else:
                return "Mid-AF-artifact"  # 0.4-0.7 range with low naf
        else:
            return "Uncertain-naf"  # 0.05-0.2 range
    return "Uncertain"


def classify_tp(to_af):
    """TP 亞型分類"""
    try:
        af = float(to_af)
    except (ValueError, TypeError):
        return "Unknown"
    if af < 0.2:
        return "Subclonal"
    elif af < 0.6:
        return "Somatic"
    else:
        return "High-AF"


def safe_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def safe_bool(val):
    return str(val).strip().lower() in ('true', '1', 'yes')


# ──────────────────────────────────────────────────────────────────
# Step 1: 讀取 Provenance Master
# ──────────────────────────────────────────────────────────────────

def load_provenance_data():
    """從 provenance_master.tsv.gz 提取 TP 和 FP 記錄"""
    print("[Step 1] 讀取 provenance_master.tsv.gz ...", flush=True)

    tp_records = []
    fp_records = []
    total = 0

    with gzip.open(PROVENANCE_GZ, 'rt', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            total += 1
            if total % 500000 == 0:
                print(f"  已讀取 {total:,} 行 ...", flush=True)

            truth_status = row.get('truth_status', '')
            primary_class = row.get('primary_class', '')
            to_baseline_pass = row.get('to_baseline_pass', '')

            # TP: 最終輸出的 TP
            if truth_status == 'TP' and safe_bool(to_baseline_pass):
                tp_records.append({
                    'variant_key': row.get('variant_key', ''),
                    'to_af': row.get('to_af', ''),
                    'to_gq': row.get('to_gq', ''),
                    'to_qual': row.get('to_qual', ''),
                    'to_dp': row.get('to_dp', ''),
                    'to_verdict_somatic': row.get('to_verdict_somatic', ''),
                    'to_verdict_subclonal': row.get('to_verdict_subclonal', ''),
                    'to_verdict_germline': row.get('to_verdict_germline', ''),
                    'paired_raw_naf': row.get('paired_raw_naf', ''),
                    'paired_raw_ndp': row.get('paired_raw_ndp', ''),
                    'paired_raw_af': row.get('paired_raw_af', ''),
                    'paired_resolution_stage': row.get('paired_resolution_stage', ''),
                })

            # FP: 殘餘 PASS FP
            elif primary_class == 'to_residual_final_fp':
                fp_records.append({
                    'variant_key': row.get('variant_key', ''),
                    'to_af': row.get('to_af', ''),
                    'to_gq': row.get('to_gq', ''),
                    'to_qual': row.get('to_qual', ''),
                    'to_dp': row.get('to_dp', ''),
                    'to_verdict_somatic': row.get('to_verdict_somatic', ''),
                    'to_verdict_subclonal': row.get('to_verdict_subclonal', ''),
                    'to_verdict_germline': row.get('to_verdict_germline', ''),
                    'paired_raw_naf': row.get('paired_raw_naf', ''),
                    'paired_raw_ndp': row.get('paired_raw_ndp', ''),
                    'paired_raw_af': row.get('paired_raw_af', ''),
                    'paired_resolution_stage': row.get('paired_resolution_stage', ''),
                })

    print(f"  完成。TP={len(tp_records):,}, FP={len(fp_records):,}, 總行數={total:,}", flush=True)
    return tp_records, fp_records


# ──────────────────────────────────────────────────────────────────
# Step 2: TP 特徵分析
# ──────────────────────────────────────────────────────────────────

def analyze_tp(tp_records):
    """TP 全量特徵分析"""
    print("[Step 2] 分析 TP 特徵 ...", flush=True)

    # AF bins（0.05 步長）
    af_bins = defaultdict(int)
    af_vals = []
    gq_vals = []
    subtype_counts = Counter()
    verdict_somatic_count = 0
    verdict_subclonal_count = 0
    verdict_germline_count = 0

    # High-AF TP 詳細分析
    highaf_tp = []

    for rec in tp_records:
        af = safe_float(rec['to_af'])
        gq = safe_float(rec['to_gq'])

        if af is not None:
            af_vals.append(af)
            bin_key = f"{int(af * 20) * 5:02d}-{int(af * 20) * 5 + 5:02d}"
            af_bins[bin_key] += 1
            subtype = classify_tp(af)
            subtype_counts[subtype] += 1

            if af >= 0.6:
                highaf_tp.append(rec)

        if gq is not None:
            gq_vals.append(gq)

        if safe_bool(rec.get('to_verdict_somatic')):
            verdict_somatic_count += 1
        if safe_bool(rec.get('to_verdict_subclonal')):
            verdict_subclonal_count += 1
        if safe_bool(rec.get('to_verdict_germline')):
            verdict_germline_count += 1

    n = len(tp_records)
    af_vals.sort()
    gq_vals.sort()

    def median(vals):
        if not vals:
            return None
        m = len(vals) // 2
        return vals[m] if len(vals) % 2 else (vals[m - 1] + vals[m]) / 2

    def pct(vals, p):
        if not vals:
            return None
        idx = int(len(vals) * p)
        return vals[min(idx, len(vals) - 1)]

    print(f"  TP 總數: {n:,}", flush=True)
    print(f"  AF 中位數: {median(af_vals):.3f}, Q1={pct(af_vals, 0.25):.3f}, Q3={pct(af_vals, 0.75):.3f}", flush=True)
    print(f"  亞型: Subclonal={subtype_counts['Subclonal']:,}, Somatic={subtype_counts['Somatic']:,}, High-AF={subtype_counts['High-AF']:,}", flush=True)
    print(f"  Verdict: somatic={verdict_somatic_count:,}, subclonal={verdict_subclonal_count:,}, germline={verdict_germline_count:,}", flush=True)

    return {
        'tp_records': tp_records,
        'af_vals': af_vals,
        'af_bins': af_bins,
        'gq_vals': gq_vals,
        'subtype_counts': subtype_counts,
        'verdict_somatic': verdict_somatic_count,
        'verdict_subclonal': verdict_subclonal_count,
        'verdict_germline': verdict_germline_count,
        'highaf_tp': highaf_tp,
        'median_af': median(af_vals),
        'q1_af': pct(af_vals, 0.25),
        'q3_af': pct(af_vals, 0.75),
        'median_gq': median(gq_vals),
        'n': n,
    }


# ──────────────────────────────────────────────────────────────────
# Step 3: FP 分類分析
# ──────────────────────────────────────────────────────────────────

def analyze_fp(fp_records):
    """FP Germline/LOH/Artifact 分類分析"""
    print("[Step 3] 分析 FP 分類 ...", flush=True)

    classified = defaultdict(list)

    for rec in fp_records:
        fp_type = classify_fp(rec['to_af'], rec['paired_raw_naf'], rec['paired_raw_ndp'])
        classified[fp_type].append(rec)

    # 各組統計
    group_stats = {}
    for group_name, recs in classified.items():
        afs = sorted([f for f in [safe_float(r['to_af']) for r in recs] if f is not None])
        nafs = sorted([f for f in [safe_float(r['paired_raw_naf']) for r in recs] if f is not None])
        gqs = sorted([f for f in [safe_float(r['to_gq']) for r in recs] if f is not None])
        ndps = sorted([f for f in [safe_float(r['paired_raw_ndp']) for r in recs] if f is not None])

        vg = sum(1 for r in recs if safe_bool(r.get('to_verdict_germline')))
        vs = sum(1 for r in recs if safe_bool(r.get('to_verdict_somatic')))
        vsc = sum(1 for r in recs if safe_bool(r.get('to_verdict_subclonal')))

        def med(vals):
            if not vals:
                return None
            m = len(vals) // 2
            return vals[m] if len(vals) % 2 else (vals[m - 1] + vals[m]) / 2

        group_stats[group_name] = {
            'n': len(recs),
            'median_af': med(afs),
            'q1_af': afs[len(afs) // 4] if afs else None,
            'q3_af': afs[3 * len(afs) // 4] if afs else None,
            'median_naf': med(nafs),
            'median_gq': med(gqs),
            'median_ndp': med(ndps),
            'verdict_germline_rate': vg / len(recs) if recs else 0,
            'verdict_somatic_rate': vs / len(recs) if recs else 0,
            'verdict_subclonal_rate': vsc / len(recs) if recs else 0,
            'records': recs,
            'af_vals': afs,
        }
        print(f"  {group_name}: N={len(recs):,}, median_af={med(afs) if afs else 'N/A':.3f if afs else ''}, "
              f"median_naf={med(nafs) if nafs else 'N/A':.3f if nafs else ''}, "
              f"vg_rate={vg/len(recs)*100:.1f}%", flush=True)

    return classified, group_stats


# ──────────────────────────────────────────────────────────────────
# Step 4: ISM 候選池 FP 甲基化特徵差異
# ──────────────────────────────────────────────────────────────────

def analyze_ism_fp_methylation(fp_records, tp_records, allow_unversioned_v1=False):
    """合併 rescue_joined_features.tsv 與 provenance naf，比較甲基化特徵"""
    print("[Step 4] 分析 ISM 候選池 FP 甲基化特徵 ...", flush=True)

    if not os.path.exists(RESCUE_FEATURES):
        print(f"  [警告] rescue_joined_features.tsv 不存在：{RESCUE_FEATURES}", flush=True)
        return None

    # 建立 provenance naf lookup
    prov_fp_lookup = {}
    for rec in fp_records:
        prov_fp_lookup[rec['variant_key']] = rec

    prov_tp_lookup = {}
    for rec in tp_records:
        prov_tp_lookup[rec['variant_key']] = rec

    # 讀取 rescue features
    ism_fp = []
    ism_tp = []
    unmatched = 0

    with open(RESCUE_FEATURES, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        source_rows = list(reader)
        source_columns = reader.fieldnames
        source_rows, verification_metadata = attach_historical_legacy_records(
            source_rows,
            source_columns,
            allow_unversioned_v1,
        )
        for row in source_rows:
            region_key = row.get('region_key', '')
            # region_key 格式: chrom:pos:ref:alt 或 chr:pos:ref>alt
            # variant_key 格式可能略有不同，需嘗試匹配
            truth_status = row.get('truth_status', '')

            # 嘗試多種 key 格式
            candidate_keys = [region_key]
            # 如果含 > 則轉換
            if '>' in region_key:
                parts = region_key.split(':')
                if len(parts) == 4:
                    chrom, pos, ref_alt = parts[0], parts[1], parts[2]
                    # 可能是 chr:pos:ref>alt 格式
                    pass
            # 也嘗試 chrom:pos:ref:alt
            sample = row.get('sample', 'HCC1395')
            chrom_field = row.get('chrom', '')
            pos_field = row.get('pos', '')

            # 重建 variant_key (chrom:pos:ref:alt)
            ref = row.get('ref', '')
            alt = row.get('alt', '')
            if chrom_field and pos_field and ref and alt:
                variant_key = f"{chrom_field}:{pos_field}:{ref}:{alt}"
                candidate_keys.append(variant_key)
            # Also try without 'chr' prefix
            if chrom_field.startswith('chr'):
                candidate_keys.append(f"{chrom_field[3:]}:{pos_field}:{ref}:{alt}")

            # 讀取 ISM 特徵
            ism_record = {
                'region_key': region_key,
                'truth_status': truth_status,
                'QS': safe_float(row.get('Quality_Score', '')),
                'PMD': safe_float(row.get('PairwiseMedianDist', '')),
                'CV': safe_float(row.get('CramersV', '')),
                'AD': safe_float(row.get('AlleleDelta', '')),
                'hp_assign': safe_float(row.get('hp_assign_rate', '')),
                'allele_assign': safe_float(row.get('allele_assign_rate', '')),
                'VerificationClass_Legacy_Selected': row['VerificationClass_Legacy_Selected'],
                'verification_selection_field': row['verification_selection_field'],
                'verification_schema_status': row['verification_schema_status'],
                'af': safe_float(row.get('af', '')),
                'gq': safe_float(row.get('gq', '')),
                'naf': None,
                'ndp': None,
                'fp_type': None,
                'prov_matched': False,
            }

            # 嘗試與 provenance 匹配
            matched = False
            for k in candidate_keys:
                if k in prov_fp_lookup:
                    prov = prov_fp_lookup[k]
                    ism_record['naf'] = safe_float(prov['paired_raw_naf'])
                    ism_record['ndp'] = safe_float(prov['paired_raw_ndp'])
                    ism_record['fp_type'] = classify_fp(
                        prov['to_af'], prov['paired_raw_naf'], prov['paired_raw_ndp']
                    )
                    ism_record['prov_matched'] = True
                    matched = True
                    break
                elif k in prov_tp_lookup:
                    ism_record['prov_matched'] = True
                    ism_record['fp_type'] = 'TP_from_prov'
                    matched = True
                    break

            if not matched:
                unmatched += 1

            if truth_status == 'caller_removed_fp':
                ism_fp.append(ism_record)
            elif truth_status == 'caller_lost_tp':
                ism_tp.append(ism_record)

    print(f"  ISM FP={len(ism_fp)}, TP={len(ism_tp)}, unmatched={unmatched}", flush=True)

    # 按 fp_type 分組
    fp_groups = defaultdict(list)
    for rec in ism_fp:
        fp_type = rec['fp_type'] if rec['fp_type'] else 'Uncertain'
        fp_groups[fp_type].append(rec)

    print(f"  FP 分組：", flush=True)
    for k, v in sorted(fp_groups.items(), key=lambda x: -len(x[1])):
        print(f"    {k}: N={len(v)}", flush=True)

    # 特徵比較
    def get_feature_vals(records, feature):
        return [r[feature] for r in records if r.get(feature) is not None]

    def summarize_features(records, label):
        result = {'label': label, 'n': len(records)}
        for feat in ['QS', 'PMD', 'CV', 'AD', 'hp_assign', 'allele_assign', 'af', 'gq']:
            vals = sorted(get_feature_vals(records, feat))
            if vals:
                n = len(vals)
                result[f'median_{feat}'] = vals[n // 2]
                result[f'q1_{feat}'] = vals[n // 4]
                result[f'q3_{feat}'] = vals[3 * n // 4]
            else:
                result[f'median_{feat}'] = None
                result[f'q1_{feat}'] = None
                result[f'q3_{feat}'] = None
        # Explicit historical legacy-four-state distribution.
        vc_counts = Counter(
            r['VerificationClass_Legacy_Selected']
            for r in records
            if r.get('VerificationClass_Legacy_Selected')
        )
        total_vc = sum(vc_counts.values())
        result['noise_pct'] = vc_counts.get('Noise', 0) / total_vc * 100 if total_vc else None
        result['strong_pct'] = vc_counts.get('Strong', 0) / total_vc * 100 if total_vc else None
        result['subclone_pct'] = vc_counts.get('Subclone', 0) / total_vc * 100 if total_vc else None
        result['vc_counts'] = vc_counts
        return result

    tp_summary = summarize_features(ism_tp, 'TP')
    fp_summaries = {}
    for group_name, recs in fp_groups.items():
        fp_summaries[group_name] = summarize_features(recs, group_name)

    # Mann-Whitney U test vs TP
    def mw_test(fp_vals, tp_vals):
        if len(fp_vals) < 3 or len(tp_vals) < 3:
            return None
        try:
            stat, p = mannwhitneyu(fp_vals, tp_vals, alternative='two-sided')
            return p
        except Exception:
            return None

    print("\n  [ISM 特徵比較表] QS, PMD, CramersV（TP vs FP groups）:", flush=True)
    header = f"  {'Group':<30} | {'N':>5} | {'med_QS':>7} | {'med_PMD':>7} | {'med_CV':>7} | {'Noise%':>7} | {'Strong%':>7}"
    print(header, flush=True)
    print("  " + "-" * 80, flush=True)
    print(f"  {'TP':<30} | {tp_summary['n']:>5} | "
          f"{tp_summary['median_QS']:>7.1f} | {tp_summary['median_PMD']:>7.3f} | "
          f"{tp_summary['median_CV']:>7.3f} | {tp_summary['noise_pct']:>6.1f}% | "
          f"{tp_summary['strong_pct']:>6.1f}%", flush=True)
    for gname, fs in sorted(fp_summaries.items(), key=lambda x: -x[1]['n']):
        qs_p = mw_test(
            get_feature_vals(fp_groups[gname], 'QS'),
            get_feature_vals(ism_tp, 'QS')
        )
        cv_p = mw_test(
            get_feature_vals(fp_groups[gname], 'CV'),
            get_feature_vals(ism_tp, 'CV')
        )
        qs_str = f"{fs['median_QS']:>7.1f}" if fs['median_QS'] is not None else "    N/A"
        pmd_str = f"{fs['median_PMD']:>7.3f}" if fs['median_PMD'] is not None else "    N/A"
        cv_str = f"{fs['median_CV']:>7.3f}" if fs['median_CV'] is not None else "    N/A"
        noise_str = f"{fs['noise_pct']:>6.1f}%" if fs['noise_pct'] is not None else "    N/A"
        strong_str = f"{fs['strong_pct']:>6.1f}%" if fs['strong_pct'] is not None else "    N/A"
        p_info = ""
        if qs_p is not None:
            p_info += f" QS_p={qs_p:.3f}"
        if cv_p is not None:
            p_info += f" CV_p={cv_p:.3f}"
        print(f"  {gname:<30} | {fs['n']:>5} | {qs_str} | {pmd_str} | {cv_str} | {noise_str} | {strong_str} |{p_info}", flush=True)

    return {
        'ism_tp': ism_tp,
        'ism_fp': ism_fp,
        'fp_groups': dict(fp_groups),
        'tp_summary': tp_summary,
        'fp_summaries': fp_summaries,
        'verification_metadata': verification_metadata,
    }


# ──────────────────────────────────────────────────────────────────
# Step 5: 高 AF TP 的 LOH 特徵探索
# ──────────────────────────────────────────────────────────────────

def analyze_highaf_tp(tp_analysis, ism_data):
    """高 AF TP 的 LOH 特徵"""
    print("[Step 5] 高 AF TP LOH 特徵探索 ...", flush=True)

    highaf_tp = tp_analysis['highaf_tp']
    n_highaf = len(highaf_tp)
    print(f"  高 AF TP (AF >= 0.6): N={n_highaf:,} ({n_highaf/tp_analysis['n']*100:.1f}%)", flush=True)

    # 按 naf 分類
    loh_confirmed = []  # naf ≈ 0 AND paired_raw_af >= 0.6 → LOH
    germline_in_cn = []  # naf > 0 → germline in CN-altered region
    uncertain_highaf = []  # ndp = 0

    for rec in highaf_tp:
        naf = safe_float(rec['paired_raw_naf'])
        ndp = safe_float(rec['paired_raw_ndp'])
        paired_af = safe_float(rec['paired_raw_af'])

        if ndp is None or ndp == 0:
            uncertain_highaf.append(rec)
        elif naf is not None and naf < 0.05:
            if paired_af is not None and paired_af >= 0.6:
                loh_confirmed.append(rec)  # high AF in both TO and paired tumor
            else:
                loh_confirmed.append(rec)  # naf=0 anyway
        elif naf is not None and naf >= 0.2:
            germline_in_cn.append(rec)
        else:
            uncertain_highaf.append(rec)

    print(f"  高 AF TP 分類:", flush=True)
    print(f"    LOH-confirmed (naf≈0): N={len(loh_confirmed):,} ({len(loh_confirmed)/n_highaf*100:.1f}%)", flush=True)
    print(f"    Germline-in-CN (naf>0.2): N={len(germline_in_cn):,} ({len(germline_in_cn)/n_highaf*100:.1f}%)", flush=True)
    print(f"    Uncertain: N={len(uncertain_highaf):,} ({len(uncertain_highaf)/n_highaf*100:.1f}%)", flush=True)

    # ISM 特徵子集（高 AF TP 中有 ISM 特徵的）
    if ism_data:
        ism_tp = ism_data['ism_tp']
        highaf_ism = [r for r in ism_tp if r.get('af') is not None and r['af'] >= 0.6]
        print(f"  高 AF TP 在 ISM 候選池中: N={len(highaf_ism)}", flush=True)
        if highaf_ism:
            vc_counts = Counter(
                r['VerificationClass_Legacy_Selected']
                for r in highaf_ism
                if r.get('VerificationClass_Legacy_Selected')
            )
            print(f"  ISM legacy VerificationClass 分佈: {dict(vc_counts)}", flush=True)
            cv_vals = sorted([r['CV'] for r in highaf_ism if r.get('CV') is not None])
            if cv_vals:
                print(f"  CramersV: median={cv_vals[len(cv_vals)//2]:.3f}, "
                      f"Q1={cv_vals[len(cv_vals)//4]:.3f}, "
                      f"Q3={cv_vals[3*len(cv_vals)//4]:.3f}", flush=True)

    return {
        'n_highaf': n_highaf,
        'loh_confirmed': loh_confirmed,
        'germline_in_cn': germline_in_cn,
        'uncertain_highaf': uncertain_highaf,
    }


# ──────────────────────────────────────────────────────────────────
# Step 6: 輸出 TSV 統計表
# ──────────────────────────────────────────────────────────────────

def write_outputs(tp_analysis, fp_classified, fp_group_stats, ism_data, highaf_analysis):
    """寫出統計 TSV 檔案"""
    print("[Step 6] 寫出統計表 ...", flush=True)

    # ── tp_af_bins.tsv ──
    af_bin_path = os.path.join(OUTPUT_DIR, "tp_af_bins.tsv")
    with open(af_bin_path, 'w') as f:
        f.write("af_bin_lower\taf_bin_upper\tcount\tfraction\n")
        total = tp_analysis['n']
        for lower in range(0, 100, 5):
            bin_key = f"{lower:02d}-{lower + 5:02d}"
            count = tp_analysis['af_bins'].get(bin_key, 0)
            f.write(f"{lower/100:.2f}\t{(lower+5)/100:.2f}\t{count}\t{count/total:.4f}\n")
    print(f"  寫出: {af_bin_path}", flush=True)

    # ── tp_subtype_summary.tsv ──
    tp_sub_path = os.path.join(OUTPUT_DIR, "tp_subtype_summary.tsv")
    with open(tp_sub_path, 'w') as f:
        f.write("subtype\tcount\tfraction\taf_threshold_description\n")
        total = tp_analysis['n']
        for subtype, count in sorted(tp_analysis['subtype_counts'].items(), key=lambda x: -x[1]):
            desc = {
                'Subclonal': 'AF < 0.2 (亞克隆體細胞突變)',
                'Somatic': '0.2 ≤ AF < 0.6 (主克隆體細胞突變)',
                'High-AF': 'AF ≥ 0.6 (LOH 或 CN gain 區域)',
            }.get(subtype, '')
            f.write(f"{subtype}\t{count}\t{count/total:.4f}\t{desc}\n")
    print(f"  寫出: {tp_sub_path}", flush=True)

    # ── fp_naf_classification.tsv ──
    fp_class_path = os.path.join(OUTPUT_DIR, "fp_naf_classification.tsv")
    with open(fp_class_path, 'w') as f:
        f.write("fp_type\tcount\tfraction\tdescription\n")
        total = sum(len(v) for v in fp_classified.values())
        type_order = ['Germline-escaped', 'LOH-region', 'Mid-AF-artifact',
                      'Tumor-specific', 'Uncertain-naf', 'Uncertain']
        for fp_type in type_order:
            count = len(fp_classified.get(fp_type, []))
            desc = {
                'Germline-escaped': 'paired_raw_naf >= 0.2，逃過 PON 的 germline 變異',
                'LOH-region': 'naf < 0.05 AND AF >= 0.7，LOH 區域高 AF FP',
                'Mid-AF-artifact': 'naf < 0.05 AND 0.4 ≤ AF < 0.7，中等 AF artifact',
                'Tumor-specific': 'naf < 0.05 AND AF < 0.4，低 AF tumor-specific artifact',
                'Uncertain-naf': 'naf 0.05-0.2，介於 germline 和 somatic 之間',
                'Uncertain': 'paired_raw_ndp = 0，無 normal coverage 記錄（主要為 raw_absent）',
            }.get(fp_type, '')
            f.write(f"{fp_type}\t{count}\t{count/total:.4f}\t{desc}\n")
    print(f"  寫出: {fp_class_path}", flush=True)

    # ── fp_group_feature_stats.tsv ──
    fp_stats_path = os.path.join(OUTPUT_DIR, "fp_group_feature_stats.tsv")
    with open(fp_stats_path, 'w') as f:
        f.write("fp_type\tcount\tmedian_af\tq1_af\tq3_af\tmedian_naf\tmedian_gq\t"
                "verdict_germline_rate\tverdict_somatic_rate\tverdict_subclonal_rate\n")
        for fp_type in ['Germline-escaped', 'LOH-region', 'Mid-AF-artifact',
                        'Tumor-specific', 'Uncertain-naf', 'Uncertain']:
            stats = fp_group_stats.get(fp_type, {})
            if not stats:
                f.write(f"{fp_type}\t0\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\n")
                continue
            n = stats['n']
            maf = f"{stats['median_af']:.3f}" if stats['median_af'] is not None else "NA"
            q1af = f"{stats['q1_af']:.3f}" if stats['q1_af'] is not None else "NA"
            q3af = f"{stats['q3_af']:.3f}" if stats['q3_af'] is not None else "NA"
            mnaf = f"{stats['median_naf']:.3f}" if stats['median_naf'] is not None else "NA"
            mgq = f"{stats['median_gq']:.1f}" if stats['median_gq'] is not None else "NA"
            vg = f"{stats['verdict_germline_rate']:.3f}"
            vs = f"{stats['verdict_somatic_rate']:.3f}"
            vsc = f"{stats['verdict_subclonal_rate']:.3f}"
            f.write(f"{fp_type}\t{n}\t{maf}\t{q1af}\t{q3af}\t{mnaf}\t{mgq}\t{vg}\t{vs}\t{vsc}\n")
    print(f"  寫出: {fp_stats_path}", flush=True)

    # ── ism_fp_methylation_comparison.tsv ──
    if ism_data:
        ism_path = os.path.join(OUTPUT_DIR, "ism_fp_methylation_comparison.tsv")
        verification_metadata = ism_data['verification_metadata']
        with open(ism_path, 'w') as f:
            f.write("group\tverification_view\tverification_selection_field\tverification_schema_status\t"
                    "count\tmedian_QS\tq1_QS\tq3_QS\tmedian_PMD\tmedian_CV\t"
                    "median_hp_assign\tnoise_pct\tstrong_pct\tsubclone_pct\n")
            # TP first
            ts = ism_data['tp_summary']
            f.write(
                f"TP\t{verification_metadata['verification_view']}\t"
                f"{verification_metadata['verification_selection_field']}\t"
                f"{verification_metadata['verification_schema_status']}\t{ts['n']}\t"
                f"{ts['median_QS']:.1f}\t{ts['q1_QS']:.1f}\t{ts['q3_QS']:.1f}\t"
                f"{ts['median_PMD']:.3f}\t{ts['median_CV']:.3f}\t"
                f"{ts['median_hp_assign']:.3f}\t"
                f"{ts['noise_pct']:.1f}\t{ts['strong_pct']:.1f}\t{ts['subclone_pct']:.1f}\n"
            )
            for group_name, fs in sorted(ism_data['fp_summaries'].items(), key=lambda x: -x[1]['n']):
                def fmt(v, fmt_str=".3f"):
                    return f"{v:{fmt_str}}" if v is not None else "NA"
                f.write(
                    f"{group_name}\t{verification_metadata['verification_view']}\t"
                    f"{verification_metadata['verification_selection_field']}\t"
                    f"{verification_metadata['verification_schema_status']}\t{fs['n']}\t"
                    f"{fmt(fs['median_QS'], '.1f')}\t{fmt(fs['q1_QS'], '.1f')}\t{fmt(fs['q3_QS'], '.1f')}\t"
                    f"{fmt(fs['median_PMD'])}\t{fmt(fs['median_CV'])}\t"
                    f"{fmt(fs['median_hp_assign'])}\t"
                    f"{fmt(fs['noise_pct'], '.1f')}\t{fmt(fs['strong_pct'], '.1f')}\t{fmt(fs['subclone_pct'], '.1f')}\n"
                )
        print(f"  寫出: {ism_path}", flush=True)

    # ── tp_highaf_loh_summary.tsv ──
    highaf_path = os.path.join(OUTPUT_DIR, "tp_highaf_loh_summary.tsv")
    with open(highaf_path, 'w') as f:
        f.write("highaf_type\tcount\tfraction_of_highaf\tdescription\n")
        n_highaf = highaf_analysis['n_highaf']
        rows = [
            ('LOH-confirmed', len(highaf_analysis['loh_confirmed']),
             'naf ≈ 0，LOH 區域的體細胞突變（alt allele copy gain）'),
            ('Germline-in-CN', len(highaf_analysis['germline_in_cn']),
             'naf > 0.2，CN 改變區域的 germline SNP'),
            ('Uncertain', len(highaf_analysis['uncertain_highaf']),
             'ndp = 0，無 normal coverage 記錄'),
        ]
        for name, count, desc in rows:
            frac = count / n_highaf if n_highaf > 0 else 0
            f.write(f"{name}\t{count}\t{frac:.4f}\t{desc}\n")
    print(f"  寫出: {highaf_path}", flush=True)


# ──────────────────────────────────────────────────────────────────
# 主程式
# ──────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize the fixed historical rescue table's raw VerificationClass "
            "as a legacy-four-state input"
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    print("=" * 70, flush=True)
    print("TO TP 全量特徵畫像 + FP Germline/LOH/Artifact 分類分析", flush=True)
    print("=" * 70, flush=True)
    print(f"Baseline: TP={BASELINE_TP:,}, FP={BASELINE_FP:,}, FN={BASELINE_FN:,}, F1={BASELINE_F1:.6f}", flush=True)
    print()

    # Step 1: 讀取資料
    tp_records, fp_records = load_provenance_data()

    # Step 2: TP 分析
    tp_analysis = analyze_tp(tp_records)
    print()

    # Step 3: FP 分類
    fp_classified, fp_group_stats = analyze_fp(fp_records)
    print()

    # Step 4: ISM 甲基化特徵比較
    ism_data = analyze_ism_fp_methylation(
        fp_records,
        tp_records,
        args.allow_unversioned_v1,
    )
    print()

    # Step 5: 高 AF TP LOH 分析
    highaf_analysis = analyze_highaf_tp(tp_analysis, ism_data)
    print()

    # Step 6: 輸出 TSV
    write_outputs(tp_analysis, fp_classified, fp_group_stats, ism_data, highaf_analysis)

    print()
    print("=" * 70, flush=True)
    print("分析完成！輸出目錄：", OUTPUT_DIR, flush=True)
    print("=" * 70, flush=True)

    # 回傳摘要供報告使用
    return {
        'tp_analysis': tp_analysis,
        'fp_classified': fp_classified,
        'fp_group_stats': fp_group_stats,
        'ism_data': ism_data,
        'highaf_analysis': highaf_analysis,
    }


if __name__ == '__main__':
    main()
