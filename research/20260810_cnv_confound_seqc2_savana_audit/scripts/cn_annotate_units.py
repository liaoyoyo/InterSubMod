#!/usr/bin/env python3
"""Annotate HCC1395 canonical exact-PS units with copy-number state.

Produces the evidence layer for the CN-confound audit: for every frozen unit in
the 2026-07-24 canonical cohort, record what copy-number ground it sits on,
independently from SEQC2 benchmark truth and from the SAVANA WGS fit.

Design decisions that differ from the 2026-07-10 pass (and why)
---------------------------------------------------------------
* Integer CN is read from ngs_benchmark_cnv_{gain,loss}_cn.bed.  Those two
  files are non-overlapping and carry real copy numbers (gain 3-11.5, loss
  0-1); the earlier pass used only the categorical gain/loss/loh file and so
  could not distinguish CN=3 from CN=11.
* LOH is a separate axis.  The categorical file contains 110 overlapping
  records because gain and loh co-occur (GAINLOH).  Collapsing it into one
  label per region loses that; here total-CN and LOH are annotated separately.
* A unit is assigned by base-pair overlap across its whole span, not by its
  midpoint, so units straddling a CN breakpoint are reported as mixed instead
  of being silently forced into one state.
* Both a span-level view (the region) and a site-level view (each sSNV) are
  emitted, because "how much of the region is CN-altered" and "how many sSNVs
  sit on CN-altered ground" are different denominators.

Outputs (all under ../data/):
  hcc1395_unit_cn_annotation.jsonl   one record per canonical unit
  hcc1395_site_cn_annotation.jsonl   one record per active sSNV site
  cn_annotation_summary.json         aggregate counts used by the report
"""

from __future__ import annotations

import bisect
import json
import os
from collections import Counter, defaultdict

SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_gain_cn.bed"
LOSS_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_loss_cn.bed"
CATEGORICAL_BED = f"{SEQC2_DIR}/ngs_benchmark_cnvs_gain_loss_loh.bed"
EXCLUSION_BED = f"{SEQC2_DIR}/exclusion.bed"

SAVANA_CN_TSV = (
    "/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna/"
    "HCC1395_segmented_absolute_copy_number.tsv"
)

CANON = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/HCC1395"
)
TOPOLOGY_JSONL = f"{CANON}/HCC1395.topology.jsonl"
CENSUS_JSONL = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1/HCC1395.census.jsonl"
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

NEUTRAL_CN = 2.0


class IntervalIndex:
    """Per-chromosome sorted intervals with base-pair overlap queries."""

    def __init__(self):
        self._by_chrom = defaultdict(list)
        self._starts = {}
        self._frozen = False

    def add(self, chrom, start, end, payload):
        if self._frozen:
            raise RuntimeError("cannot add after freeze")
        self._by_chrom[chrom].append((int(start), int(end), payload))

    def freeze(self):
        for chrom, items in self._by_chrom.items():
            items.sort(key=lambda t: (t[0], t[1]))
            self._starts[chrom] = [t[0] for t in items]
        self._frozen = True
        return self

    def overlaps(self, chrom, start, end):
        """Return [(overlap_bp, payload), ...] for every intersecting segment."""
        items = self._by_chrom.get(chrom)
        if not items:
            return []
        starts = self._starts[chrom]
        # Segments are sorted by start; scan left from the insertion point far
        # enough to catch long segments that begin well before `start`.
        idx = bisect.bisect_right(starts, end)
        hits = []
        for i in range(idx - 1, -1, -1):
            s, e, payload = items[i]
            if e <= start:
                # Keep scanning: a longer earlier segment may still overlap.
                # Bound the scan by the widest segment seen on this chromosome.
                if s + self._max_len(chrom) < start:
                    break
                continue
            ov = min(end, e) - max(start, s)
            if ov > 0:
                hits.append((ov, payload))
        return hits

    def _max_len(self, chrom):
        key = ("_maxlen", chrom)
        cached = getattr(self, "_maxlen_cache", None)
        if cached is None:
            cached = {}
            self._maxlen_cache = cached
        if key not in cached:
            cached[key] = max((e - s) for s, e, _ in self._by_chrom[chrom])
        return cached[key]

    def point(self, chrom, pos):
        """Return payloads of segments containing `pos` (half-open [s, e))."""
        return [p for _, p in self.overlaps(chrom, pos, pos + 1)]


def load_bed(path, value_fn=None, name=""):
    idx = IntervalIndex()
    n = 0
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            payload = value_fn(parts) if value_fn else True
            idx.add(chrom, start, end, payload)
            n += 1
    idx.freeze()
    print(f"  loaded {name or path}: {n} segments")
    return idx


def load_savana(path):
    """SAVANA absolute CN segments -> IntervalIndex of dicts."""
    idx = IntervalIndex()
    n = 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue

            def num(key):
                raw = p[col[key]] if col.get(key) is not None and col[key] < len(p) else ""
                try:
                    return float(raw)
                except (TypeError, ValueError):
                    return None

            idx.add(
                p[col["chromosome"]],
                int(p[col["start"]]),
                int(p[col["end"]]),
                {
                    "segment_id": p[col["segment_id"]] if "segment_id" in col else None,
                    "total_cn": num("copyNumber"),
                    "minor_cn": num("minorAlleleCopyNumber"),
                    "mean_baf": num("meanBAF"),
                },
            )
            n += 1
    idx.freeze()
    print(f"  loaded SAVANA CN: {n} segments")
    return idx


# ---------------------------------------------------------------------------
# SEQC2 truth resolution
# ---------------------------------------------------------------------------

def seqc2_total_cn_at(gain_idx, loss_idx, chrom, start, end):
    """Base-pair-weighted total CN over [start, end) using SEQC2 integer CN.

    Anything not covered by a gain or loss segment is treated as CN = 2.
    Returns (weighted_cn, state_bp_map, segment_cn_values).
    """
    length = max(end - start, 1)
    state_bp = Counter()
    cn_values = []
    covered = 0
    weighted = 0.0

    for ov, cn in gain_idx.overlaps(chrom, start, end):
        state_bp["gain"] += ov
        covered += ov
        weighted += ov * cn
        cn_values.append(cn)
    for ov, cn in loss_idx.overlaps(chrom, start, end):
        state_bp["loss"] += ov
        covered += ov
        weighted += ov * cn
        cn_values.append(cn)

    neutral_bp = max(length - covered, 0)
    if neutral_bp:
        state_bp["cn2"] += neutral_bp
        weighted += neutral_bp * NEUTRAL_CN

    return weighted / length, dict(state_bp), cn_values


def dominant(state_bp):
    if not state_bp:
        return None, 0.0
    total = sum(state_bp.values())
    state, bp = max(state_bp.items(), key=lambda kv: kv[1])
    return state, bp / total if total else 0.0


def classify(total_cn_state, loh_frac, tol=1e-9):
    """Combine the total-CN axis and the LOH axis into one reported class."""
    has_loh = loh_frac > 0.5
    if total_cn_state == "gain":
        return "gain_loh" if has_loh else "gain"
    if total_cn_state == "loss":
        return "loss_loh" if has_loh else "loss"
    return "cnloh" if has_loh else "neutral"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading CN references...")
    gain_idx = load_bed(GAIN_CN_BED, lambda p: float(p[3]), "SEQC2 gain CN")
    loss_idx = load_bed(LOSS_CN_BED, lambda p: float(p[3]), "SEQC2 loss CN")
    loh_idx = load_bed(
        CATEGORICAL_BED,
        lambda p: p[3],
        "SEQC2 categorical (LOH axis)",
    )
    excl_idx = load_bed(EXCLUSION_BED, lambda p: True, "SEQC2 exclusion")
    savana_idx = load_savana(SAVANA_CN_TSV)

    # LOH-only sub-index: the categorical file mixes gain/loss/loh records.
    loh_only = IntervalIndex()
    for chrom, items in loh_idx._by_chrom.items():
        for s, e, label in items:
            if label == "loh":
                loh_only.add(chrom, s, e, label)
    loh_only.freeze()

    print("Loading canonical census (resolution_class)...")
    resolution = {}
    with open(CENSUS_JSONL) as fh:
        for ln in fh:
            r = json.loads(ln)
            resolution[r["region_id"]] = {
                "resolution_class": r.get("resolution_class"),
                "topology_signature_count": r.get("topology_signature_count"),
                "enumerated_best_tree_count": r.get("enumerated_best_tree_count"),
                "coarse_class_count": r.get("coarse_class_count"),
            }
    print(f"  census rows: {len(resolution)}")

    print("Annotating units...")
    unit_out = open(os.path.join(OUT_DIR, "hcc1395_unit_cn_annotation.jsonl"), "w")
    site_out = open(os.path.join(OUT_DIR, "hcc1395_site_cn_annotation.jsonl"), "w")

    n_units = 0
    n_sites = 0
    n_unlocatable = 0
    unit_class = Counter()
    unit_class_by_status = defaultdict(Counter)
    unit_class_by_resolution = defaultdict(Counter)
    unit_class_by_chrom = defaultdict(Counter)
    unit_class_by_k = defaultdict(Counter)
    site_class = Counter()
    site_class_by_chrom = defaultdict(Counter)
    mixed_units = 0
    excluded_units = 0
    savana_vs_seqc2 = Counter()
    cn_value_hist = Counter()
    af_by_class = defaultdict(list)
    status_counter = Counter()

    with open(TOPOLOGY_JSONL) as fh:
        for ln in fh:
            u = json.loads(ln)
            n_units += 1
            chrom = u.get("chrom")
            positions = u.get("active_positions") or []
            status = u.get("unit_status")
            status_counter[status] += 1
            k = u.get("active_bit_count", 0)

            rec = {
                "unit_id": u.get("unit_id"),
                "region_id": u.get("region_id"),
                "block_id": u.get("block_id"),
                "chrom": chrom,
                "hp_family": u.get("hp_family"),
                "phase_set": u.get("phase_set"),
                "active_bit_count": k,
                "unit_status": status,
                "family_status": u.get("family_status"),
                "read_af_status": u.get("read_af_status"),
                "recurrence_required": u.get("recurrence_required"),
                "best_tree_unique": u.get("best_tree_unique"),
                "total_tree_count": u.get("total_tree_count"),
            }
            rec.update(resolution.get(u.get("region_id"), {}))

            if not positions:
                rec["cn_status"] = "UNLOCATABLE_NO_ACTIVE_POSITION"
                n_unlocatable += 1
                unit_out.write(json.dumps(rec) + "\n")
                continue

            span_start, span_end = min(positions), max(positions) + 1

            wcn, state_bp, cn_vals = seqc2_total_cn_at(
                gain_idx, loss_idx, chrom, span_start, span_end
            )
            dom_state, dom_frac = dominant(state_bp)

            loh_hits = loh_only.overlaps(chrom, span_start, span_end)
            loh_bp = sum(ov for ov, _ in loh_hits)
            loh_frac = loh_bp / max(span_end - span_start, 1)

            excl_bp = sum(ov for ov, _ in excl_idx.overlaps(chrom, span_start, span_end))
            excl_frac = excl_bp / max(span_end - span_start, 1)

            cls = classify(dom_state, loh_frac)
            is_mixed = dom_frac < 1.0 - 1e-9
            if is_mixed:
                mixed_units += 1
            if excl_frac > 0.5:
                excluded_units += 1

            sav_hits = savana_idx.overlaps(chrom, span_start, span_end)
            sav_cn = None
            sav_minor = None
            if sav_hits:
                tot_ov = sum(ov for ov, _ in sav_hits)
                vals = [(ov, p) for ov, p in sav_hits if p["total_cn"] is not None]
                if vals and tot_ov:
                    sav_cn = sum(ov * p["total_cn"] for ov, p in vals) / sum(
                        ov for ov, _ in vals
                    )
                mvals = [(ov, p) for ov, p in sav_hits if p["minor_cn"] is not None]
                if mvals:
                    sav_minor = sum(ov * p["minor_cn"] for ov, p in mvals) / sum(
                        ov for ov, _ in mvals
                    )

            rec.update(
                {
                    "cn_status": "ANNOTATED",
                    "span_start": span_start,
                    "span_end": span_end,
                    "span_bp": span_end - span_start,
                    "n_active_positions": len(positions),
                    "seqc2_total_cn_weighted": round(wcn, 4),
                    "seqc2_state_bp": state_bp,
                    "seqc2_dominant_state": dom_state,
                    "seqc2_dominant_fraction": round(dom_frac, 6),
                    "seqc2_cn_segment_values": cn_vals,
                    "seqc2_loh_fraction": round(loh_frac, 6),
                    "seqc2_class": cls,
                    "spans_cn_breakpoint": is_mixed,
                    "seqc2_exclusion_fraction": round(excl_frac, 6),
                    "savana_total_cn_weighted": (
                        round(sav_cn, 4) if sav_cn is not None else None
                    ),
                    "savana_minor_cn_weighted": (
                        round(sav_minor, 4) if sav_minor is not None else None
                    ),
                }
            )

            unit_class[cls] += 1
            unit_class_by_status[status][cls] += 1
            unit_class_by_chrom[chrom][cls] += 1
            unit_class_by_k[min(k, 6)][cls] += 1
            rc = rec.get("resolution_class")
            if rc:
                unit_class_by_resolution[rc][cls] += 1

            if sav_cn is not None:
                sav_state = (
                    "gain" if sav_cn > 2.5 else "loss" if sav_cn < 1.5 else "cn2"
                )
                seq_state = dom_state
                savana_vs_seqc2[(seq_state, sav_state)] += 1
                rec["savana_state"] = sav_state

            unit_out.write(json.dumps(rec) + "\n")

            # ---- site level ----
            af_map = {}
            for c in u.get("af_coverage") or []:
                pos = c.get("position")
                a, r = c.get("alt_reads"), c.get("ref_reads")
                if pos is not None and a is not None and r is not None and (a + r) > 0:
                    af_map[pos] = (a, r, a / (a + r))

            for pos in positions:
                n_sites += 1
                s_gain = gain_idx.point(chrom, pos)
                s_loss = loss_idx.point(chrom, pos)
                s_loh = bool(loh_only.point(chrom, pos))
                if s_gain:
                    s_state, s_cn = "gain", max(s_gain)
                elif s_loss:
                    s_state, s_cn = "loss", min(s_loss)
                else:
                    s_state, s_cn = "cn2", NEUTRAL_CN
                s_cls = classify(s_state, 1.0 if s_loh else 0.0)
                site_class[s_cls] += 1
                site_class_by_chrom[chrom][s_cls] += 1
                cn_value_hist[s_cn] += 1

                af = af_map.get(pos)
                if af:
                    af_by_class[s_cls].append(af[2])

                site_out.write(
                    json.dumps(
                        {
                            "unit_id": u.get("unit_id"),
                            "region_id": u.get("region_id"),
                            "chrom": chrom,
                            "position": pos,
                            "hp_family": u.get("hp_family"),
                            "unit_status": status,
                            "seqc2_total_cn": s_cn,
                            "seqc2_state": s_state,
                            "seqc2_loh": s_loh,
                            "seqc2_class": s_cls,
                            "alt_reads": af[0] if af else None,
                            "ref_reads": af[1] if af else None,
                            "read_af": round(af[2], 6) if af else None,
                        }
                    )
                    + "\n"
                )

    unit_out.close()
    site_out.close()

    def pct(n, d):
        return round(100.0 * n / d, 4) if d else None

    def af_stats(vals):
        if not vals:
            return None
        v = sorted(vals)
        n = len(v)
        return {
            "n": n,
            "median": round(v[n // 2], 4),
            "q1": round(v[n // 4], 4),
            "q3": round(v[(3 * n) // 4], 4),
            "mean": round(sum(v) / n, 4),
        }

    located = n_units - n_unlocatable
    altered = sum(v for k, v in unit_class.items() if k != "neutral")
    site_altered = sum(v for k, v in site_class.items() if k != "neutral")

    summary = {
        "generated_by": os.path.basename(__file__),
        "inputs": {
            "topology_jsonl": TOPOLOGY_JSONL,
            "census_jsonl": CENSUS_JSONL,
            "seqc2_gain_cn_bed": GAIN_CN_BED,
            "seqc2_loss_cn_bed": LOSS_CN_BED,
            "seqc2_categorical_bed": CATEGORICAL_BED,
            "seqc2_exclusion_bed": EXCLUSION_BED,
            "savana_cn_tsv": SAVANA_CN_TSV,
        },
        "unit_totals": {
            "units_in_topology_jsonl": n_units,
            "units_annotated": located,
            "units_unlocatable": n_unlocatable,
            "units_spanning_cn_breakpoint": mixed_units,
            "units_majority_in_exclusion": excluded_units,
        },
        "unit_status_counts": dict(status_counter),
        "unit_class_counts": dict(unit_class),
        "unit_class_percent": {k: pct(v, located) for k, v in unit_class.items()},
        "unit_cn_altered": {
            "count": altered,
            "denominator": located,
            "percent": pct(altered, located),
        },
        "site_totals": {"sites": n_sites},
        "site_class_counts": dict(site_class),
        "site_class_percent": {k: pct(v, n_sites) for k, v in site_class.items()},
        "site_cn_altered": {
            "count": site_altered,
            "denominator": n_sites,
            "percent": pct(site_altered, n_sites),
        },
        "unit_class_by_status": {k: dict(v) for k, v in unit_class_by_status.items()},
        "unit_class_by_resolution_class": {
            k: dict(v) for k, v in unit_class_by_resolution.items()
        },
        "unit_class_by_chrom": {k: dict(v) for k, v in unit_class_by_chrom.items()},
        "unit_class_by_active_bit_count": {
            str(k): dict(v) for k, v in unit_class_by_k.items()
        },
        "site_class_by_chrom": {k: dict(v) for k, v in site_class_by_chrom.items()},
        "seqc2_total_cn_histogram_sites": {
            str(k): v for k, v in sorted(cn_value_hist.items())
        },
        "read_af_by_cn_class": {k: af_stats(v) for k, v in af_by_class.items()},
        "savana_vs_seqc2_unit_confusion": {
            f"seqc2={a}|savana={b}": c for (a, b), c in sorted(savana_vs_seqc2.items())
        },
    }

    out_path = os.path.join(OUT_DIR, "cn_annotation_summary.json")
    with open(out_path, "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(f"\nunits={n_units} annotated={located} unlocatable={n_unlocatable}")
    print(f"unit CN-altered: {altered}/{located} = {pct(altered, located)}%")
    print(f"site CN-altered: {site_altered}/{n_sites} = {pct(site_altered, n_sites)}%")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
