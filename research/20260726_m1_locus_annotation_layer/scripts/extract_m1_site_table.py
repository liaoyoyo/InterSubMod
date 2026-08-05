#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stream the authority M1 stable-multigroup assignments into a compact per-site table.

Read-only against the authority artifact. Emits one row per M1 site with:
  identity      : dataset, chrom, pos, ref, alt, truth_label
  methyl screen : K (distinct coarse labels), group_sizes, core_read_n, min_group_n
  HP presence   : distinct latest_hp families at the locus -> HP1_ONLY / HP2_ONLY / BOTH / OTHER
  PS            : distinct latest_ps values
  association   : hp_exact / hp_family Cramer V + permutation p + aligned flag

Fail-closed: any record missing dataset/chrom/pos aborts. No value is defaulted or invented.
"""
import argparse
import collections
import gzip
import hashlib
import json
import os
import sys

AUTHORITY = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/"
    "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
)

COLUMNS = [
    "dataset", "chrom", "pos", "ref", "alt", "truth_label",
    "methyl_K", "group_sizes", "core_read_n", "min_group_n", "min_group_share",
    "hp_presence", "hp_families", "n_ps", "ps_values",
    "assoc_hp_exact_v", "assoc_hp_exact_p", "assoc_hp_exact_aligned",
    "assoc_hp_family_v", "assoc_hp_family_p", "assoc_hp_family_aligned",
]


def hp_family(tag):
    """Map a raw latest_hp tag to its primary family, matching the exact-PS contract.

    HP1 <- {1, 1-1, 1-2};  HP2 <- {2, 2-1, 2-2};  anything else is kept verbatim.
    """
    if tag is None:
        return "NA"
    t = str(tag).strip()
    if t in ("", "None", "nan"):
        return "NA"
    head = t.split("-", 1)[0]
    if head == "1":
        return "HP1"
    if head == "2":
        return "HP2"
    return f"OTHER({head})"


def assoc_fields(assoc, key):
    node = (assoc or {}).get(key) or {}
    return node.get("v"), node.get("p_perm"), node.get("aligned")


def fmt(v):
    return "" if v is None else str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=AUTHORITY)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-receipt", required=True)
    ap.add_argument("--limit", type=int, default=0, help="debug: stop after N records")
    a = ap.parse_args()

    if not os.path.exists(a.input):
        sys.exit(f"ERROR: authority artifact missing: {a.input}")

    os.makedirs(os.path.dirname(a.out_tsv), exist_ok=True)

    n = 0
    skipped_no_labels = 0
    per_dataset = collections.Counter()
    per_hp_presence = collections.Counter()
    per_K = collections.Counter()
    contracts = collections.Counter()

    with gzip.open(a.input, "rt") as fh, open(a.out_tsv, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)

            dataset = rec.get("dataset")
            chrom = rec.get("chrom")
            pos = rec.get("pos")
            if dataset is None or chrom is None or pos is None:
                sys.exit(f"ERROR: record {n} missing dataset/chrom/pos -> fail closed")

            contracts[rec.get("screen_contract", "MISSING")] += 1
            post = rec.get("posthoc") or {}

            labels = rec.get("coarse_labels_all_after_peel") or rec.get("labels") or []
            if not labels:
                skipped_no_labels += 1
            counts = collections.Counter(labels)
            K = len(counts)
            core_n = len(labels)
            min_n = min(counts.values()) if counts else 0
            min_share = (min_n / core_n) if core_n else 0.0
            group_sizes = ",".join(str(counts[k]) for k in sorted(counts))

            hps = rec.get("latest_hp") or []
            fams = sorted({hp_family(h) for h in hps} - {"NA"})
            if fams == ["HP1"]:
                presence = "HP1_ONLY"
            elif fams == ["HP2"]:
                presence = "HP2_ONLY"
            elif "HP1" in fams and "HP2" in fams:
                presence = "BOTH"
            elif not fams:
                presence = "NO_HP"
            else:
                presence = "OTHER_ONLY"

            ps_vals = sorted({p for p in (rec.get("latest_ps") or []) if p is not None})

            assoc = rec.get("associations") or {}
            ev, ep, ea = assoc_fields(assoc, "hp_exact")
            fv, fp, fa = assoc_fields(assoc, "hp_family")

            row = [
                dataset, chrom, pos, post.get("ref"), post.get("alt"), post.get("truth_label"),
                K, group_sizes, core_n, min_n, f"{min_share:.6f}",
                presence, "|".join(fams), len(ps_vals),
                ",".join(str(p) for p in ps_vals[:8]),
                ev, ep, ea, fv, fp, fa,
            ]
            out.write("\t".join(fmt(x) for x in row) + "\n")

            n += 1
            per_dataset[dataset] += 1
            per_hp_presence[presence] += 1
            per_K[K] += 1
            if a.limit and n >= a.limit:
                break

    sha = hashlib.sha256()
    with open(a.out_tsv, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)

    receipt = {
        "schema_name": "intersubmod.m1_locus_annotation.site_table_receipt",
        "schema_version": "1.0.0",
        "source_authority": a.input,
        "source_bytes": os.path.getsize(a.input),
        "rows_written": n,
        "records_without_labels": skipped_no_labels,
        "screen_contracts_seen": dict(contracts),
        "per_dataset": dict(per_dataset),
        "per_hp_presence": dict(per_hp_presence),
        "per_methyl_K": {str(k): v for k, v in sorted(per_K.items())},
        "out_tsv": a.out_tsv,
        "out_tsv_bytes": os.path.getsize(a.out_tsv),
        "out_tsv_sha256": sha.hexdigest(),
        "columns": COLUMNS,
    }
    with open(a.out_receipt, "w") as fh:
        json.dump(receipt, fh, ensure_ascii=False, indent=1)
    print(json.dumps(receipt, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
