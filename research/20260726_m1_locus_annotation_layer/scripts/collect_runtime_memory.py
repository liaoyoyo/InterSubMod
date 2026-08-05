#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Collect measured runtime and memory for the C++ lineage stages.

Two independent measurements, kept separate because they mean different things:

  solver   per-unit `solver_elapsed_microseconds` and `search_nodes` recorded by the
           exact-PS topology binary itself -> the cost of solving one region's lineage
  process  /usr/bin/time -v around the signature-census binary -> wall clock, CPU and
           peak RSS for a whole-sample pass

Per-unit solver time excludes I/O and JSON serialisation, so it must never be presented
as the end-to-end cost of a sample; both are reported side by side for that reason.
"""
import argparse
import json
import os
import re
import statistics
import sys

TOPOLOGY = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples")
CENSUS = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
          "20260724_exact_ps_cpp_topology_signature_census/all7_v1")
SAMPLES = ("HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829")


def parse_time_v(path):
    """read GNU time -v output; returns None when the file is absent"""
    if not os.path.exists(path):
        return None
    out = {}
    for line in open(path):
        line = line.strip()
        m = re.match(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)$", line)
        if m:
            parts = [float(x) for x in m.group(1).split(":")]
            secs = 0.0
            for part in parts:
                secs = secs * 60 + part
            out["wall_seconds"] = round(secs, 3)
        for key, label in (("User time (seconds)", "user_seconds"),
                           ("System time (seconds)", "sys_seconds"),
                           ("Maximum resident set size (kbytes)", "max_rss_kb")):
            if line.startswith(key):
                out[label] = float(line.split(":")[-1].strip())
        if line.startswith("Percent of CPU this job got"):
            out["cpu_percent"] = line.split(":")[-1].strip()
        if line.startswith("Exit status"):
            out["exit_status"] = int(line.split(":")[-1].strip())
    return out or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    per_sample = {}
    cohort_us = 0
    cohort_units = 0
    cohort_nodes = 0
    slowest = {"micros": -1}
    for sample in SAMPLES:
        topo = os.path.join(TOPOLOGY, sample, f"{sample}.topology.jsonl")
        if not os.path.exists(topo):
            sys.exit(f"FAIL CLOSED: missing {topo}")
        micros = []
        nodes = []
        with open(topo) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                us = rec.get("solver_elapsed_microseconds")
                if us is None:
                    continue
                us = int(us)
                micros.append(us)
                nodes.append(int(rec.get("search_nodes") or 0))
                if us > slowest["micros"]:
                    slowest = {"micros": us, "sample": sample,
                               "region_id": rec.get("region_id"),
                               "active_bit_count": rec.get("active_bit_count"),
                               "search_nodes": rec.get("search_nodes"),
                               "unit_status": rec.get("unit_status"),
                               "family_status": rec.get("family_status")}
        if not micros:
            sys.exit(f"FAIL CLOSED: {sample} topology has no solver timings")
        micros.sort()
        total_us = sum(micros)
        cohort_us += total_us
        cohort_units += len(micros)
        cohort_nodes += sum(nodes)
        proc = parse_time_v(os.path.join(CENSUS, f"{sample}.time.txt"))
        per_sample[sample] = {
            "units": len(micros),
            "solver_total_seconds": round(total_us / 1e6, 3),
            "solver_mean_us": round(total_us / len(micros), 1),
            "solver_median_us": micros[len(micros) // 2],
            "solver_p95_us": micros[int(len(micros) * 0.95)],
            "solver_p99_us": micros[int(len(micros) * 0.99)],
            "solver_max_us": micros[-1],
            "search_nodes_total": sum(nodes),
            "search_nodes_max": max(nodes) if nodes else 0,
            "census_process": proc,
        }

    census_rss = [v["census_process"]["max_rss_kb"] for v in per_sample.values()
                  if v["census_process"] and "max_rss_kb" in v["census_process"]]
    census_wall = [v["census_process"]["wall_seconds"] for v in per_sample.values()
                   if v["census_process"] and "wall_seconds" in v["census_process"]]

    out = {
        "schema_name": "intersubmod.runtime_memory.v1",
        "schema_version": "1.0.0",
        "solver_measurement": "per-unit solver_elapsed_microseconds from the exact-PS "
                              "topology binary; excludes I/O and JSON serialisation",
        "process_measurement": "/usr/bin/time -v around the signature-census binary; "
                               "whole-sample pass including I/O",
        "per_sample": per_sample,
        "cohort": {
            "units_with_timing": cohort_units,
            "solver_total_seconds": round(cohort_us / 1e6, 3),
            "solver_mean_us": round(cohort_us / cohort_units, 1),
            "search_nodes_total": cohort_nodes,
            "census_wall_seconds_total": round(sum(census_wall), 3) if census_wall else None,
            "census_wall_seconds_max": max(census_wall) if census_wall else None,
            "census_max_rss_kb": max(census_rss) if census_rss else None,
            "census_max_rss_mb": round(max(census_rss) / 1024, 1) if census_rss else None,
        },
        "slowest_unit": slowest,
        "interpretation_guard": (
            "solver_total_seconds is the sum of per-region solve time on one core; it is "
            "not wall-clock for the pipeline, which also pays BAM/JSON I/O. Peak RSS is "
            "measured on the census binary, the most memory-hungry C++ stage in this run."),
    }
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
