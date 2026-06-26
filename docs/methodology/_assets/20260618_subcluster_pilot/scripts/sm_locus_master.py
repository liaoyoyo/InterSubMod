"""[L0 per-locus 主紀錄] join 所有 sSNV 狀態到單一 SoT 表。
每 sSNV: chrom/pos/src/normal_ref/normal_alt/somatic/tumor_ref/tumor_alt/vaf/HP/
         region/tree_role/cn/seqc2_class/n_partners/n_links。
來源: sm_linkage_genomewide.json(census+pairs) + sm_region_integration.json + CN bed + SEQC2 3 檔。
輸出 報告夾/data/sm_locus_master.{tsv,json}。sum-check: row 數 == union。
"""
import json
import gzip
from bisect import bisect_right
from collections import defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
RPT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report"
S = "/big8_disk/data/HCC1395/SEQC2"
CNBED = f"{S}/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"


def load_intervals(path, lab_col=3):
    iv = defaultdict(list)
    for ln in open(path):
        p = ln.split()
        if len(p) < 3 or not p[1].isdigit():
            continue
        lab = p[lab_col] if len(p) > lab_col else ""
        iv[p[0]].append((int(p[1]), int(p[2]), lab))
    for c in iv:
        iv[c].sort()
    return iv


def in_iv(iv, chrom, pos):
    arr = iv.get(chrom)
    if not arr:
        return None
    i = bisect_right([x[0] for x in arr], pos) - 1
    if i >= 0 and arr[i][1] >= pos:
        return arr[i][2] if arr[i][2] else True
    return None


def load_vcf_set(path):
    s = set()
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.split("\t")
            s.add((p[0], int(p[1])))
    return s


def main():
    d = json.load(open(f"{A}/sm_linkage_genomewide.json"))
    cen = d["census"]
    pairs = d["pairs"]
    regions = json.load(open(f"{A}/sm_region_integration.json"))["regions"]
    cn = load_intervals(CNBED)
    hc = load_vcf_set(f"{S}/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz")
    ss = load_vcf_set(f"{S}/sSNV.MSDUKT.superSet.v1.2.vcf.gz")
    hcreg = load_intervals(f"{S}/High-Confidence_Regions_v1.2.bed", lab_col=99)  # no label col

    # per-sSNV dominant HP from pairs (任一 pair 取得)
    hp = {}
    for p in pairs:
        hp.setdefault(f"{p['chrom']}:{p['a']}", p["hp_a"])
        hp.setdefault(f"{p['chrom']}:{p['b']}", p["hp_b"])

    # region containment + tree_role
    reg_by_chrom = defaultdict(list)
    for r in regions:
        reg_by_chrom[r["chrom"]].append(r)
    for c in reg_by_chrom:
        reg_by_chrom[c].sort(key=lambda r: r["start"])
    # precompute per-region role sets
    for r in regions:
        out = set(u for u, v in r["nested_edges"])
        inc = set(v for u, v in r["nested_edges"])
        sib = set()
        for a, b in r["sibling_pairs"]:
            sib.add(a); sib.add(b)
        r["_out"] = out; r["_inc"] = inc; r["_sib"] = sib

    def find_region_role(chrom, pos):
        arr = reg_by_chrom.get(chrom, [])
        starts = [r["start"] for r in arr]
        i = bisect_right(starts, pos) - 1
        if i < 0 or arr[i]["end"] < pos:
            return None, "isolated"
        r = arr[i]
        key = f"{chrom}:{pos}"
        anc = key in r["_out"]
        des = key in r["_inc"]
        if anc and des:
            role = "intermediate"
        elif anc:
            role = "ancestor"
        elif des:
            role = "descendant"
        elif key in r["_sib"]:
            role = "sibling"
        else:
            role = "region_member"  # co_linked-merged or unlinked-in-region
        return r["region"], role

    def seqc2_class(chrom, pos):
        if (chrom, pos) in hc:
            return "highconf"
        if (chrom, pos) in ss:
            return "superset_only"
        if in_iv(hcreg, chrom, pos):
            return "in_HC_not_called"
        return "outside_HC"

    rows = []
    for key, c in cen.items():
        chrom, pos = c["chrom"], c["pos"]
        region, role = find_region_role(chrom, pos)
        rows.append({
            "locus": key, "chrom": chrom, "pos": pos, "src": c["src"],
            "normal_ref": c["normal_ref"], "normal_alt": c["normal_alt"], "somatic": c["somatic"],
            "tumor_ref": c["tumor_ref"], "tumor_alt": c["tumor_alt"], "vaf": c["vaf"],
            "hp": hp.get(key), "region": region, "tree_role": role,
            "cn": in_iv(cn, chrom, pos) or "neutral", "seqc2_class": seqc2_class(chrom, pos),
            "n_partners": c["n_partners_le50k"], "n_links": c["n_realized_links"],
        })

    # outputs
    cols = ["locus", "chrom", "pos", "src", "normal_ref", "normal_alt", "somatic",
            "tumor_ref", "tumor_alt", "vaf", "hp", "region", "tree_role", "cn",
            "seqc2_class", "n_partners", "n_links"]
    with open(f"{RPT}/data/sm_locus_master.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in cols) + "\n")
    json.dump({"n_loci": len(rows), "columns": cols, "rows": rows},
              open(f"{RPT}/data/sm_locus_master.json", "w"), ensure_ascii=False)

    # summary stats for the record
    from collections import Counter
    role_c = Counter(r["tree_role"] for r in rows)
    role_som = Counter(r["tree_role"] for r in rows if r["somatic"])
    hp_c = Counter(r["hp"] for r in rows if r["somatic"] and r["n_links"] >= 1)
    summ = {"n_loci": len(rows), "n_union": d["universe"]["n_union"],
            "sum_check_ok": len(rows) == d["universe"]["n_union"],
            "tree_role_all": dict(role_c), "tree_role_somatic": dict(role_som),
            "hp_dist_linked_somatic": dict(hp_c),
            "src": dict(Counter(r["src"] for r in rows)),
            "seqc2_class": dict(Counter(r["seqc2_class"] for r in rows)),
            "cn": dict(Counter(r["cn"] for r in rows))}
    json.dump(summ, open(f"{RPT}/data/sm_locus_master_summary.json", "w"), ensure_ascii=False, indent=1)
    print(f"MASTER n_loci={len(rows)} sum_check_ok={summ['sum_check_ok']}")
    print("tree_role (somatic):", dict(role_som))
    print("hp (linked somatic):", dict(hp_c))
    print("seqc2_class:", dict(Counter(r['seqc2_class'] for r in rows)))
    # chr17 spot-check
    print("\nchr17 known loci:")
    for pos in (48357368, 48362515, 48365089, 48365161):
        r = next((x for x in rows if x["locus"] == f"chr17:{pos}"), None)
        if r:
            print(f"  {r['locus']}: src={r['src']} somatic={r['somatic']} vaf={r['vaf']} hp={r['hp']} "
                  f"role={r['tree_role']} cn={r['cn']} seqc2={r['seqc2_class']} region={r['region']}")


if __name__ == "__main__":
    main()
