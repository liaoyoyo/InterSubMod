#!/usr/bin/env python3
"""[P2 sSNV read-level 連鎖驗證] 對 multi-sSNV-linkable 候選: pysam 抽 tumor BAM per-read 在各 sSNV 的等位(REF/ALT)
+ HP tag → 共現 2×2(每 sSNV 對) + 互斥檢定 + HP 一致 + normal BAM REF 確認(somatic)。chr2:18M 式非循環 lineage 驗證。
輸出 p2_linkage.json。用法: p2_snv_linkage.py [chr17:48360161 ...] (無參數=全 linkable)。"""
import json, gzip, sys
from collections import Counter, defaultdict
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
VD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"


def load_vcf_refalt(chrom):
    d = {}
    with gzip.open(f"{VD}/filtered_snv_tp_{chrom}.vcf.gz", "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.split("\t")
            d[int(p[1])] = (p[3], p[4])
    return d


def per_read_allele(bam, chrom, snvs):
    """snvs=[(pos,ref,alt)]. 回 {read_name: {pos: 'REF'/'ALT'/'OTHER'}} + {read_name: hp}。"""
    calls = defaultdict(dict); hp = {}
    for pos, ref, alt in snvs:
        for col in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0, stepper="samtools"):
            if col.reference_pos != pos - 1:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                aln = pr.alignment
                if aln.mapping_quality < 20:
                    continue
                base = aln.query_sequence[pr.query_position].upper()
                rn = aln.query_name
                calls[rn][pos] = "REF" if base == ref else ("ALT" if base == alt else "OTHER")
                if aln.has_tag("HP"):
                    hp[rn] = aln.get_tag("HP")
    return calls, hp


def analyze(chrom, anchor, near, refalt):
    snvs = sorted(set([anchor] + [p for p in near]))
    snvs = [(p, refalt[p][0], refalt[p][1]) for p in snvs if p in refalt]
    if len(snvs) < 2:
        return {"status": "LT2_SNV_IN_VCF", "n_snv": len(snvs)}
    tb = pysam.AlignmentFile(TBAM, "rb")
    nb = pysam.AlignmentFile(NBAM, "rb")
    tcalls, thp = per_read_allele(tb, chrom, snvs)
    ncalls, _ = per_read_allele(nb, chrom, snvs)
    tb.close(); nb.close()
    poss = [p for p, _, _ in snvs]
    # 覆蓋 ≥2 sSNV 的 tumor reads
    multi = {rn: c for rn, c in tcalls.items() if sum(1 for p in poss if p in c) >= 2}
    # 每 sSNV 對的 2×2 共現 + 關係分類
    def classify(rr, ra, ar, aa):
        if aa < 2:
            return "mutual_excl" if (ra + ar) >= 2 else "sparse"   # 互斥=sibling lineage
        if ra == 0 and ar == 0:
            return "co_linked"          # 完美共連 = 同 haplotype/lineage
        if ar == 0:
            return "nested_a_in_b"      # a-ALT ⊂ b-ALT (b 祖先, a 衍生)
        if ra == 0:
            return "nested_b_in_a"      # b-ALT ⊂ a-ALT (a 祖先, b 衍生)
        return "independent"            # 無乾淨結構
    pairs = {}
    for i in range(len(poss)):
        for j in range(i + 1, len(poss)):
            a, b = poss[i], poss[j]
            tab = Counter()
            for rn, c in multi.items():
                if a in c and b in c and c[a] in ("REF", "ALT") and c[b] in ("REF", "ALT"):
                    tab[(c[a], c[b])] += 1
            n = sum(tab.values())
            if n >= 4:
                rr, ra, ar, aa = tab[("REF", "REF")], tab[("REF", "ALT")], tab[("ALT", "REF")], tab[("ALT", "ALT")]
                pairs[f"{a}_{b}"] = {"n_coread": n, "REF_REF": rr, "REF_ALT": ra, "ALT_REF": ar, "ALT_ALT": aa,
                                     "rel": classify(rr, ra, ar, aa)}
    # 每 sSNV 的 tumor ALT 數 + normal REF 確認 (somatic)
    snv_stat = {}
    for p, ref, alt in snvs:
        tc = Counter(c.get(p) for c in tcalls.values() if p in c)
        nc = Counter(c.get(p) for c in ncalls.values() if p in c)
        n_norm = sum(nc.values())
        snv_stat[p] = {"ref": ref, "alt": alt, "tumor_REF": tc["REF"], "tumor_ALT": tc["ALT"],
                       "normal_REF": nc["REF"], "normal_ALT": nc["ALT"], "n_normal": n_norm,
                       "somatic_confirmed": n_norm >= 4 and nc["ALT"] == 0}
    # HP 一致: ALT reads 是否共享 HP
    alt_hp = defaultdict(Counter)
    for rn, c in tcalls.items():
        for p in poss:
            if c.get(p) == "ALT" and rn in thp:
                alt_hp[p][thp[rn]] += 1
    # 🔴 HP 一致性: 互斥/嵌套只在「同 germline HP」才是 subclone; 不同 HP = allelic(等位非 subclone)
    def dom_hp(pos):
        h = alt_hp.get(pos, {})
        return max(h, key=h.get) if h else None
    subclone_pairs = 0; allelic_pairs = 0
    for pk, v in pairs.items():
        if v["rel"] in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl"):
            a, b = pk.split("_"); ha, hb = dom_hp(int(a)), dom_hp(int(b))
            v["hp_a"] = ha; v["hp_b"] = hb; v["same_hp"] = (ha == hb and ha is not None)
            if v["same_hp"]:
                subclone_pairs += 1
            else:
                allelic_pairs += 1
    rels = [v["rel"] for v in pairs.values()]
    n_structured = sum(1 for r in rels if r in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl"))
    n_som = sum(1 for v in snv_stat.values() if v["somatic_confirmed"])
    has_excl = any(v["rel"] == "mutual_excl" and v.get("same_hp") for v in pairs.values())
    has_nest = any(v["rel"] in ("nested_a_in_b", "nested_b_in_a") and v.get("same_hp") for v in pairs.values())
    if subclone_pairs >= 1 and n_som >= 2:
        lv = "subclone_confirmed_sibling" if has_excl else ("subclone_confirmed_nested" if has_nest else "subclone_confirmed_linked")
    elif subclone_pairs >= 1:
        lv = "subclone_structure_but_<2_somatic"
    elif allelic_pairs >= 1:
        lv = "allelic_not_subclone"
    else:
        lv = "no_clean_linkage"
    return {"status": "OK", "n_snv": len(snvs), "snv_positions": poss, "n_multi_read": len(multi),
            "n_somatic_confirmed": n_som, "n_structured_pairs": n_structured, "rels": dict(Counter(rels)),
            "lineage_verdict": lv, "pairs": pairs, "snv_stat": snv_stat, "alt_hp": {p: dict(v) for p, v in alt_hp.items()}}


def main():
    feas = json.load(open(f"{A}/p2_feasibility.json"))
    per = {r["key"]: r for r in feas["per_candidate"]}
    if len(sys.argv) > 1:
        keys = sys.argv[1:]
    else:
        link = feas["summary"]["linkable_9"] + feas["summary"]["linkable_30"]
        keys = sorted({l["key"] for l in link})
    print(f"linkage loci: {len(keys)} → {keys}", flush=True)
    out = {}
    vcache = {}
    for key in keys:
        r = per.get(key)
        if not r:
            continue
        chrom = key.split(":")[0]; anchor = r["anchor"]
        near = [p for p in r["near_snvs"] if abs(p - anchor) <= 10000]
        if chrom not in vcache:
            vcache[chrom] = load_vcf_refalt(chrom)
        res = analyze(chrom, anchor, near, vcache[chrom])
        res["set"] = r["set"]; res["nact"] = r["nact"]; res["in9"] = r["in9"]
        out[key] = res
        st = res.get("status")
        if st == "OK":
            print(f"  {key} [{r['set']},in9={r['in9']}]: {res['n_snv']} sSNV, {res['n_multi_read']} multi-read, "
                  f"rels={res['rels']}, {res['n_somatic_confirmed']}/{res['n_snv']} somatic → {res['lineage_verdict']}", flush=True)
        else:
            print(f"  {key}: {st}", flush=True)
    json.dump(out, open(f"{A}/p2_linkage.json", "w"), ensure_ascii=False, indent=1)
    print("[-> p2_linkage.json]", flush=True)


if __name__ == "__main__":
    main()
