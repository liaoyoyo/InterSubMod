#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
每區基因註釋(2026-06-29 多樣本任務):region(chrom:start-end) → 基因名/啟動子/癌症相關/oncogene-TSG/用藥。
來源:GENCODE v46 GTF(基因 body + TSS±2kb 啟動子) + DGIdb(可用藥) + COSMIC CGC(oncogene/TSG/tier,optional)。
用法: SM_DATA=<workdir> python3 region_gene_annotation.py   (讀 SM_DATA/topology_per_region.json → 寫 SM_DATA/region_gene_annotation.json)
§13-A: 標籤可溯源(每基因→來源)。
"""
import json, os, gzip, csv
from collections import defaultdict
from bisect import bisect_right

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(HERE, "..", "data")))
ANN = os.environ.get("SM_ANNOT_DIR", "/big7_disk/liaoyoyo2001/gene_annotation")
GTF = f"{ANN}/gencode.v46.basic.annotation.gtf.gz"
DGIDB = f"{ANN}/dgidb_interactions.tsv"
COSMIC = os.environ.get("SM_COSMIC", f"{ANN}/Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz")  # optional; gz/tsv 皆可
PROMOTER_BP = 2000


def parse_gtf():
    """回 per-chrom sorted [(start,end,strand,name,gtype)] + tss list。"""
    genes = defaultdict(list)
    with gzip.open(GTF, "rt") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.split("\t")
            if p[2] != "gene":
                continue
            attr = p[8]
            name = _attr(attr, "gene_name"); gtype = _attr(attr, "gene_type")
            genes[p[0]].append((int(p[3]), int(p[4]), p[6], name, gtype))
    for c in genes:
        genes[c].sort()
    return genes


def _attr(s, key):
    i = s.find(key + ' "')
    if i < 0:
        return ""
    j = s.find('"', i + len(key) + 2)
    return s[i + len(key) + 2:j]


def load_dgidb():
    """gene_name → set(drug)。"""
    d = defaultdict(set)
    if not os.path.exists(DGIDB):
        return d
    with open(DGIDB) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = (r.get("gene_name") or "").strip().upper()
            dr = (r.get("drug_name") or r.get("drug_claim_name") or "").strip()
            if g and dr:
                d[g].add(dr)
    return d


def load_cosmic():
    """gene → {role, tier, tumour}。COSMIC CGC v104(欄 GENE_SYMBOL/ROLE_IN_CANCER/TIER/TUMOUR_TYPES_SOMATIC,tab-sep,可 gz)。optional。"""
    c = {}
    if not os.path.exists(COSMIC):
        return c
    op = gzip.open(COSMIC, "rt") if COSMIC.endswith(".gz") else open(COSMIC)
    with op as f:
        rd = csv.DictReader(f, delimiter="\t")  # CGC v104 = TSV
        for r in rd:
            g = (r.get("GENE_SYMBOL") or r.get("Gene Symbol") or r.get("gene_symbol") or "").strip().upper()
            if not g:
                continue
            c[g] = {"role": (r.get("ROLE_IN_CANCER") or r.get("Role in Cancer") or "").strip(),
                    "tier": (r.get("TIER") or r.get("Tier") or "").strip(),
                    "tumour": (r.get("TUMOUR_TYPES_SOMATIC") or r.get("Tumour Types(Somatic)") or "").strip()[:80]}
    return c


def overlap_genes(genes_c, s, e):
    """回 body-overlap genes + promoter-overlap genes。genes_c=sorted by start。"""
    body, promo = [], []
    for gs, ge, strand, name, gtype in genes_c:
        if gs > e:
            break
        if ge < s:
            continue
        if gs <= e and ge >= s:  # body overlap
            body.append((name, gtype))
        # promoter = TSS±PROMOTER_BP
        tss = gs if strand == "+" else ge
        if (tss - PROMOTER_BP) <= e and (tss + PROMOTER_BP) >= s:
            promo.append(name)
    return body, promo


def main():
    print("loading GENCODE...", flush=True)
    genes = parse_gtf()
    dgidb = load_dgidb()
    cosmic = load_cosmic()
    print(f"genes chroms={len(genes)} dgidb={len(dgidb)} cosmic={len(cosmic)}", flush=True)
    det = json.load(open(os.path.join(DATA, "topology_per_region.json")))["detail"]
    out = {}
    nhit_cancer = nhit_drug = nhit_promo = 0
    for r in det:
        chrom = r["chrom"]; s = r["start"]; e = int(r["region"].split("-")[-1])
        body, promo = overlap_genes(genes.get(chrom, []), s, e)
        names = sorted(set(n for n, t in body if n))
        protein = sorted(set(n for n, t in body if t == "protein_coding" and n))
        cancer = {n: cosmic[n.upper()] for n in names if n.upper() in cosmic}
        drug = {n: sorted(dgidb[n.upper()])[:8] for n in names if n.upper() in dgidb}
        rec = {"genes": names[:30], "n_genes": len(names), "protein_coding": protein[:30],
               "promoter_genes": sorted(set(promo))[:20], "has_promoter": bool(promo),
               "cancer_genes": cancer, "druggable_genes": drug}
        out[r["region"]] = rec
        if cancer:
            nhit_cancer += 1
        if drug:
            nhit_drug += 1
        if promo:
            nhit_promo += 1
    summary = {"n_regions": len(det), "regions_with_cancer_gene": nhit_cancer,
               "regions_with_druggable": nhit_drug, "regions_with_promoter": nhit_promo,
               "sources": {"gencode": os.path.basename(GTF), "dgidb": os.path.basename(DGIDB),
                           "cosmic": os.path.basename(COSMIC) if cosmic else "(none)"},
               "promoter_def": f"TSS±{PROMOTER_BP}bp"}
    json.dump({"summary": summary, "regions": out}, open(os.path.join(DATA, "region_gene_annotation.json"), "w"), ensure_ascii=False)
    print("GENE ANNOTATION DONE", json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
