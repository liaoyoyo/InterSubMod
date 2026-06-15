#!/usr/bin/env python3
"""
Analyze extracted matrix JSON. Pure read-back + tabulation (no fabrication).
Outputs structured verification tables for claim groups B/C/D.
Usage: analyze_matrix.py <matrix_json>
"""
import sys, json
from collections import Counter, defaultdict

j = json.load(open(sys.argv[1]))
SNV = j["snv_1based"]
CPG = j["cpg_1based"]
ORDER = ["1", "2", "3", "4", "5", "6"]
CPG_ORDER = ["2.1", "2.2", "3.1", "3.2", "3.3", "3.4", "3.5", "4.1", "5.1", "5.2"]

print(f"### {j['label']}  region={j['region']}")
print(f"tumor_reads={len(j['tumor_reads'])}  normal_reads={len(j['normal_reads'])}\n")


def allele_counts(reads, label_set):
    out = {}
    for lab in label_set:
        c = Counter()
        for r in reads:
            b = r["bases"].get(lab)
            if b is not None:
                c[b] += 1
        out[lab] = c
    return out


def fmt_counts(c):
    tot = sum(c.values())
    return f"n={tot}  " + ", ".join(f"{k}:{v}({100*v/tot:.0f}%)" for k, v in c.most_common())


# ---------- Section 1: per-SNV allele counts ----------
print("=" * 70)
print("S1. PER-SNV ALLELE COUNTS")
print("=" * 70)
tum = j["tumor_reads"]; nor = j["normal_reads"]
tum_hp = defaultdict(list)
for r in tum:
    tum_hp[r["HP"]].append(r)
print(f"\n[TUMOR HP distribution] " + ", ".join(f"HP{k}:{len(v)}" for k, v in sorted(tum_hp.items(), key=lambda x: str(x[0]))))

for lab in ORDER:
    info = SNV[lab]
    print(f"\n--- ({lab}) chr2:{info['pos']} ref={info['ref']} alt={info['alt']} ---")
    ct = allele_counts(tum, [lab])[lab]
    cn = allele_counts(nor, [lab])[lab]
    print(f"  TUMOR  all : {fmt_counts(ct)}")
    for hp in sorted(tum_hp, key=lambda x: str(x)):
        chp = allele_counts(tum_hp[hp], [lab])[lab]
        if sum(chp.values()) > 0:
            print(f"  TUMOR  HP{hp}: {fmt_counts(chp)}")
    print(f"  NORMAL all : {fmt_counts(cn)}")


# ---------- Section 2: linkage reads (cover >=2 SNV) ----------
print("\n" + "=" * 70)
print("S2. LINKAGE READS (cover >=2 of the 6 SNV positions)  [TUMOR]")
print("=" * 70)
def cov_count(r):
    return sum(1 for lab in ORDER if r["bases"].get(lab) is not None)

link = [r for r in tum if cov_count(r) >= 2]
link.sort(key=lambda r: (str(r["HP"]), r["ref_start"]))
print(f"linkage reads (>=2 cov): {len(link)} / {len(tum)}\n")
hdr = "  ".join(f"({l})" for l in ORDER)
print(f"{'HP':>3} {'strand':>3} {'cov':>3}  {hdr}   name")
for r in link:
    bs = []
    for lab in ORDER:
        b = r["bases"].get(lab)
        info = SNV[lab]
        if b is None:
            bs.append(" . ")
        elif b == info["ref"]:
            bs.append(f" {b} ")        # ref
        else:
            bs.append(f"*{b}*")        # non-ref (variant)
    print(f"{str(r['HP']):>3} {r['strand']:>3} {cov_count(r):>3}  " + " ".join(f"{x:>4}" for x in bs) + f"   {r['name'][:16]}")


# ---------- Section 3: methylation at user CpG sites (linkage reads) ----------
print("\n" + "=" * 70)
print("S3. METHYLATION (5mC prob) at user CpG sites  [linkage reads, TUMOR]")
print("=" * 70)
print("H = 5mC prob>=0.5 ; L = <0.5 ; . = no call ; show m_prob")
hdr2 = " ".join(f"{l:>5}" for l in CPG_ORDER)
print(f"{'HP':>3} (3)base  {hdr2}")
for r in link:
    # key discriminator: base at position 3 (A vs G) splits user's groups
    b3 = r["bases"].get("3") or "."
    cells = []
    for cl in CPG_ORDER:
        m = r["meth"].get(cl)
        if m is None or m.get("m_prob") is None:
            cells.append("  .  ")
        else:
            p = m["m_prob"]
            tag = "H" if p >= 0.5 else "L"
            cells.append(f"{tag}{p:.2f}")
    print(f"{str(r['HP']):>3} {b3:>6}  " + " ".join(f"{c:>5}" for c in cells))


# ---------- Section 4: read group signature (pattern over 6 positions) ----------
print("\n" + "=" * 70)
print("S4. READ GROUP SIGNATURES (variant pattern over covered SNVs)")
print("=" * 70)
sig = Counter()
for r in tum:
    s = []
    for lab in ORDER:
        b = r["bases"].get(lab)
        info = SNV[lab]
        if b is None:
            s.append("-")
        elif b == info["ref"]:
            s.append("R")          # ref
        else:
            s.append(b if b in ("A", "C", "G", "T") else ("D" if b == "DEL" else "?"))
    sig[" ".join(s)] += 1
print("pattern over (1)(2)(3)(4)(5)(6):  R=ref(matches SNV ref)  letter=base  D=DEL  -=nocov")
for pat, n in sig.most_common(30):
    print(f"  {n:>4}x   {pat}")

# ---------- Section 5: 'no complete REF read' check ----------
print("\n" + "=" * 70)
print("S5. INFERENCE-2 CHECK: any tumor read REF at ALL covered SNVs (cov>=4)?")
print("=" * 70)
allref = 0; checked = 0
for r in tum:
    if cov_count(r) >= 4:
        checked += 1
        if all((r["bases"][lab] == SNV[lab]["ref"]) for lab in ORDER if r["bases"].get(lab) is not None):
            allref += 1
print(f"reads cov>=4: {checked}; all-REF reads: {allref}")
