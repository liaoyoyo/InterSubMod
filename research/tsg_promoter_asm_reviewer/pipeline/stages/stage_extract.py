"""
stage_extract.py — P1 BAM-direct mod-type-aware extractor (Level-1-plus) + cache.

WHY (root cause from gap analysis): MSA Level-1 type-erases 5mC/5hmC (two unlabeled
rows) and is deleted after pivot. This emits a persistent Level-1-PLUS with an explicit
mod_code column ('m'=5mC / 'h'=5hmC), so any mod-type policy is a re-read not a re-run.

CRITICAL fix (PoC 36): window-restrict CpGs to ±window of the variant. ONT reads span
10-30kb and pysam read.modified_bases reports the WHOLE read's mod positions; without the
restriction, far-away CpGs dilute the focal signal (BRCA2 −0.122 -> −0.036 artifact).

Output columns (Level-1-plus):
  chrom somatic_pos methyl_pos read_id bam_source haplotype_tag somatic_allele_type mod_code meth_call strand

Guardrail: characterization extraction only; no TP/FP score. Cache writes go to a
gitignored cache dir (small per-locus; *.tsv.gz is gitignored except the one fixture).
"""
import os
import gzip
import hashlib
import pysam

WINDOW = 1000
COLS = ["chrom", "somatic_pos", "methyl_pos", "read_id", "bam_source",
        "haplotype_tag", "somatic_allele_type", "mod_code", "meth_call", "strand"]


def _cache_key(chrom, spos, tumor_bam, normal_bam, window):
    h = hashlib.sha1()
    for b in (tumor_bam, normal_bam):
        try:
            st = os.stat(b)
            h.update(f"{os.path.basename(b)}:{st.st_size}:{int(st.st_mtime)}".encode())
        except OSError:
            h.update(f"{b}:missing".encode())
    h.update(f"{chrom}:{spos}:{window}".encode())
    return h.hexdigest()[:16]


def extract_locus(chrom, spos, ref, alt, tumor_bam, normal_bam, window=WINDOW):
    """Walk tumor+normal reads in ±window; emit Level-1-plus rows (mod-type labeled,
    window-restricted). Returns list[dict]."""
    rows = []
    for src, bampath in (("tumor", tumor_bam), ("normal", normal_bam)):
        bam = pysam.AlignmentFile(bampath)
        for r in bam.fetch(chrom, spos - window, spos + window):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            try:
                hp = str(r.get_tag("HP"))
            except KeyError:
                hp = "untag"
            q2r = dict((q, rp) for q, rp in r.get_aligned_pairs(matches_only=True))
            r2q = {v: k for k, v in q2r.items()}
            al = "na"
            if (spos - 1) in r2q:
                b = r.query_sequence[r2q[spos - 1]].upper()
                al = "alt" if b == alt.upper() else ("ref" if b == ref.upper() else "oth")
            mb = r.modified_bases
            if not mb:
                continue
            for (canon, strand, mod), calls in mb.items():
                mc = "m" if str(mod) in ("m", "0") else ("h" if str(mod) in ("h", "17") else None)
                if mc is None:
                    continue
                for qpos, qual in calls:
                    rp = q2r.get(qpos)
                    if rp is None or abs(rp - (spos - 1)) > window:  # FIX: window-restrict
                        continue
                    rows.append({
                        "chrom": chrom, "somatic_pos": spos, "methyl_pos": rp,
                        "read_id": r.query_name, "bam_source": src, "haplotype_tag": hp,
                        "somatic_allele_type": al, "mod_code": mc,
                        "meth_call": round(qual / 255.0, 4), "strand": "+" if strand == 0 else "-",
                    })
        bam.close()
    return rows


def _default_cache_dir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cache", "level1")


def extract_or_cache(chrom, spos, ref, alt, tumor_bam, normal_bam,
                     cache_dir=None, window=WINDOW):
    """Persist Level-1-plus per locus keyed by BAM size/mtime+window. Returns (path, status).
    status: 'HIT' (re-read, no BAM walk) or 'MISS' (extracted+wrote)."""
    cache_dir = cache_dir or _default_cache_dir()
    os.makedirs(cache_dir, exist_ok=True)
    key = _cache_key(chrom, spos, tumor_bam, normal_bam, window)
    path = os.path.join(cache_dir, f"{chrom}_{spos}_{key}.tsv.gz")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path, "HIT"
    rows = extract_locus(chrom, spos, ref, alt, tumor_bam, normal_bam, window)
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt") as f:
        f.write("\t".join(COLS) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in COLS) + "\n")
    os.replace(tmp, path)
    return path, "MISS"
