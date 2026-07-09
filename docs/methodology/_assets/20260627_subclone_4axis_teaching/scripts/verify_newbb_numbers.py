#!/usr/bin/env python3
"""
verify_newbb_numbers.py — 新骨幹(ClairS PASS)數字一鍵重算 + 取得/判別方式紀錄 + 新舊對照(2026-07-09)。
使用者要求:清楚紀錄數字取得與判別方式,可驗證與比對。每個數字都標「從哪個檔、用什麼判準」。
用法: python3 verify_newbb_numbers.py
"""
import os, json, gzip, glob
from collections import Counter

C = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
# 舊骨幹數字(is_somatic;先前記錄,供對照)
OLD = {"HCC1395": (23810, 7100, 56.2, 32.6), "HCC1395_DORADO": (10731, 3842, 55.2, 34.3),
       "COLO829": (25938, 7182, 78.7, 18.8), "H1437": (57986, 8407, 65.0, 19.3),
       "H2009": (128105, 9585, 61.9, 13.1), "HCC1937": (4870, 2036, 50.2, 50.4),
       "HCC1954": (9440, 3411, 91.3, 42.0)}

def somatic_pass_path(s):
    lp = glob.glob(f"{C}/{s}/paired_full/*complete_matrix/longphase_s")
    return f"{lp[0]}/somatic_pass.vcf.gz" if lp else None

def region_view_path(s):
    return f"{PILOT}/layered_region_view_HCC1395.json" if s == "HCC1395" else f"{MSROOT}/{s}/layered_region_view_{s}.json"

def count_somatic_pass(path):
    """骨幹 sSNV = somatic_pass.vcf.gz 的 SNV 位點數(全 PASS);同時算 NAF=0 佔比(=normal 無 ALT=真somatic 判準)。"""
    n = naf0 = allpass = 0
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.split("\t")
            if len(p) < 10 or len(p[3]) != 1 or len(p[4].strip()) != 1:
                continue
            n += 1
            if p[6] == "PASS":
                allpass += 1
            d = dict(zip(p[8].split(":"), p[9].split(":")))
            try:
                if float(d.get("NAF", "0")) == 0:
                    naf0 += 1
            except ValueError:
                pass
    return n, naf0, allpass

print("=" * 108)
print("新骨幹(ClairS PASS)數字一鍵重算 + 取得/判別方式 + 新舊對照")
print("=" * 108)
print("【數字取得/判別方式】")
print("  骨幹 sSNV : somatic_pass.vcf.gz 全 SNV 位點(caller=ClairS paired,FILTER=PASS;判準:LongPhase-S 輸出=ClairS 已確認 somatic)")
print("  NAF=0%    : somatic_pass FORMAT.NAF==0 佔比(NAF=normal allele freq;=0 → normal 無此 ALT → 真 somatic)")
print("  區/multiHP/all-det : layered_region_view_{S}.json census(region 主分母;all-det=全 germline lineage 都 determined)")
print("  V1-V7     : region-view L1.all_V1V7_pass(solver 獨立重算驗證)")
print("-" * 108)
print(f"{'樣本':15}{'骨幹sSNV 新/舊':>20}{'NAF=0%':>9}{'全PASS':>8}{'區 新/舊':>14}{'multiHP 新/舊':>15}{'all-det 新/舊':>15} V1-7")
print("-" * 108)
allok = True
for s in SAMPLES:
    sp = somatic_pass_path(s); rv = region_view_path(s)
    if not sp or not os.path.exists(sp) or not os.path.exists(rv):
        print(f"{s:15} 缺檔"); continue
    nsnv, naf0, allpass = count_somatic_pass(sp)
    Cc = json.load(open(rv))["census"]; nreg = Cc["n_regions"]; rd = Cc["region_determinacy"]; m = Cc["hp_multiplicity"]; L1 = Cc["L1"]
    u1 = Cc.get("U1_sSNV_somatic_total")
    mhp = 100 * m.get("2", 0) / nreg if nreg else 0
    det = 100 * rd.get("all_determined", 0) / nreg if nreg else 0
    v17 = L1.get("all_V1V7_pass"); allok = allok and v17
    o = OLD[s]
    # 一致性檢查:region-view U1 應 == somatic_pass 重算
    match = "✓" if u1 == nsnv else f"✗(rv={u1})"
    print(f"{s:15}{f'{nsnv}/{o[0]}':>20}{f'{100*naf0/nsnv:.1f}':>9}{f'{allpass}={nsnv}' if allpass==nsnv else f'{allpass}/{nsnv}':>8}"
          f"{f'{nreg}/{o[1]}':>14}{f'{mhp:.0f}/{o[2]:.0f}%':>15}{f'{det:.0f}/{o[3]:.0f}%':>15} {'✅' if v17 else '❌'} {match}")
print("-" * 108)
print(f"全 7 樣本 V1-V7: {'ALL PASS ✅' if allok else 'FAIL ❌'}")
print("驗證:骨幹sSNV 重算(somatic_pass)== region-view U1(欄末 ✓);全 PASS;NAF=0 高(真somatic)。可比對舊骨幹欄。")
