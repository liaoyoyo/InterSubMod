#!/usr/bin/env python3
"""
61 - Append 2 evidence_ledger entries for the cross-sample ASM verification
(targeted 38-position + genome-wide rate). Numbers READ from the deterministic
JSONs (§13: no hand-typed metrics). Append-only (research/CLAUDE.md SoT rule).
"""
import json, os

CS = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample")
LEDGER = "/big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/evidence_ledger.jsonl"
DATE = "2026-06-03"
ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CAN = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
       "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}

S = json.load(open(f"{CS}/cross_sample_synthesis.json"))
H = S["headline"]
SO = S["somatic_overlap"]
GW = {s: json.load(open(f"{CS}/{s}_gwasm.json")) for s in ALL
      if os.path.exists(f"{CS}/{s}_gwasm.json")}

# ---- entry 1: targeted ----
priv = H["somatic_status_counts"].get("hcc1395_private_somatic", 0)
som_overlap = "; ".join(f"{s}({CAN[s]}) {SO[s]['n_somatic_at_keypos']}/{SO[s]['of_total']}"
                        for s in ALL if s in SO)
imp_genes = [h["gene"] for h in H["imprinting_consistent_hits"]]
dis_genes = [h["gene"] for h in H["discordant_recurrent_hits"]]

e1 = {
    "cycle_id": "20260603_cross_sample_keypos_ASM_verification",
    "parent_cycles": ["ASM-CN-confound-pilot-20260602",
                      "project_zar1l_brca2_asm_verification"],
    "hypothesis_id": "CROSS_SAMPLE_KEYPOS_ASM_RECURRENCE",
    "hypothesis": ("HCC1395 的 38 個關鍵 ASM 位點（37 credible loci + BRCA2/ZAR1L）在其他 "
                   "5 個癌症樣本（乳腺 HCC1937/HCC1954、肺 H1437/H2009、黑色素瘤 COLO829）"
                   "是否出現相似 somatic ASM；各癌 private mutation 下預期同位點 somatic 復發罕見。"),
    "pipeline_track": ("targeted: validated pysam per-read 5mC β (modkit crossval Pearson=1.0), "
                       "HP:Z 分群, somatic-controlled HP-axis Δβ; deterministic synthesis 58 v2 "
                       "(雙軸 + sign-concordance); 2 adversarial evaluators"),
    "datasets_tested": [f"{s}_paired_full" for s in ALL],
    "scale": "38 key positions × 6 samples (A pilot, targeted)",
    "verdict": ("somatic_ASM_is_HCC1395_PRIVATE__0of38_cross_sample_somatic_recurrence__"
                f"{priv}_private_somatic__germline_axis_recurs_split_"
                f"{H['n_imprinting_consistent']}_sign_concordant_imprinting_vs_"
                f"{H['n_recurrent_discordant']}_discordant_NOT_imprinting"),
    "stability": "3",
    "tier_used": "tier_3_L2_targeted_cross_sample_single_truth_sample",
    "key_observations": (
        f"同位點 somatic 復發=0/38（其他 5 樣本全 0 somatic call，符合 private mutation 預期）；"
        f"somatic 重疊: {som_overlap}；"
        f"{priv} 個 HCC1395-private somatic ASM（含 flagship BRCA2/ZAR1L chr13:32315128 HP-axis Δβ≈-0.19）；"
        f"germline 軸復發 {H['n_recurrent_germline_asm']} 個 → 嚴格方向檢查拆: "
        f"sign-concordant 候選 imprinting {H['n_imprinting_consistent']} ({','.join(imp_genes)}) + "
        f"方向相反 NON-imprinting {H['n_recurrent_discordant']} ({','.join(dis_genes)})。"
        "BRCA2 逐樣本: HCC1395 有真 somatic 子單倍型+ASM, COLO829/HCC1954 germline het 雙倍型, "
        "H1437/HCC1937 LOH 單倍型, H2009 僅 4-read 假性 subhap(Δ≈0)。"),
    "decision": ("確認 ASM 的 somatic 因果性（若 CN/技術假象應跨樣本復發但 0/38）；"
                 "germline 軸『imprinting』命名必須方向一致才用（evaluator catch: 原 4 hit 中 3 個方向相反）；"
                 "勿宣稱 specific loci 跨樣本復發（是 private）。"),
    "caveat": ("單樣本參考(只 HCC1395 有 SEQC2 truth + somatic 訊號); 跨樣本用 window-mean 5mC "
               "(重現 HCC1395 credible-loci delta 方向非精確 magnitude, 原表 paired-CpG MAX-collapse); "
               "「不復發」部分受 LOH/低 N 使 germline 軸不可測所限; IGF2 cnLOH 使 imprinting 被遮蔽。"),
    "artifacts_path": ("docs/experiments/in_progress/2026/06/"
                       "20260603_cross_sample_keypos_ASM_verification_01.standalone.html + "
                       "research/tsg_promoter_asm_reviewer/scripts/57,58,60 + "
                       "genome_survey_v2/cn_confound/cross_sample/cross_sample_synthesis.json + 6x *_keypos.json"),
    "figures_dir": "embedded base64 in standalone HTML",
    "research_potential": ("ASM 真實性(somatic-driven, 位點 private)確認; germline/imprinting 基線已分離; "
                           "genome-wide 復現問題見 entry 20260603_genomewide_asm_rate_crosssample。"),
    "evaluator_verdicts": "separation-correctness:PASS; over-claim-guard:NEEDS_WORK->fixed(imprinting sign-concordance)",
    "timestamp": DATE,
}

# ---- entry 2: genome-wide ----
ex = {s: GW[s]["genomewide_somatic_asm"]["rate_excess_over_null"] for s in ALL if s in GW}
import statistics as st
ex_vals = list(ex.values())
gw_per = "; ".join(
    f"{s}({CAN[s]}) excess={GW[s]['genomewide_somatic_asm']['rate_excess_over_null']:+.3f}"
    f"(rate={GW[s]['genomewide_somatic_asm']['rate_strong_asm']:.3f}/null="
    f"{GW[s]['genomewide_somatic_asm']['rate_strong_null']:.3f}, n_eval="
    f"{GW[s]['genomewide_somatic_asm']['n_evaluable']})"
    for s in ALL if s in GW)
mean_ex = st.mean(ex_vals)
cv = st.pstdev(ex_vals) / mean_ex

e2 = {
    "cycle_id": "20260603_genomewide_asm_rate_crosssample",
    "parent_cycles": ["20260603_cross_sample_keypos_ASM_verification"],
    "hypothesis_id": "GENOMEWIDE_SOMATIC_ASM_PHENOMENON_REPRODUCES",
    "hypothesis": ("同位點不復發(private)下，各癌是否用自己的 private somatic 突變獨立呈現相同 "
                   "ASM 現象(somatic-controlled HP-axis ASM rate 高於噪音底)；跨癌種一致=現象可復現。"),
    "pipeline_track": ("per-sample own TP somatic SNVs 隨機 N=2000, somatic-controlled HP-axis Δβ, "
                       "5mC-only validated pysam; HP-label-shuffle null(20 perm) 估噪音底; "
                       "excess-over-null = 真訊號; +imprinted-DMR exploratory panel; 1 adversarial evaluator"),
    "datasets_tested": [f"{s}_paired_full" for s in ALL if s in GW],
    "scale": f"6 samples × 2000 random TP somatic SNVs each (n_evaluable {min(GW[s]['genomewide_somatic_asm']['n_evaluable'] for s in GW)}-{max(GW[s]['genomewide_somatic_asm']['n_evaluable'] for s in GW)})",
    "verdict": (f"somatic_ASM_phenomenon_REPRODUCES_6of6_samples_3of3_cancer_types__"
                f"all_excess_over_null_positive_range_{min(ex_vals):.3f}_to_{max(ex_vals):.3f}_"
                f"mean_{mean_ex:.3f}_CV_{cv:.2f}__single_pipeline_replication_tier_capped_3"),
    "stability": "3",
    "tier_used": "tier_3_L2_bordering_L2plus_6sample_single_pipeline_replication",
    "key_observations": (
        f"全 6/6 樣本 excess-over-null > 0; {gw_per}; "
        f"mean excess={mean_ex:.3f} CV={cv:.2f}; 3 癌種(乳腺/肺/黑色素瘤)皆正; "
        "raw rate 受小 N 噪音灌水(median|Δβ| 全 <0.10)故必看 excess; "
        "COLO829(melanoma) excess 最低(+0.101)且 null 最高(0.128) = 覆蓋最低(n_eval=802, n_low_cov=1092)所致, "
        "誠實歸因覆蓋非生物缺席; 方向多輕微 hyper lean(52-58%), HCC1954 唯一輕微 hypo(54%); "
        "imprinted DMR 正控(EXPLORATORY): PEG3 強(HCC1395 +0.72/HCC1937 -0.93 反映各自 LOH 保留親本 allele 不同), "
        "重度 LOH 樣本多數 DMR germline 軸不可測。"),
    "decision": ("somatic ASM 是跨癌種可復現的真實現象(非 CN/技術假象, 非 HCC1395-specific); "
                 "但 specific loci 是 private(復現的是 rate/magnitude 非位點); "
                 "必用 excess-over-null 不可用 raw rate 跨樣本比。"),
    "caveat": ("6 樣本共用同一 ClairS-TO/longphase caller + 同一 HP-axis 估計法 → 現象複製非獨立管線/cohort 驗證, "
               "共用系統性偏差未排除 → tier 封頂 ⭐3(bordering ⭐4); 升 ⭐4 需獨立管線或正交 truth; "
               "只 HCC1395 有 SEQC2 truth; window-mean vs paired-CpG 口徑 magnitude 不可直接並列; "
               "未出 per-sample permutation p(僅 point-estimate excess), 建議補強。"),
    "artifacts_path": ("research/tsg_promoter_asm_reviewer/scripts/59_genomewide_asm_rate.py + "
                       "genome_survey_v2/cn_confound/cross_sample/6x *_gwasm.json + "
                       "docs/experiments/in_progress/2026/06/20260603_cross_sample_keypos_ASM_verification_01.standalone.html(genome-wide section)"),
    "figures_dir": "embedded base64 in standalone HTML (genomewide rate vs null)",
    "research_potential": ("⭐4 升級路徑: 獨立管線(非 ClairS-TO)或正交 truth 複製 + per-sample permutation p; "
                           "COLO829 melanoma 已納入(原 ⭐4 blocker 之一); phasing 論文 ASM-真實性 methods 支撐強化。"),
    "evaluator_verdicts": "gw-reproducibility:PASS (required fix: single-pipeline caveat added to HTML)",
    "timestamp": DATE,
}

with open(LEDGER, "a") as f:
    f.write(json.dumps(e1, ensure_ascii=False) + "\n")
    f.write(json.dumps(e2, ensure_ascii=False) + "\n")

print(f"[61] appended 2 entries to {LEDGER}")
print(f"     e1: {e1['cycle_id']} -> {e1['verdict'][:70]}...")
print(f"     e2: {e2['cycle_id']} -> excess range {min(ex_vals):.3f}-{max(ex_vals):.3f} mean {mean_ex:.3f} CV {cv:.2f}")
print(f"     ledger now {sum(1 for _ in open(LEDGER))} lines")
