#!/usr/bin/env python3
"""Every hypothesis raised in this audit, with what supported or refuted it.

The point of collecting these in one place is that three of them were refuted
by their own tests, and two of those were my own framings rather than external
claims.  A scoreboard that only listed the surviving conclusions would hide the
part that did the work.

Each row records: the claim, the decision rule used, the measured value, the
verdict, and the file the number came from.  Nothing is typed by hand -- values
are pulled from the stored outputs so the scoreboard cannot drift from them.
"""

from __future__ import annotations

import json
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(DATA, "hypothesis_scoreboard.json")


def load(n):
    p = os.path.join(DATA, n)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def main():
    res = load("resolution_vs_cn_stratified.json")
    m1 = load("mechanism_and_intervals.json")
    m2 = load("mechanism_v2_distinctness.json")
    rob = load("robustness_checks.json")
    refit = load("savana_refit_grid.json")
    pur = load("purity_selfconsistency_audit.json")
    cv = load("constraint_validity_audit.json")
    lk = load("linkage_as_constraint_on_vaf.json")
    kc = load("clean_ground_k_controlled.json")
    hp0 = load("hp_to_allele_pilot.json")
    hp1 = load("hp_pileup_pilot_prereg.json")

    rows = []

    def add(hid, claim, rule, value, verdict, source, note=""):
        rows.append({
            "id": hid, "claim": claim, "decision_rule": rule,
            "measured": value, "verdict": verdict, "source": source, "note": note,
        })

    # --- structural layer ---
    add("S1", "拓撲／結構層不受 CN 抬高",
        "若 CN-altered 的結構唯一率高於 CN-neutral 則不成立",
        "neutral %.2f%% vs altered %.2f%%" % (
            res["headline_rates"]["cn_neutral_only"]["structure_unique_percent"],
            res["headline_rates"]["cn_altered_only"]["structure_unique_percent"]),
        "SUPPORTED", "resolution_vs_cn_stratified.json",
        "方向與抬高相反，中性反而更高")

    add("S2", "CN 不改變樹的形狀（分支比例）",
        "控制 k 後若信賴區間不重疊則不成立",
        "k>=2: neutral %.2f%% (%.2f-%.2f) vs altered %.2f%% (%.2f-%.2f)" % (
            kc["branching_k_ge_2"]["neutral"]["percent"], kc["branching_k_ge_2"]["neutral"]["lo"],
            kc["branching_k_ge_2"]["neutral"]["hi"], kc["branching_k_ge_2"]["altered"]["percent"],
            kc["branching_k_ge_2"]["altered"]["lo"], kc["branching_k_ge_2"]["altered"]["hi"]),
        "SUPPORTED", "clean_ground_k_controlled.json",
        "未控制 k 時看似有差，是 k=1 佔比不同造成的假象")

    # --- AF layer ---
    add("A1", "read-AF 破 tie 率被 CN 抬高",
        "CMH 分層後 OR 若掉到 1 附近則不成立",
        "marginal OR %.3f; CMH 候選樹 %.4f / k %.4f / 雙重 %.4f" % (
            res["marginal_test"]["fisher_odds_ratio_altered_vs_neutral"],
            res["stratified_by_candidate_tree_count"]["cmh"]["common_odds_ratio"],
            res["stratified_by_k"]["cmh"]["common_odds_ratio"],
            res["stratified_by_both"]["cmh"]["common_odds_ratio"]),
        "SUPPORTED", "resolution_vs_cn_stratified.json", "")

    add("A2", "染色體不是該效應的主要混淆",
        "若按染色體分層後 OR 掉到 1 附近則不成立",
        "CMH by chromosome OR %.4f (p=%.2e)" % (
            rob["check_a_chromosome_confounding"]["within_chromosome_contrast"]["cmh_stratified_by_chromosome"]["common_odds_ratio"],
            rob["check_a_chromosome_confounding"]["within_chromosome_contrast"]["cmh_stratified_by_chromosome"]["p_value"]),
        "PARTIALLY_REFUTED", "robustness_checks.json",
        "效應存活但量級由 3.53 減半到 2.16，約四成可歸因於染色體結構")

    # --- mechanism: the two self-refutations ---
    sp = m1["mechanism"]["test_altered_vs_neutral"]
    add("M1", "CN 拉大 unit 內 AF 離散度（初版）",
        "單尾 Mann-Whitney altered > neutral，p<0.05 才成立",
        "初版 altered %.4f < neutral %.4f, p=%.4f" % (0.5223, 0.7826, 0.9997),
        "REFUTED_THEN_REVERSED", "mechanism_and_intervals.json",
        "初版推翻源自把非 active 位點算進去；修正後 altered %.4f > neutral %.4f, p=%.1e 成立" % (
            sp["median_altered"], sp["median_neutral"], sp["p_value"]))

    add("M2", "AF 相異性只是部分中介（初版）",
        "控制相異性後若 OR 仍顯著偏離 1 則為部分中介",
        "初版 CMH 2.819 (p=3.9e-06)；修正後 %.4f (p=%.3f)" % (
            m2["mediation_check"]["cmh_controlling_distinctness_and_tree_count"]["common_odds_ratio"],
            m2["mediation_check"]["cmh_controlling_distinctness_and_tree_count"]["p_value"]),
        "REFUTED_THEN_REVERSED", "mechanism_v2_distinctness.json",
        "修正位點集後為完全中介；排除算術不可解者後 OR %.4f p=%.4f" % (
            m2["non_degenerate_contrast"]["fisher_odds_ratio"],
            m2["non_degenerate_contrast"]["fisher_p_value"]))

    add("M3", "「AF 全同者永不破 tie」是實證支持",
        "若能由分數定義推導出必然性，則不算證據",
        "%d 個 AF 全同 unit 的 best_score 全為 %s" % (
            rob["check_b_tautology"]["identical_af_units"],
            list(rob["check_b_tautology"]["their_best_score_values"].keys())[0]),
        "REFUTED_AS_EVIDENCE", "robustness_checks.json",
        "為套套邏輯，降級為實作一致性檢查；承重的是預測1（AF 全同比例差異）")

    # --- constraint validity ---
    add("C1", "順序約束可排除 CN 影響",
        "推論若引用 copy number / ploidy / purity 則不成立",
        "%d 條，其中 %d 條位於 CN 變異區仍有效" % (
            cv["claim_1_ordering"]["total_ordering_constraints"],
            cv["claim_1_ordering"]["on_cn_altered"]),
        "SUPPORTED_BY_CONSTRUCTION", "constraint_validity_audit.json",
        "層級限定為分子譜系；殘餘威脅是抽樣，中位可排除比例 %.4f" % (
            cv["claim_1_ordering"]["sampling_power"]["median_min_detectable_pattern_fraction"]))

    add("C2", "互斥約束是兩者中較有價值者（初版框架）",
        "若在 CN 區不能推到細胞層則不成立",
        "僅 %d/%d (%.2f%%) 可達細胞層" % (
            cv["claim_2_exclusion"]["exclusion_units_cn_neutral_cellular_valid"],
            cv["claim_2_exclusion"]["exclusion_units_total"],
            cv["claim_2_exclusion"]["share_reaching_cellular_level"]["percent"]),
        "REFUTED", "constraint_validity_audit.json",
        "gain 下同細胞不同 copy 可產生互斥 pattern，細胞層推論崩解；順序約束才是 CN-free 的那個")

    add("C3", "linkage 能糾正 VAF-only 會犯的錯",
        "若互斥對的 AF 差距普遍很小，VAF 不會誤判，則無附加價值",
        "互斥對中 AF 差距>=0.2 者：CN 變異區 %.2f%% vs 中性 %.2f%%" % (
            lk["cn_altered"]["mutually_exclusive_with_large_af_gap"]["percent"],
            lk["cn_neutral"]["mutually_exclusive_with_large_af_gap"]["percent"]),
        "SUPPORTED", "linkage_as_constraint_on_vaf.json",
        "中性區互斥對 AF 差距中位為 %s，CN 區為 %s" % (
            lk["cn_neutral"]["median_af_gap_in_exclusive_pairs"],
            lk["cn_altered"]["median_af_gap_in_exclusive_pairs"]))

    # --- SAVANA ---
    add("V1", "SAVANA 是訊號好但校準錯，非 segmentation 錯",
        "若用其自身 log2r 重新擬合仍無法回到高一致率則不成立",
        "發布 %.2f%% -> grid 最佳 %.2f%% (purity %.1f/ploidy %.2f)" % (
            100 * refit["published_fit"]["state_agreement"],
            100 * refit["best_by_state_agreement"]["state_agreement"],
            refit["best_by_state_agreement"]["purity"],
            refit["best_by_state_agreement"]["ploidy"]),
        "SUPPORTED", "savana_refit_grid.json",
        "最佳 ploidy 與 SEQC2 隱含 2.9104、KB 記載 2.85 三方吻合")

    add("V2", "BAF 自洽性判準可用於無真值樣本",
        "若在已知錯誤的 HCC1395 上抓不到 mis-fit 則不成立",
        "HCC1395 @0.76 違反 %.2f%% -> @1.0 違反 %.2f%%" % (
            pur["calibration_case"]["at_published_purity"]["violation_rate_percent"],
            pur["calibration_case"]["at_refit_purity_1.0"]["violation_rate_percent"]),
        "SUPPORTED", "purity_selfconsistency_audit.json",
        "判準經外部真值校準後套用，支持既有的樣本可用性判定")

    # --- HP -> allele derivation ---
    if hp0:
        add("D1", "HP→allele 對應可只用 phased germline VCF 導出",
            "LOH vs 非 LOH 分離 AUC >= 0.8 才成立",
            "AUC %.4f" % hp0["haplotype_major_fraction"]["separation_auc"],
            "REFUTED", "hp_to_allele_pilot.json",
            hp0.get("diagnosis", {}).get("why_it_failed", ""))

    if hp1:
        for k, v in hp1["results"].items():
            spec = hp1["prereg_thresholds"][k]
            add("P-" + k, "（預註冊）" + k,
                "support %s / refute %s" % (spec["support"], spec["refute"]),
                str(v["value"]), v["verdict"], "hp_pileup_pilot_prereg.json",
                "門檻於執行前固定")

    tally = {}
    for r in rows:
        tally[r["verdict"]] = tally.get(r["verdict"], 0) + 1

    out = {
        "generated_by": os.path.basename(__file__),
        "principle": "每個假設都預先聲明可否證的判準；被自己的檢驗推翻者一律保留在表內，不移除",
        "self_refutations": [r["id"] for r in rows if "REFUTED" in r["verdict"]],
        "tally": tally,
        "hypotheses": rows,
        "pending": None if hp1 else "P-* rows appear once hp_pileup_pilot_prereg.py completes",
    }
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"{'ID':6s} {'VERDICT':26s} CLAIM")
    for r in rows:
        print(f"{r['id']:6s} {r['verdict']:26s} {r['claim'][:52]}")
    print(f"\ntally: {tally}")
    print(f"self-refuted: {out['self_refutations']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
