import { createHash } from "node:crypto";
import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { DatabaseSync } from "node:sqlite";

const repoRoot = resolve(import.meta.dirname, "../../..");
const outputPath = resolve(import.meta.dirname, "artifact.json");
const reportDatabasePath = resolve(import.meta.dirname, "data/report.sqlite");
const reportDatabaseSourcePath =
    "InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/data/report.sqlite";
const strictPath = resolve(
    repoRoot,
    "research/20260723_hcc1395_crosssource_topology_resolution/strict_pair_validation/results/strict_pairwise_metrics_01.json",
);
const vafArtifactPath = resolve(
    repoRoot,
    "research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/artifact.json",
);
const bridgePath = resolve(
    repoRoot,
    "research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/summary.json",
);
const exactRegionArtifactPath = resolve(
    repoRoot,
    "research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/artifact.json",
);

const expectedInputHashes = new Map([
    [strictPath, "aaadfa8751fbe6b78efb4b5070aeef9f6bbdfc14e67755bad52b6b9d443cfa41"],
    [vafArtifactPath, "1094b39be9d08a5314224a7ee9a56619a36fb22a6e4aacbf5def95b4088109ad"],
    [bridgePath, "8a8b489ea5b77bd19e6d1160374b77fdddac0eac3cdb760347603bb1f83ba3c8"],
    [exactRegionArtifactPath, "6497253851cd3fb3a6f043d90aee3685ae33ecd4f0cc084cd268325df04154a6"],
]);

const readVerifiedJson = (path) => {
    const bytes = readFileSync(path);
    const observedHash = createHash("sha256").update(bytes).digest("hex");
    const expectedHash = expectedInputHashes.get(path);
    if (observedHash !== expectedHash) {
        throw new Error(`Input hash mismatch for ${path}: ${observedHash}`);
    }
    return JSON.parse(bytes.toString("utf8"));
};

const strict = readVerifiedJson(strictPath);
const vafArtifact = readVerifiedJson(vafArtifactPath);
const bridge = readVerifiedJson(bridgePath);
const exactRegionArtifact = readVerifiedJson(exactRegionArtifactPath);

if (!strict.all_validations_pass) {
    throw new Error("Strict pairwise source did not pass all validations.");
}

const generatedAt = "2026-07-23T17:00:48+08:00";
const title = "HCC1395 biological-ID specificity 與多位點融合樹驗證";
const pct = (value) => value * 100;
const ratio = (numerator, denominator) => (denominator === 0 ? null : numerator / denominator);

const datasetQueries = {
    headline: 'SELECT * FROM "headline";',
    specificity: 'SELECT * FROM "specificity" ORDER BY "order";',
    specificity_chart:
        'SELECT * FROM "specificity_chart" ORDER BY "order", "pair_class";',
    shared_structure: 'SELECT * FROM "shared_structure" ORDER BY "order";',
    shared_structure_chart:
        'SELECT * FROM "shared_structure_chart" ORDER BY "order", "measure";',
    resolution: 'SELECT * FROM "resolution" ORDER BY "order";',
    evidence_denominators:
        'SELECT * FROM "evidence_denominators" ORDER BY "order";',
    denominator_chart:
        'SELECT * FROM "denominator_chart" ORDER BY "order", "measure";',
    evidence_ladder: 'SELECT * FROM "evidence_ladder" ORDER BY "order";',
    seqc2_comparison: 'SELECT * FROM "seqc2_comparison" ORDER BY "order";',
    fusion_gates: 'SELECT * FROM "fusion_gates" ORDER BY "gate_order";',
    tree_metrics: 'SELECT * FROM "tree_metrics" ORDER BY "order";',
    claim_tiers: 'SELECT * FROM "claim_tiers" ORDER BY "level";',
};

const metricLabels = {
    candidate_sites: "Candidate sSNV",
    active_sites: "Active sSNV",
    w_sites: "W 內 sSNV",
    primary_edges: "無向 direct edge",
    alt_informative_edges: "ALT-informative edge",
    aa_edges: "AA edge",
    exact_components: "Exact component",
};

const specificityRows = Object.entries(strict.target_ranks).map(([key, row], index) => ({
    order: index + 1,
    layer: metricLabels[key],
    target_jaccard: row.target_jaccard,
    target_jaccard_pct: pct(row.target_jaccard),
    target_rank: row.target_rank_of_21,
    next_unrelated_pair: row.next_unrelated_pair.join(" × "),
    next_unrelated_jaccard: row.next_unrelated_jaccard,
    next_unrelated_jaccard_pct: pct(row.next_unrelated_jaccard),
    target_to_next_ratio: row.target_to_next_ratio,
}));

const specificityChartRows = specificityRows.flatMap((row) => [
    {
        order: row.order,
        layer: row.layer,
        pair_class: "HCC1395 same biological ID",
        jaccard: row.target_jaccard,
        rank: row.target_rank,
        comparator_pair: "HCC1395 × HCC1395_DORADO",
    },
    {
        order: row.order,
        layer: row.layer,
        pair_class: "次佳不同 biological ID",
        jaccard: row.next_unrelated_jaccard,
        rank: null,
        comparator_pair: row.next_unrelated_pair,
    },
]);

const target = strict.target_pair;
const sharedProjection = target.shared_candidate_projection;
const sharedRows = [
    {
        order: 1,
        layer: "Candidate sSNV",
        unit: "exact allele",
        hcc_n: target.candidate_sites.left_n,
        dorado_n: target.candidate_sites.right_n,
        intersection_n: target.candidate_sites.intersection_n,
        jaccard: target.candidate_sites.jaccard,
        dorado_in_hcc: target.candidate_sites.right_recall,
        universe: "完整 strict site catalog",
    },
    {
        order: 2,
        layer: "Active sSNV",
        unit: "exact allele",
        hcc_n: sharedProjection.active_sites.left_n,
        dorado_n: sharedProjection.active_sites.right_n,
        intersection_n: sharedProjection.active_sites.intersection_n,
        jaccard: sharedProjection.active_sites.jaccard,
        dorado_in_hcc: sharedProjection.active_sites.right_recall,
        universe: "76,721 shared candidate alleles",
    },
    {
        order: 3,
        layer: "W 內 sSNV",
        unit: "exact allele",
        hcc_n: sharedProjection.w_sites.left_n,
        dorado_n: sharedProjection.w_sites.right_n,
        intersection_n: sharedProjection.w_sites.intersection_n,
        jaccard: sharedProjection.w_sites.jaccard,
        dorado_in_hcc: sharedProjection.w_sites.right_recall,
        universe: "76,721 shared candidate alleles",
    },
    {
        order: 4,
        layer: "無向 direct edge",
        unit: "allele pair",
        hcc_n: sharedProjection.primary_edges.left_n,
        dorado_n: sharedProjection.primary_edges.right_n,
        intersection_n: sharedProjection.primary_edges.intersection_n,
        jaccard: sharedProjection.primary_edges.jaccard,
        dorado_in_hcc: sharedProjection.primary_edges.right_recall,
        universe: "兩端皆在 shared candidates",
    },
    {
        order: 5,
        layer: "ALT-informative edge",
        unit: "allele pair",
        hcc_n: sharedProjection.alt_informative_edges.left_n,
        dorado_n: sharedProjection.alt_informative_edges.right_n,
        intersection_n: sharedProjection.alt_informative_edges.intersection_n,
        jaccard: sharedProjection.alt_informative_edges.jaccard,
        dorado_in_hcc: sharedProjection.alt_informative_edges.right_recall,
        universe: "兩端皆在 shared candidates",
    },
    {
        order: 6,
        layer: "AA edge",
        unit: "allele pair",
        hcc_n: sharedProjection.aa_edges.left_n,
        dorado_n: sharedProjection.aa_edges.right_n,
        intersection_n: sharedProjection.aa_edges.intersection_n,
        jaccard: sharedProjection.aa_edges.jaccard,
        dorado_in_hcc: sharedProjection.aa_edges.right_recall,
        universe: "兩端皆在 shared candidates",
    },
    {
        order: 7,
        layer: "Component co-membership",
        unit: "allele pair",
        hcc_n: sharedProjection.co_membership_pairs.left_n,
        dorado_n: sharedProjection.co_membership_pairs.right_n,
        intersection_n: sharedProjection.co_membership_pairs.intersection_n,
        jaccard: sharedProjection.co_membership_pairs.jaccard,
        dorado_in_hcc: sharedProjection.co_membership_pairs.right_recall,
        universe: "兩端皆在 shared candidates",
    },
].map((row) => ({
    ...row,
    jaccard_pct: pct(row.jaccard),
    dorado_in_hcc_pct: pct(row.dorado_in_hcc),
}));

const sharedChartRows = sharedRows.flatMap((row) => [
    {
        order: row.order,
        layer: row.layer,
        measure: "對稱 Jaccard",
        value: row.jaccard,
        numerator: row.intersection_n,
        denominator_note: row.universe,
    },
    {
        order: row.order,
        layer: row.layer,
        measure: "DORADO → HCC containment",
        value: row.dorado_in_hcc,
        numerator: row.intersection_n,
        denominator_note: `${row.dorado_n} DORADO ${row.unit}`,
    },
]);

const datasetSummary = exactRegionArtifact.snapshot.datasets.dataset_summary;
const hccRegion = datasetSummary.find((row) => row.dataset === "HCC1395");
const doradoRegion = datasetSummary.find((row) => row.dataset === "HCC1395_DORADO");
if (!hccRegion || !doradoRegion) {
    throw new Error("Exact-region dataset summary is missing the HCC pair.");
}

const resolutionRows = [
    {
        order: 1,
        metric: "Exact chromosome×PS",
        unit: "count",
        hcc: hccRegion.exact_chromosome_PS,
        dorado: doradoRegion.exact_chromosome_PS,
        interpretation: "DORADO phase blocks 較碎",
    },
    {
        order: 2,
        metric: "Exact PS×HP containers",
        unit: "count",
        hcc: hccRegion.exact_PS_HP_containers,
        dorado: doradoRegion.exact_PS_HP_containers,
        interpretation: "DORADO hard boundaries 較多",
    },
    {
        order: 3,
        metric: "Median W span",
        unit: "bp",
        hcc: hccRegion.median_W_span_bp,
        dorado: doradoRegion.median_W_span_bp,
        interpretation: "DORADO linkage reach 較短",
    },
    {
        order: 4,
        metric: "Median direct-edge distance",
        unit: "bp",
        hcc: hccRegion.median_direct_edge_distance_bp,
        dorado: doradoRegion.median_direct_edge_distance_bp,
        interpretation: "DORADO 長距離 edge 較少",
    },
    {
        order: 5,
        metric: "Canonical molecule rows",
        unit: "count",
        hcc: hccRegion.canonical_molecule_rows,
        dorado: doradoRegion.canonical_molecule_rows,
        interpretation: "差異不是單純總 molecule 數不足",
    },
].map((row) => ({
    ...row,
    dorado_over_hcc: ratio(row.dorado, row.hcc),
    dorado_over_hcc_pct: pct(ratio(row.dorado, row.hcc)),
}));

const direction = bridge.join_coverage;
const multiCluster = bridge.k_gt_1_informative_both_multicluster;
const evidenceDenominatorRows = [
    {
        order: 1,
        evidence: "共同 candidate 中同時進 W",
        parent_universe: "shared candidate alleles",
        parent_n: sharedProjection.shared_candidate_sites_n,
        informative_n: target.w_sites.intersection_n,
        informative_coverage: target.joint_w_coverage_of_shared_candidates,
        conditional_agreement: sharedProjection.w_sites.right_recall,
        agreement_definition: "DORADO W allele 被 HCC W 包含",
        limitation: "只涵蓋共同 candidate 的 13.23%",
    },
    {
        order: 2,
        evidence: "兩來源皆可判定方向",
        parent_universe: "read-directed relation pairs",
        parent_n: direction.both_sources_determinate_direction.denominator,
        informative_n: direction.both_sources_determinate_direction.n,
        informative_coverage: direction.both_sources_determinate_direction.share,
        conditional_agreement: direction.same_direction_among_both_determinate.share,
        agreement_definition: "596/598 jointly determinate pairs 同向",
        limitation: "可判定方向 coverage 只有 7.39%",
    },
    {
        order: 3,
        evidence: "兩側皆 multi-cluster region",
        parent_universe: "fixed multi-locus regions",
        parent_n: bridge.k_gt_1_region_pattern.fixed_regions_in_stratum,
        informative_n: multiCluster.denominator,
        informative_coverage: ratio(
            multiCluster.denominator,
            bridge.k_gt_1_region_pattern.fixed_regions_in_stratum,
        ),
        conditional_agreement: multiCluster.exact_share,
        agreement_definition: "21/34 informative regions partition exact",
        limitation: "n=34；不可用 5,028/5,189 取代",
    },
].map((row) => ({
    ...row,
    informative_coverage_pct: pct(row.informative_coverage),
    conditional_agreement_pct: pct(row.conditional_agreement),
}));

const denominatorChartRows = evidenceDenominatorRows.flatMap((row) => [
    {
        order: row.order,
        evidence: row.evidence,
        measure: "條件式一致",
        value: row.conditional_agreement,
        numerator: row.informative_n,
        denominator: row.parent_n,
    },
    {
        order: row.order,
        evidence: row.evidence,
        measure: "可辨識覆蓋",
        value: row.informative_coverage,
        numerator: row.informative_n,
        denominator: row.parent_n,
    },
]);

const hccMetrics = Object.fromEntries(
    vafArtifact.snapshot.datasets.hcc_metrics.map((row) => [row.metric, row]),
);
const vafHeadline = vafArtifact.snapshot.datasets.headline_metrics[0];
const pycloneMain = bridge.global_concordance;
const evidenceLadderRows = [
    {
        order: 1,
        layer: "Exact allele identity",
        primary_metric: "Jaccard",
        estimate: strict.target_pair.candidate_sites.jaccard,
        estimate_text: "0.9276（76,721 shared alleles）",
        interpretation: "同一 biological ID 底層位點高度重現",
        claim_level: "L1 runtime",
    },
    {
        order: 2,
        layer: "Raw VAF backbone",
        primary_metric: "CCC",
        estimate: hccMetrics.vaf_ccc.value,
        estimate_text: `0.9339（95% block-bootstrap ${hccMetrics.vaf_ccc.ci95_low_1mb_block_bootstrap.toFixed(4)}–${hccMetrics.vaf_ccc.ci95_high_1mb_block_bootstrap.toFixed(4)}）`,
        interpretation: "主要頻率結構高度相近",
        claim_level: "L2 technical reproducibility",
    },
    {
        order: 3,
        layer: "Independent PyClone clustering",
        primary_metric: "ARI",
        estimate: pycloneMain.ari,
        estimate_text: "0.5389",
        interpretation: "全體 cluster agreement 只有中度",
        claim_level: "L3 conditional model",
    },
    {
        order: 4,
        layer: "Minor/subclonal membership",
        primary_metric: "Subclonal Jaccard",
        estimate: pycloneMain.subclonal_jaccard,
        estimate_text: "0.3810（277/727）",
        interpretation: "精細 minor-clone membership 尚不高度重現",
        claim_level: "L3 conditional model",
    },
    {
        order: 5,
        layer: "Informative multi-cluster regions",
        primary_metric: "Partition exact",
        estimate: multiCluster.exact_share,
        estimate_text: "0.6176（21/34）",
        interpretation: "多位點分群相容性中度且樣本數小",
        claim_level: "L3 regional bridge",
    },
    {
        order: 6,
        layer: "Global directed fusion tree",
        primary_metric: "Strict candidate-tree set",
        estimate: null,
        estimate_text: "尚未產出 DORADO strict topology",
        interpretation: "完整 ancestry 仍不可驗證",
        claim_level: "L5 absent",
    },
];

const seqc2ComparisonRows = [
    {
        order: 1,
        item: "Node universe",
        seqc2_reference: "Bulk WGS/WES somatic SNV + CNA",
        intersubmod_current: "76,721 exact shared candidate sSNVs；strict read-linkage",
        valid_comparison: "投影到共同 exact alleles，另保留 CN context",
        verdict: "可對齊",
    },
    {
        order: 2,
        item: "Root / trunk",
        seqc2_reference: "S1 為 MRCA，多數 driver events 在 trunk",
        intersubmod_current: "候選位點與 raw VAF backbone 高度重現",
        valid_comparison: "比較 clonal mutation set、CCF 與 ancestor relations",
        verdict: "大方向合理；未逐邊驗證",
    },
    {
        order: 3,
        item: "Branching",
        seqc2_reference: "兩個主要分支從 S2、S8 開始",
        intersubmod_current: "局部無向骨架；DORADO 多為 HCC contraction",
        valid_comparison: "rooted-triplet / collapsed-tree compatibility",
        verdict: "尚待 strict tree",
    },
    {
        order: 4,
        item: "Clone fractions",
        seqc2_reference: "PhyloWGS CCF；例 S2=60%",
        intersubmod_current: "raw VAF + conditional PyClone-VI",
        valid_comparison: "allele-specific CN、purity、multiplicity-aware CCF",
        verdict: "部分可做；source-specific CN 不足",
    },
    {
        order: 5,
        item: "Orthogonal support",
        seqc2_reference: "SuperFreq、subHMM、1,270 tumor single cells",
        intersubmod_current: "跨技術 ONT 重現 + 20 unrelated negative pairs",
        valid_comparison: "再加入 single-cell / multi-region / independent CN",
        verdict: "technical specificity 強；biological truth 未滿",
    },
    {
        order: 6,
        item: "Truth status",
        seqc2_reference: "多模型支持的 inferred reference architecture",
        intersubmod_current: "候選區域與局部 constraint evidence",
        valid_comparison: "不可用肉眼相似或 S1–S10 逐邊硬比",
        verdict: "Figure 4 不是唯一 tree truth",
    },
];

const fusionGateRows = [
    {
        gate_order: 1,
        gate: "G0 固定 biological identity / node universe",
        required_evidence: "預註冊 same-ID pair；exact-allele common universe",
        current_value: "J=92.76%；rank 1/21",
        status: "PASS",
        claim_if_pass: "結果與 gross sample mix-up 不相容",
    },
    {
        gate_order: 2,
        gate: "G1 局部無向骨架",
        required_evidence: "shared-opportunity edge / co-membership containment",
        current_value: "edge 98.28%；co-membership 90.05%",
        status: "PASS（條件式）",
        claim_if_pass: "DORADO 可辨識的局部無向連結大多被 HCC 包含",
    },
    {
        gate_order: 3,
        gate: "G2 可辨識廣度",
        required_evidence: "joint W 與 determinate direction coverage",
        current_value: "13.23%；7.39%",
        status: "NOT PASS",
        claim_if_pass: "能把局部相容推廣到大部分共同位點",
    },
    {
        gate_order: 4,
        gate: "G3 CN/CCF-aware minor-clone stability",
        required_evidence: "獨立 fits；minor membership 與 fraction uncertainty",
        current_value: "minor J=0.381；ARI=0.539",
        status: "PARTIAL",
        claim_if_pass: "精細 subclone assignment 可跨來源重現",
    },
    {
        gate_order: 5,
        gate: "G4 Strict candidate-tree set",
        required_evidence: "兩來源獨立 strict-bound solver；同 node universe",
        current_value: "DORADO topology 尚未完成",
        status: "NOT RUN",
        claim_if_pass: "可比較有向 topology 而非無向 graph",
    },
    {
        gate_order: 6,
        gate: "G5 Consensus / fusion tree",
        required_evidence: "ancestor matrix、triplets、collapsed RF、bootstrap",
        current_value: "尚無 frozen candidate-tree set",
        status: "NOT RUN",
        claim_if_pass: "可輸出 collapsed consensus tree",
    },
    {
        gate_order: 7,
        gate: "G6 Orthogonal lineage truth",
        required_evidence: "single-cell / multi-region / independent CN-CCF",
        current_value: "SEQC2 Figure 4 僅作外部 reference",
        status: "CONTEXT ONLY",
        claim_if_pass: "升格為 biological clone genealogy",
    },
];

const treeMetricRows = [
    {
        order: 1,
        object: "共同祖先關係",
        primary_metric: "Ancestor-pair concordance",
        must_report: "jointly resolved coverage + opposite rate",
        fusion_rule: "雙側同向才進 strict consensus",
    },
    {
        order: 2,
        object: "三葉分枝",
        primary_metric: "Rooted-triplet concordance",
        must_report: "resolved-triplet coverage",
        fusion_rule: "未解 triplet 保留 polytomy",
    },
    {
        order: 3,
        object: "Clade / contraction",
        primary_metric: "Generalized RF + maximum agreement subtree",
        must_report: "common-leaf / retained-node fraction",
        fusion_rule: "先 collapse 低支持 edge 再比",
    },
    {
        order: 4,
        object: "候選樹集合",
        primary_metric: "Minimum + Hausdorff set distance",
        must_report: "兩側 candidate counts",
        fusion_rule: "至少一組相容才稱 compatible",
    },
    {
        order: 5,
        object: "Clone assignment",
        primary_metric: "Minor-only Jaccard、ARI、NMI",
        must_report: "clonal-majority fraction",
        fusion_rule: "不得用 overall accuracy 掩蓋 minor instability",
    },
    {
        order: 6,
        object: "Clone fractions",
        primary_metric: "CCF CCC、MAE、interval overlap",
        must_report: "CN / purity / multiplicity sensitivity",
        fusion_rule: "source-specific CN 優先；共享 CN 只作 sensitivity",
    },
    {
        order: 7,
        object: "穩定性與 specificity",
        primary_metric: "Chromosome×PS block bootstrap + 21-pair rank",
        must_report: "coverage-matched negative controls",
        fusion_rule: "HCC pair 的 CI 應與 unrelated maximum 分離",
    },
];

const claimTierRows = [
    {
        level: "L1",
        evidence: "執行產物與守恆檢查",
        current_result: "7 datasets × chr1–22；21 pairs；all validations PASS",
        allowed_claim: "數值與比較母群可重現",
    },
    {
        level: "L2",
        evidence: "跨資料集 descriptive specificity",
        current_result: "HCC pair 七層皆 rank 1/21",
        allowed_claim: "方法保留強 biological-ID signal",
    },
    {
        level: "L3",
        evidence: "獨立 PyClone / region bridge",
        current_result: "minor J=0.381；21/34 multi-cluster exact",
        allowed_claim: "bulk raw-VAF backbone 穩、minor 結構中度",
    },
    {
        level: "L4",
        evidence: "SEQC2 bulk + single-cell external reference",
        current_result: "HCC1395 確有 branching heterogeneity",
        allowed_claim: "目標生物學形態合理",
    },
    {
        level: "L5",
        evidence: "Exact directed clone-tree truth",
        current_result: "不存在；DORADO strict tree 未完成",
        allowed_claim: "不可宣稱完整融合樹已證實",
    },
];

const headline = [
    {
        exact_candidate_jaccard: strict.target_pair.candidate_sites.jaccard,
        target_rank: strict.target_ranks.candidate_sites.target_rank_of_21,
        shared_candidate_alleles: strict.target_pair.candidate_sites.intersection_n,
        shared_edge_containment: sharedProjection.primary_edges.right_recall,
        shared_comembership_containment: sharedProjection.co_membership_pairs.right_recall,
        raw_vaf_ccc: hccMetrics.vaf_ccc.value,
        joint_w_coverage: target.joint_w_coverage_of_shared_candidates,
        determinate_direction_coverage: direction.both_sources_determinate_direction.share,
        minor_subclone_jaccard: vafHeadline.independent_subclone_jaccard,
    },
];

const seqc2Schematic = `
<section aria-label="SEQC2 reference architecture and current InterSubMod evidence" style="padding:8px 0 2px">
  <svg viewBox="0 0 1240 430" role="img" aria-labelledby="seqc2-title seqc2-desc"
       style="width:100%;height:auto;display:block;font-family:Inter,'Noto Sans TC',sans-serif">
    <title id="seqc2-title">SEQC2 reference architecture 與目前 InterSubMod 證據層級</title>
    <desc id="seqc2-desc">左側為抽象化的 SEQC2 S1 MRCA 與 S2、S8 兩個主要分支；右側顯示 shared-opportunity 上 HCC1395 的局部無向連結與 DORADO 較稀疏 subset。此圖是比較框架，不是新推論的 clone tree。</desc>
    <rect x="18" y="18" width="564" height="346" rx="20" fill="#f7f9fb" stroke="#cbd5df" stroke-width="2"/>
    <text x="42" y="55" font-size="22" font-weight="700" fill="#17324a">SEQC2 Figure 4：外部參考架構</text>
    <text x="42" y="82" font-size="14" fill="#5e6f80">PhyloWGS inference；SuperFreq / subHMM / single-cell CNV 支持異質性</text>
    <line x1="296" y1="310" x2="296" y2="245" stroke="#35556f" stroke-width="6"/>
    <circle cx="296" cy="320" r="24" fill="#35556f"/>
    <text x="296" y="326" text-anchor="middle" font-size="17" font-weight="700" fill="#fff">S1</text>
    <text x="338" y="324" font-size="15" fill="#314658">MRCA / trunk</text>
    <line x1="296" y1="245" x2="176" y2="163" stroke="#35556f" stroke-width="5"/>
    <line x1="296" y1="245" x2="416" y2="163" stroke="#a66a09" stroke-width="5"/>
    <circle cx="176" cy="151" r="22" fill="#35556f"/>
    <circle cx="416" cy="151" r="22" fill="#a66a09"/>
    <text x="176" y="157" text-anchor="middle" font-size="16" font-weight="700" fill="#fff">S2</text>
    <text x="416" y="157" text-anchor="middle" font-size="16" font-weight="700" fill="#fff">S8</text>
    <line x1="176" y1="129" x2="112" y2="94" stroke="#6f8799" stroke-width="3"/>
    <line x1="176" y1="129" x2="234" y2="91" stroke="#6f8799" stroke-width="3"/>
    <line x1="416" y1="129" x2="360" y2="91" stroke="#c58b2c" stroke-width="3"/>
    <line x1="416" y1="129" x2="482" y2="94" stroke="#c58b2c" stroke-width="3"/>
    <circle cx="112" cy="88" r="10" fill="#93a7b7"/>
    <circle cx="234" cy="85" r="10" fill="#93a7b7"/>
    <circle cx="360" cy="85" r="10" fill="#d7ad65"/>
    <circle cx="482" cy="88" r="10" fill="#d7ad65"/>
    <text x="42" y="393" font-size="14" fill="#5e6f80">抽象化示意：只保留論文明示的 S1 trunk 與 S2 / S8 主要分支。</text>

    <rect x="620" y="18" width="602" height="346" rx="20" fill="#f7f9fb" stroke="#cbd5df" stroke-width="2"/>
    <text x="644" y="55" font-size="22" font-weight="700" fill="#17324a">InterSubMod：目前能支持的結構</text>
    <text x="644" y="82" font-size="14" fill="#5e6f80">同一 shared allele universe；solid=共同可觀測，dashed=未解方向</text>
    <line x1="760" y1="308" x2="760" y2="226" stroke="#28715a" stroke-width="8"/>
    <line x1="760" y1="226" x2="708" y2="170" stroke="#28715a" stroke-width="6"/>
    <line x1="760" y1="226" x2="820" y2="170" stroke="#28715a" stroke-width="6"/>
    <line x1="708" y1="170" x2="680" y2="126" stroke="#28715a" stroke-width="4"/>
    <line x1="820" y1="170" x2="858" y2="124" stroke="#28715a" stroke-width="4"/>
    <text x="668" y="338" font-size="16" font-weight="700" fill="#164f3d">HCC1395</text>
    <text x="655" y="359" font-size="13" fill="#47685e">較多可觀測局部無向連結</text>

    <line x1="1004" y1="308" x2="1004" y2="226" stroke="#35556f" stroke-width="8"/>
    <line x1="1004" y1="226" x2="952" y2="170" stroke="#35556f" stroke-width="6"/>
    <line x1="1004" y1="226" x2="1064" y2="170" stroke="#35556f" stroke-width="6" stroke-dasharray="12 8"/>
    <line x1="952" y1="170" x2="924" y2="126" stroke="#35556f" stroke-width="4"/>
    <text x="932" y="338" font-size="16" font-weight="700" fill="#17324a">DORADO</text>
    <text x="904" y="359" font-size="13" fill="#526779">與較稀疏 subset / contraction 相容</text>

    <line x1="866" y1="246" x2="910" y2="246" stroke="#46556b" stroke-width="2"/>
    <path d="M902,239 L914,246 L902,253" fill="none" stroke="#46556b" stroke-width="2"/>
    <text x="871" y="226" font-size="13" fill="#526779">共同骨架比較</text>
    <rect x="620" y="382" width="602" height="31" rx="10" fill="#fff6df" stroke="#d7ad65"/>
    <text x="921" y="403" text-anchor="middle" font-size="14" fill="#5b430f">這是 evidence map，不是已完成的 HCC1395 × DORADO global fusion tree。</text>
  </svg>
</section>`;

const fusionFlow = `
<section aria-label="Multi-locus fusion tree validation flow" style="padding:8px 0 2px">
  <svg viewBox="0 0 1240 360" role="img" aria-labelledby="fusion-title fusion-desc"
       style="width:100%;height:auto;display:block;font-family:Inter,'Noto Sans TC',sans-serif">
    <title id="fusion-title">多位點融合樹的合理建構與驗證流程</title>
    <desc id="fusion-desc">共同位點母群經兩來源獨立局部候選樹、CN與CCF校正、候選樹集合比較、收縮未支持邊，最後才形成帶有未解多分枝的 consensus tree。</desc>
    <defs>
      <marker id="fusion-arrow" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
        <path d="M0,0 L0,6 L9,3 z" fill="#46556b"/>
      </marker>
    </defs>
    <rect x="20" y="55" width="180" height="112" rx="16" fill="#eef4f8" stroke="#35556f" stroke-width="2"/>
    <text x="110" y="87" text-anchor="middle" font-size="18" font-weight="700" fill="#17324a">共同位點母群</text>
    <text x="110" y="117" text-anchor="middle" font-size="14" fill="#405467">exact allele</text>
    <text x="110" y="142" text-anchor="middle" font-size="14" fill="#405467">opportunity frozen</text>
    <line x1="200" y1="111" x2="240" y2="111" stroke="#46556b" stroke-width="3" marker-end="url(#fusion-arrow)"/>

    <rect x="250" y="28" width="210" height="80" rx="16" fill="#e8f5ef" stroke="#28715a" stroke-width="2"/>
    <text x="355" y="61" text-anchor="middle" font-size="17" font-weight="700" fill="#164f3d">HCC candidate-tree set</text>
    <text x="355" y="87" text-anchor="middle" font-size="13" fill="#375f53">independent strict solver</text>
    <rect x="250" y="132" width="210" height="80" rx="16" fill="#eef4f8" stroke="#35556f" stroke-width="2"/>
    <text x="355" y="165" text-anchor="middle" font-size="17" font-weight="700" fill="#17324a">DOR candidate-tree set</text>
    <text x="355" y="191" text-anchor="middle" font-size="13" fill="#405467">independent strict solver</text>

    <line x1="460" y1="68" x2="515" y2="106" stroke="#46556b" stroke-width="3" marker-end="url(#fusion-arrow)"/>
    <line x1="460" y1="172" x2="515" y2="134" stroke="#46556b" stroke-width="3" marker-end="url(#fusion-arrow)"/>
    <rect x="525" y="55" width="220" height="112" rx="16" fill="#fff6df" stroke="#a66a09" stroke-width="2"/>
    <text x="635" y="87" text-anchor="middle" font-size="18" font-weight="700" fill="#573b0b">CN / purity / CCF</text>
    <text x="635" y="117" text-anchor="middle" font-size="14" fill="#5e513b">ancestor relation matrix</text>
    <text x="635" y="142" text-anchor="middle" font-size="14" fill="#5e513b">uncertainty preserved</text>
    <line x1="745" y1="111" x2="790" y2="111" stroke="#46556b" stroke-width="3" marker-end="url(#fusion-arrow)"/>

    <rect x="800" y="55" width="190" height="112" rx="16" fill="#eef4f8" stroke="#35556f" stroke-width="2"/>
    <text x="895" y="87" text-anchor="middle" font-size="18" font-weight="700" fill="#17324a">Tree-set comparison</text>
    <text x="895" y="117" text-anchor="middle" font-size="14" fill="#405467">triplet / RF / set distance</text>
    <text x="895" y="142" text-anchor="middle" font-size="14" fill="#405467">coverage + conflict</text>
    <line x1="990" y1="111" x2="1035" y2="111" stroke="#46556b" stroke-width="3" marker-end="url(#fusion-arrow)"/>

    <rect x="1045" y="55" width="175" height="112" rx="16" fill="#e8f5ef" stroke="#28715a" stroke-width="2"/>
    <text x="1132" y="87" text-anchor="middle" font-size="18" font-weight="700" fill="#164f3d">Consensus tree</text>
    <text x="1132" y="117" text-anchor="middle" font-size="14" fill="#375f53">supported trunk</text>
    <text x="1132" y="142" text-anchor="middle" font-size="14" fill="#375f53">polytomy / dashed edges</text>

    <rect x="20" y="255" width="1200" height="74" rx="16" fill="#f7f8fa" stroke="#d9dee3"/>
    <circle cx="52" cy="282" r="7" fill="#28715a"/>
    <text x="70" y="288" font-size="14" fill="#3e4d5a">雙側同向 → strict consensus</text>
    <circle cx="326" cy="282" r="7" fill="#a66a09"/>
    <text x="344" y="288" font-size="14" fill="#3e4d5a">單側可觀測 → provisional / dashed</text>
    <circle cx="684" cy="282" r="7" fill="#9c3f48"/>
    <text x="702" y="288" font-size="14" fill="#3e4d5a">雙側反向 → conflict，不融合</text>
    <circle cx="982" cy="282" r="7" fill="#78838d"/>
    <text x="1000" y="288" font-size="14" fill="#3e4d5a">雙側未解 → polytomy</text>
    <text x="52" y="316" font-size="13" fill="#657687">跨區沒有物理 linkage 時，保留 candidate ensemble；不得直接把所有 local edges 聯集成單一樹。</text>
  </svg>
</section>`;

const sqlSource = (id, label, dataset, description, metricDefinitions = []) => ({
    id,
    label,
    path: reportDatabaseSourcePath,
    query: {
        engine: "sqlite",
        language: "sql",
        executed_at: generatedAt,
        description,
        sql: datasetQueries[dataset],
        tables_used: [dataset],
        filters: [
            "7 technical datasets",
            "GRCh38 chr1–22",
            "HCC1395 × HCC1395_DORADO is the preregistered same-ID pair",
        ],
        metric_definitions: metricDefinitions,
    },
});

const sources = [
    sqlSource(
        "source_specificity",
        "Executed SQLite snapshot query: specificity",
        "specificity",
        "Target rank and next-unrelated Jaccard across seven strict evidence layers.",
        [
            "Jaccard = intersection / union within each evidence layer.",
            "Target rank is among all 21 unordered pairs.",
        ],
    ),
    sqlSource(
        "source_specificity_chart",
        "Executed SQLite snapshot query: specificity_chart",
        "specificity_chart",
        "Tidy same-ID versus next-unrelated Jaccard rows for the grouped chart.",
    ),
    sqlSource(
        "source_shared_structure",
        "Executed SQLite snapshot query: shared_structure",
        "shared_structure",
        "Shared-candidate node, edge, and co-membership counts and containment.",
        [
            "DORADO→HCC containment = intersection / DORADO set size.",
            "Read-linkage edges are undirected and do not encode ancestry.",
        ],
    ),
    sqlSource(
        "source_shared_structure_chart",
        "Executed SQLite snapshot query: shared_structure_chart",
        "shared_structure_chart",
        "Tidy symmetric Jaccard and directional containment rows.",
    ),
    sqlSource(
        "source_resolution",
        "Executed SQLite snapshot query: resolution",
        "resolution",
        "HCC1395 and DORADO phase-fragmentation and linkage-reach descriptors.",
        ["DORADO/HCC ratios are computed within each metric and never across mixed units."],
    ),
    sqlSource(
        "source_denominator",
        "Executed SQLite snapshot query: evidence_denominators",
        "evidence_denominators",
        "Parent-universe coverage paired with conditional agreement.",
        [
            "Joint-W coverage = both-source W alleles / shared candidate alleles.",
            "Direction coverage = both-source determinate pairs / read-directed relation pairs.",
            "Multi-cluster coverage = informative both-multicluster regions / fixed regions.",
        ],
    ),
    sqlSource(
        "source_denominator_chart",
        "Executed SQLite snapshot query: denominator_chart",
        "denominator_chart",
        "Tidy conditional-agreement and informative-coverage rows.",
    ),
    sqlSource(
        "source_evidence_ladder",
        "Executed SQLite snapshot query: evidence_ladder",
        "evidence_ladder",
        "Evidence ladder from exact allele identity through global directed-tree absence.",
        [
            "Raw VAF is not CCF.",
            "PyClone-VI clusters are not a phylogenetic tree.",
        ],
    ),
    sqlSource(
        "source_seqc2_comparison",
        "Executed SQLite snapshot query: seqc2_comparison",
        "seqc2_comparison",
        "Valid and invalid comparison targets between SEQC2 Figure 4 and current InterSubMod evidence.",
    ),
    sqlSource(
        "source_fusion_gates",
        "Executed SQLite snapshot query: fusion_gates",
        "fusion_gates",
        "Sequential evidence gates from identity through orthogonal lineage truth.",
    ),
    sqlSource(
        "source_tree_metrics",
        "Executed SQLite snapshot query: tree_metrics",
        "tree_metrics",
        "Tree-to-tree and candidate-set comparison metrics with required coverage denominators.",
    ),
    sqlSource(
        "source_claim_tiers",
        "Executed SQLite snapshot query: claim_tiers",
        "claim_tiers",
        "Claim ceiling by evidence level.",
    ),
    {
        id: "input_strict_pairwise",
        label: "Input: HCC1395 cross-source strict pairwise metrics",
        path: "InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/strict_pair_validation/results/strict_pairwise_metrics_01.json",
    },
    {
        id: "input_exact_region_summary",
        label: "Input: exact PS×HP strict read-linkage report artifact",
        path: "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/artifact.json",
    },
    {
        id: "input_raw_vaf",
        label: "Input: cross-source raw VAF and PyClone report artifact",
        path: "InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/artifact.json",
    },
    {
        id: "input_clone_region_bridge",
        label: "Input: independent PyClone clone-region bridge",
        path: "InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/summary.json",
    },
    {
        id: "seqc2_paper",
        label: "Fang et al. 2021 SEQC2 HCC1395 reference study",
        href: "https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/",
    },
    {
        id: "seqc2_figure4",
        label: "Official NCBI viewer: SEQC2 HCC1395 Figure 4",
        href: "https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click%20on%20image%20to%20zoom&p=PMC3&id=8532138_nihms-1740806-f0004.jpg",
    },
    {
        id: "seqc2_truth_readme",
        label: "Official SEQC2 v1.2.1 truth package README",
        href: "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2.1/README.md",
    },
    {
        id: "seqc2_cnv",
        label: "Masood et al. 2024 HCC1395 CNV/LOH benchmark",
        href: "https://pmc.ncbi.nlm.nih.gov/articles/PMC11188507/",
    },
];

const charts = [
    {
        id: "specificity_chart",
        dataset: "specificity_chart",
        type: "bar",
        title: "HCC same-ID pair 與次佳不同 biological ID 的 Jaccard",
        subtitle: "HCC1395 × DORADO 在七個證據層皆排名 1/21；不同 cell line 的重疊接近零。",
        question: "同一 biological ID 的相似度是否顯著高於不同樣本？",
        rationale: "Grouped bars retain the same Jaccard scale and show the negative-control margin.",
        intent: "compare",
        valueFormat: "percent",
        sourceId: "source_specificity_chart",
        encodings: {
            x: { field: "layer", type: "nominal", label: "證據層" },
            y: { field: "jaccard", type: "quantitative", label: "Jaccard" },
            color: { field: "pair_class", type: "nominal", label: "配對類型" },
            tooltip: [
                { field: "comparator_pair", type: "nominal", label: "Pair" },
                { field: "jaccard", type: "quantitative", label: "Jaccard" },
            ],
        },
        options: { orientation: "vertical", grouping: "grouped" },
    },
    {
        id: "shared_structure_chart",
        dataset: "shared_structure_chart",
        type: "bar",
        title: "共同位點母群中的對稱重疊與 DORADO→HCC containment",
        subtitle: "Containment 遠高於 Jaccard；結果與 DORADO 為較稀疏 subset/contraction 相容，但不證明有向拓撲。",
        question: "低 Jaccard 是錯邊，還是解析量不對稱？",
        rationale: "Grouped bars separate symmetric overlap from directional containment.",
        intent: "compare",
        valueFormat: "percent",
        sourceId: "source_shared_structure_chart",
        encodings: {
            x: { field: "layer", type: "nominal", label: "共同母群證據層" },
            y: { field: "value", type: "quantitative", label: "比例" },
            color: { field: "measure", type: "nominal", label: "比較量" },
            tooltip: [
                { field: "numerator", type: "quantitative", label: "Intersection" },
                { field: "denominator_note", type: "nominal", label: "Denominator" },
            ],
        },
        options: { orientation: "vertical", grouping: "grouped" },
    },
    {
        id: "resolution_chart",
        dataset: "resolution",
        type: "bar",
        title: "DORADO / HCC1395 的 phase 與 linkage 解析比例",
        subtitle: "DORADO 有更多 PS boundaries，但 W span 與 direct-edge reach 較短；molecule rows 並未較少。",
        question: "兩來源的結構差異是否可由 phase fragmentation 與 linkage reach 解釋？",
        rationale: "A normalized ratio makes mixed-unit resolution descriptors comparable without mixing raw scales.",
        intent: "compare",
        sourceId: "source_resolution",
        encodings: {
            x: { field: "metric", type: "nominal", label: "解析指標" },
            y: { field: "dorado_over_hcc", type: "quantitative", label: "DORADO / HCC1395" },
            tooltip: [
                { field: "hcc", type: "quantitative", label: "HCC1395" },
                { field: "dorado", type: "quantitative", label: "DORADO" },
                { field: "unit", type: "nominal", label: "Unit" },
                { field: "interpretation", type: "nominal", label: "Interpretation" },
            ],
        },
        options: { orientation: "vertical", grouping: "single" },
    },
    {
        id: "denominator_chart",
        dataset: "denominator_chart",
        type: "bar",
        title: "條件式一致率與可辨識覆蓋率",
        subtitle: "98–99% 的局部一致同時伴隨 7–13% 的全域可辨識覆蓋；兩者必須一起解讀。",
        question: "高一致率是否涵蓋足夠多位點，可外推到 global fusion tree？",
        rationale: "Grouped bars prevent conditional fidelity from being mistaken for global coverage.",
        intent: "compare",
        valueFormat: "percent",
        sourceId: "source_denominator_chart",
        encodings: {
            x: { field: "evidence", type: "nominal", label: "證據終點" },
            y: { field: "value", type: "quantitative", label: "比例" },
            color: { field: "measure", type: "nominal", label: "量測" },
            tooltip: [
                { field: "numerator", type: "quantitative", label: "Informative n" },
                { field: "denominator", type: "quantitative", label: "Parent universe n" },
            ],
        },
        options: { orientation: "vertical", grouping: "grouped" },
    },
];

const tables = [
    {
        id: "specificity_table",
        dataset: "specificity",
        title: "21 組配對的 target rank 與次佳 unrelated 對照",
        subtitle: "所有比例使用各證據層自己的 set universe；ratio 只作 margin 描述。",
        sourceId: "source_specificity",
        density: "compact",
        defaultSort: { field: "order", direction: "asc" },
        columns: [
            { field: "order", label: "#", format: "number" },
            { field: "layer", label: "Evidence layer", type: "text" },
            { field: "target_jaccard_pct", label: "HCC pair Jaccard (%)", format: "number" },
            { field: "target_rank", label: "Rank / 21", format: "number" },
            { field: "next_unrelated_pair", label: "Next unrelated pair", type: "text" },
            { field: "next_unrelated_jaccard_pct", label: "Next Jaccard (%)", format: "number" },
            { field: "target_to_next_ratio", label: "Target / next", format: "number" },
        ],
    },
    {
        id: "shared_structure_table",
        dataset: "shared_structure",
        title: "共同位點母群中的 node、edge 與 co-membership 明細",
        subtitle: "Candidate row 使用完整 strict catalog；其餘投影到 76,721 shared candidates。",
        sourceId: "source_shared_structure",
        density: "compact",
        defaultSort: { field: "order", direction: "asc" },
        columns: [
            { field: "order", label: "#", format: "number" },
            { field: "layer", label: "Layer", type: "text" },
            { field: "unit", label: "Unit", type: "text" },
            { field: "hcc_n", label: "HCC1395 n", format: "number" },
            { field: "dorado_n", label: "DORADO n", format: "number" },
            { field: "intersection_n", label: "Intersection n", format: "number" },
            { field: "jaccard_pct", label: "Jaccard (%)", format: "number" },
            { field: "dorado_in_hcc_pct", label: "DOR→HCC (%)", format: "number" },
            { field: "universe", label: "Universe", type: "text" },
        ],
    },
    {
        id: "resolution_table",
        dataset: "resolution",
        title: "HCC1395 與 DORADO 的 phase / linkage 解析差異",
        subtitle: "不同單位不在原始尺度直接相加；ratio 只表示同一指標內兩來源的相對量。",
        sourceId: "source_resolution",
        density: "compact",
        defaultSort: { field: "order", direction: "asc" },
        columns: [
            { field: "order", label: "#", format: "number" },
            { field: "metric", label: "Metric", type: "text" },
            { field: "unit", label: "Unit", type: "text" },
            { field: "hcc", label: "HCC1395", format: "number" },
            { field: "dorado", label: "DORADO", format: "number" },
            { field: "dorado_over_hcc", label: "DOR / HCC", format: "number" },
            { field: "interpretation", label: "Interpretation", type: "text" },
        ],
    },
    {
        id: "denominator_table",
        dataset: "evidence_denominators",
        title: "局部一致率的 parent universe 與 informative coverage",
        subtitle: "每個條件式百分比均附上其真正可評估分母。",
        sourceId: "source_denominator",
        density: "compact",
        defaultSort: { field: "order", direction: "asc" },
        columns: [
            { field: "order", label: "#", format: "number" },
            { field: "evidence", label: "Evidence", type: "text" },
            { field: "parent_universe", label: "Parent universe", type: "text" },
            { field: "parent_n", label: "Parent n", format: "number" },
            { field: "informative_n", label: "Informative n", format: "number" },
            { field: "informative_coverage_pct", label: "Coverage (%)", format: "number" },
            { field: "conditional_agreement_pct", label: "Conditional agreement (%)", format: "number" },
            { field: "agreement_definition", label: "Agreement definition", type: "text" },
            { field: "limitation", label: "Limit", type: "text" },
        ],
    },
    {
        id: "evidence_ladder_table",
        dataset: "evidence_ladder",
        title: "從位點到全域演化樹的證據階梯",
        subtitle: "各列不是同一 estimand；用途是顯示隨推論層級增加而擴大的 identifiability gap。",
        sourceId: "source_evidence_ladder",
        density: "compact",
        defaultSort: { field: "order", direction: "asc" },
        columns: [
            { field: "order", label: "#", format: "number" },
            { field: "layer", label: "Layer", type: "text" },
            { field: "primary_metric", label: "Primary metric", type: "text" },
            { field: "estimate_text", label: "Estimate", type: "text" },
            { field: "interpretation", label: "Interpretation", type: "text" },
            { field: "claim_level", label: "Claim level", type: "text" },
        ],
    },
    {
        id: "seqc2_table",
        dataset: "seqc2_comparison",
        title: "SEQC2 Figure 4 與 InterSubMod 的有效對照方式",
        subtitle: "對照生物學層級與可量化關係，不以畫面位置或枝長作肉眼驗證。",
        sourceId: "source_seqc2_comparison",
        density: "compact",
        defaultSort: { field: "order", direction: "asc" },
        columns: [
            { field: "order", label: "#", format: "number" },
            { field: "item", label: "Comparison item", type: "text" },
            { field: "seqc2_reference", label: "SEQC2 reference", type: "text" },
            { field: "intersubmod_current", label: "InterSubMod current", type: "text" },
            { field: "valid_comparison", label: "Valid comparison", type: "text" },
            { field: "verdict", label: "Verdict", type: "text" },
        ],
    },
    {
        id: "fusion_gate_table",
        dataset: "fusion_gates",
        title: "從局部骨架升格為融合 clone tree 的 gates",
        subtitle: "PASS 只授權該列 claim；前一 gate 通過不會自動使後一 gate 通過。",
        sourceId: "source_fusion_gates",
        density: "compact",
        defaultSort: { field: "gate_order", direction: "asc" },
        columns: [
            { field: "gate_order", label: "#", format: "number" },
            { field: "gate", label: "Gate", type: "text" },
            { field: "required_evidence", label: "Required evidence", type: "text" },
            { field: "current_value", label: "Current value", type: "text" },
            { field: "status", label: "Status", type: "text" },
            { field: "claim_if_pass", label: "Claim unlocked", type: "text" },
        ],
    },
    {
        id: "tree_metric_table",
        dataset: "tree_metrics",
        title: "Strict tree-to-tree 與 candidate-set 比較指標",
        subtitle: "所有一致度必須同報 jointly evaluable coverage；先 collapse 低支持 edge 再比較。",
        sourceId: "source_tree_metrics",
        density: "compact",
        defaultSort: { field: "order", direction: "asc" },
        columns: [
            { field: "order", label: "#", format: "number" },
            { field: "object", label: "Object", type: "text" },
            { field: "primary_metric", label: "Primary metric", type: "text" },
            { field: "must_report", label: "Must co-report", type: "text" },
            { field: "fusion_rule", label: "Fusion rule", type: "text" },
        ],
    },
    {
        id: "claim_tier_table",
        dataset: "claim_tiers",
        title: "證據層級與可主張範圍",
        subtitle: "L1–L4 目前各自支持不同層級；L5 exact tree truth 仍缺。",
        sourceId: "source_claim_tiers",
        density: "compact",
        defaultSort: { field: "level", direction: "asc" },
        columns: [
            { field: "level", label: "Level", type: "text" },
            { field: "evidence", label: "Evidence", type: "text" },
            { field: "current_result", label: "Current result", type: "text" },
            { field: "allowed_claim", label: "Allowed claim", type: "text" },
        ],
    },
];

const blocks = [
    {
        id: "title",
        type: "markdown",
        body: `# ${title}\n\nTask Type B comprehensive validation · 7 technical datasets × chr1–22 · 21 pairwise comparisons · same-ID positive pair = HCC1395 × HCC1395_DORADO · SEQC2 Figure 4 僅作 external reference architecture`,
    },
    {
        id: "technical_summary",
        type: "markdown",
        body: "## 技術摘要\n\n**結論：在 exact shared-allele common universe、bulk raw-VAF profile 與 phase-invariant undirected strict-region endpoints 上，HCC1395 × HCC1395_DORADO 具有很強的 biological-ID specificity；兩者都遠比任何不同 cell line 相近。這與 gross biological-ID mismatch 或 all-pairs homogenization 不相容，但不能排除較細微的 processing bias，也不驗證 directed clone topology。**（影響：高；信心：底層位點與局部無向連結高、global ancestry 低）\n\n最強證據是 candidate Jaccard 92.76%、raw VAF CCC 0.9339，以及共同 candidate 母群內 DORADO→HCC edge containment 98.28%；而 20 組不同 biological ID 的 next-best edge Jaccard 僅 0.146%。相反地，joint-W coverage 只有 13.23%、兩來源皆可判定方向只有 7.39%，independent PyClone minor-set Jaccard 只有 0.381。正確敘述因此是：**candidate-site 與 bulk raw-VAF backbone 高度相似；shared-opportunity 中可辨識的局部無向連結高度非衝突。精細 minor branch、clone fraction 與 parent→child edge 尚未高度重現或完成 strict 驗證。**",
    },
    {
        id: "specificity_heading",
        type: "markdown",
        body: "## 只有 HCC1395 same-ID pair 明顯相似\n\nHCC pair 在 candidate、active、W、direct edge、ALT edge、AA edge 與 exact component 七層都排名 1/21。結果與「演算法把所有配對都做成相似結構」或 gross sample mix-up 不相容，但不能排除較細微的 processing bias。這仍只有一組 positive biological ID，屬很強的個案 specificity，而不是多 cell-line replicate 的普遍敏感度估計。",
        sourceId: "source_specificity",
    },
    {
        id: "specificity_chart_block",
        type: "chart",
        chartId: "specificity_chart",
    },
    {
        id: "specificity_table_block",
        type: "table",
        tableId: "specificity_table",
        layout: "full",
    },
    {
        id: "skeleton_heading",
        type: "markdown",
        body: "## DORADO 局部訊號與 HCC1395 的低解析 subset / contraction 解釋相容\n\n限制在 76,721 個共同 candidate alleles 後，DORADO→HCC 的 active、W、direct-edge、ALT-edge 與 co-membership containment 分別為 99.48%、97.45%、98.28%、98.16% 與 90.05%。這表示 DORADO 已觀測到的局部關係多數可在 HCC 找到，且 HCC 額外解析出更多 edge／component；結果與較稀疏 subset / contraction 解釋相容，但不證明兩者有相同有向拓撲。這一層只比較 phase-invariant、shared-opportunity 上的無向 read-linkage，不是 evolutionary parent→child edge。",
        sourceId: "source_shared_structure",
    },
    {
        id: "shared_structure_chart_block",
        type: "chart",
        chartId: "shared_structure_chart",
    },
    {
        id: "shared_structure_table_block",
        type: "table",
        tableId: "shared_structure_table",
        layout: "full",
    },
    {
        id: "resolution_heading",
        type: "markdown",
        body: "## Phase fragmentation 與解析差異相容，是主要候選貢獻因子之一\n\nDORADO 的 exact chromosome×PS 數是 HCC 的 2.90 倍，median W span 與 median direct-edge distance 卻只有 HCC 的約 32% 與 40%。Canonical molecule rows 反而略多，因此差異不能簡化為總 read evidence 不足；phase block 切碎、read-linkage reach、endpoint observability 與 edge identifiability 是合理的主要候選因子。未完成 matched downsampling、PS/opportunity matching 與 ablation 前，不能量化 phase fragmentation 的因果占比。",
        sourceId: "source_resolution",
    },
    {
        id: "resolution_chart_block",
        type: "chart",
        chartId: "resolution_chart",
    },
    {
        id: "resolution_table_block",
        type: "table",
        tableId: "resolution_table",
        layout: "full",
    },
    {
        id: "denominator_heading",
        type: "markdown",
        body: "## 高條件式一致不能取代全域可辨識覆蓋\n\nDORADO→HCC edge containment 98.28% 與兩來源同向 99.67% 都是真的，但它們回答的是「已被兩側觀測到的子集合是否衝突」。共同 candidate 中只有 13.23% 同時進 W；8,096 個 read-directed relation pairs 中只有 598 個由兩側共同判定方向。融合樹若不把 coverage 與 agreement 一起報告，會把 conditional fidelity 誤寫成 global topology accuracy。",
        sourceId: "source_denominator",
    },
    {
        id: "denominator_chart_block",
        type: "chart",
        chartId: "denominator_chart",
    },
    {
        id: "denominator_table_block",
        type: "table",
        tableId: "denominator_table",
        layout: "full",
    },
    {
        id: "vaf_minor_heading",
        type: "markdown",
        body: "## Bulk raw-VAF profile 高度相近，但 minor subclone membership 只有中度\n\nExact shared alleles 的 raw VAF Pearson=0.9343、CCC=0.9339、MAE=0.0624，且已用 chromosome+1 Mb block bootstrap 計算不確定度，支持 bulk allele-frequency profile 高度相似；VAF 高相關本身不能直接等同演化 trunk。兩側 independent PyClone-VI fits 的 ARI=0.539、clonal/subclonal κ=0.544、minor-set Jaccard=0.381；真正兩側皆 multi-cluster 的 regions 只有 34 個，其中 21 個 partition exact（61.76%）。因此不能把整體 clonal-majority agreement 當成精細 subclone 已穩定。",
        sourceId: "source_evidence_ladder",
    },
    {
        id: "evidence_ladder_table_block",
        type: "table",
        tableId: "evidence_ladder_table",
        layout: "full",
    },
    {
        id: "seqc2_heading",
        type: "markdown",
        body: "## SEQC2 Figure 4 支持 HCC1395 存在 branching heterogeneity，但不是本地融合樹的唯一真值\n\nFang et al. 以 bulk WGS/WES 的 PhyloWGS 推論 Figure 4a：S1 是 MRCA，兩個主要分支從 S2 與 S8 開始；S3、S6 因 CCF<10% 未畫。SuperFreq、subHMM 與 10x single-cell CNV（1,270 HCC1395 tumor cells；638 HCC1395BL cells）支持存在 clonal/subclonal CNA 與異質性，但沒有逐邊確認 S1–S10 的每條 ancestry edge。因此 Figure 4 適合作為輸出層級與粗粒度生物學 reference，不應當作可直接比對的 exact tree truth。\n\n[開啟 SEQC2 主論文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/) · [開啟官方 Figure 4](https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click%20on%20image%20to%20zoom&p=PMC3&id=8532138_nihms-1740806-f0004.jpg)",
        sourceId: "seqc2_paper",
    },
    {
        id: "seqc2_schematic_block",
        type: "html",
        body: seqc2Schematic,
    },
    {
        id: "seqc2_table_block",
        type: "table",
        tableId: "seqc2_table",
        layout: "full",
    },
    {
        id: "fusion_method_heading",
        type: "markdown",
        body: "## 多位點融合樹可以合理做到，但輸出應先是 candidate ensemble 與 collapsed consensus\n\n兩來源的 edges 不可直接聯集；那會把單側 technical edge 製造成 ancestry，甚至引入 cycle。合理做法是兩側各自獨立跑 strict-bound solver，投影到共同 exact allele universe，轉成 ancestor / descendant / same-clone / branched / unresolved 關係，再用 CN/purity/multiplicity-aware CCF 排序。雙側同向才進 strict consensus；單側可觀測用虛線；反向列 conflict；未解關係收縮成 polytomy。只有 candidate set 唯一且 bootstrap 穩定時，才輸出單一精細樹。",
    },
    {
        id: "fusion_flow_block",
        type: "html",
        body: fusionFlow,
    },
    {
        id: "fusion_gate_table_block",
        type: "table",
        tableId: "fusion_gate_table",
        layout: "full",
    },
    {
        id: "tree_metrics_heading",
        type: "markdown",
        body: "## 正式比較應以 candidate-tree set、partial order 與 contraction compatibility 為主\n\n同一 biological ID 的正式 endpoint 不應是任選一棵 best tree 的畫面相似，而是兩側候選集合是否至少存在相容樹、共同 leaves 上的 ancestor relation 與 rooted triplets 是否一致、低支持 edges 收縮後 clades 是否相容，以及 HCC pair 是否在 coverage-matched、PS-matched 的 20 組 negative controls 中仍排名第一。所有統計以 chromosome×PS block bootstrap 估計不確定度。",
    },
    {
        id: "tree_metric_table_block",
        type: "table",
        tableId: "tree_metric_table",
        layout: "full",
    },
    {
        id: "method_verdict_heading",
        type: "markdown",
        body: "## 方法學判定：strict-region contract 與 ID specificity 有支持，clone-tree 尚未驗證\n\n可以說 strict-region receipts / contract 已通過，且 common-universe 比較保留很強的 biological-ID signal；HCC 與 DORADO 的 pairwise receipts 使用同一 pre-hotfix builder，可降低兩者間的 code-version confounding。這些證據與 gross sample mix-up、all-pairs homogenization 或大量局部無向 edge 衝突不相容，但不能據此說演算法沒有重大問題，也不能排除較細微的 processing bias。\n\n跨區融合、clone 數、minor membership、clone fraction、rooted branching 與 exact parent→child edge 尚未過關。DORADO strict partition / candidate-tree set 尚未完成，shared-CN 也不是 source-specific CN measurement。故目前最強且正確的結論是：**HCC1395 兩來源的 candidate-site 與 bulk raw-VAF profile 高度相似；在目前可辨識的 shared-opportunity 無向局部連結中，DORADO 訊號大多能在 HCC 找到。完整融合樹相同仍是假說，而非已驗證結果。**",
    },
    {
        id: "claim_tier_table_block",
        type: "table",
        tableId: "claim_tier_table",
        layout: "full",
    },
    {
        id: "limitations_heading",
        type: "markdown",
        body: "## 限制、不確定度與 robustness checks\n\n- biological positive pair 只有 HCC1395 × DORADO 一組；其餘 20 組是不同 cell lines，不能估計普遍的 same-ID sensitivity。\n- 同一 cell line 不保證相同 aliquot、passage、library 或 molecules；殘差不能全歸因於演算法。\n- HCC 與 DORADO pairwise receipts 同用 pre-hotfix builder `7260a763…`、graph core `df3d6f…` 與 extractor `2ca7ccb…`，可降低 pair 內版本混雜；但 all-7 cohort 並非完全同 builder，HCC1954 使用 post-hotfix `912721f9…`。Data-specific no-trigger equivalence 支持本次數值未受該修補影響，不等於演算法已正確。\n- Exact-PS boundary 是合理的 fail-closed contract，但跨 pipeline PS fragmentation 會降低 region／edge recall。\n- Joint-W、direction 與 multi-cluster analyses 都有強 selection；高條件式一致不能外推到未觀測區。\n- DORADO 使用共享 HCC1395 CN 只能作 sensitivity，不是獨立 source-specific CN corroboration。\n- PyClone-VI 產生 CCF clusters，不產 phylogenetic tree，也不是 gold truth。\n- SEQC2 Figure 4 是 inferred reference architecture；single-cell CNV 支持 heterogeneity，不是每條 directed edge 的 truth set。\n- 下一輪必做：DORADO strict-bound solver、common-universe tree-set comparison、coverage/PS/read-length matching、block bootstrap 與 orthogonal CN/lineage validation。",
    },
    {
        id: "reproducibility_heading",
        type: "markdown",
        body: "## 資料來源與重現方式\n\n本報告依各 artifact 的既定 claim ceiling 重組，沒有重新跑長時間 topology solver，也未將 conditional / PARTIAL evidence 升級為 topology validation。四個輸入以固定 SHA-256 fail-closed 綁定；HTML 由同一 canonical artifact 產生。外部資料只引用 PMC/NCBI 與官方 SEQC2 truth package，不把網路圖片當成本地結果。\n\n完整重建與 HTML 封裝命令記錄於 companion source note；portable builder 的驗證 receipt 亦回寫至該文件。",
    },
    {
        id: "further_questions_heading",
        type: "markdown",
        body: "## 尚待回答\n\n- DORADO strict-bound candidate-tree set 完成後，collapsed rooted-triplet 與 ancestor-pair concordance 是否仍在 21 pairs 中排名第一？\n- 在 source-specific allele CN 下，major-clone CCF CCC 與 minor membership 是否改善？\n- 若將 HCC read evidence 稀釋／切碎到 DORADO 的 PS、span、edge-opportunity 分佈，Jaccard 與 containment 差距剩多少？\n- 目前 2/598 opposite-direction pairs 是 mapping/phase artifact、模型不確定，還是真實 biological drift？\n- 哪些 high-support local relations 能由 SEQC2 single-cell CNV 或獨立多區域資料直接檢驗？",
    },
];

const rawDatasets = {
    headline,
    specificity: specificityRows,
    specificity_chart: specificityChartRows,
    shared_structure: sharedRows,
    shared_structure_chart: sharedChartRows,
    resolution: resolutionRows,
    evidence_denominators: evidenceDenominatorRows,
    denominator_chart: denominatorChartRows,
    evidence_ladder: evidenceLadderRows,
    seqc2_comparison: seqc2ComparisonRows,
    fusion_gates: fusionGateRows,
    tree_metrics: treeMetricRows,
    claim_tiers: claimTierRows,
};

const quoteIdentifier = (identifier) => `"${String(identifier).replaceAll('"', '""')}"`;
const sqliteTypeForColumn = (rows, column) => {
    const values = rows.map((row) => row[column]).filter((value) => value !== null && value !== undefined);
    if (values.length > 0 && values.every((value) => typeof value === "number")) {
        return "REAL";
    }
    if (values.length > 0 && values.every((value) => typeof value === "boolean")) {
        return "INTEGER";
    }
    return "TEXT";
};

mkdirSync(dirname(reportDatabasePath), { recursive: true });
const reportDatabase = new DatabaseSync(reportDatabasePath);
const reviewedDatasets = {};

try {
    for (const [dataset, rows] of Object.entries(rawDatasets)) {
        if (rows.length === 0) {
            throw new Error(`Dataset ${dataset} has no rows.`);
        }

        const columns = [...new Set(rows.flatMap((row) => Object.keys(row)))];
        const tableName = quoteIdentifier(dataset);
        const columnDefinitions = columns
            .map((column) => `${quoteIdentifier(column)} ${sqliteTypeForColumn(rows, column)}`)
            .join(", ");

        reportDatabase.exec(`DROP TABLE IF EXISTS ${tableName};`);
        reportDatabase.exec(`CREATE TABLE ${tableName} (${columnDefinitions});`);

        const placeholders = columns.map(() => "?").join(", ");
        const insert = reportDatabase.prepare(
            `INSERT INTO ${tableName} (${columns.map(quoteIdentifier).join(", ")}) VALUES (${placeholders});`,
        );

        reportDatabase.exec("BEGIN;");
        try {
            for (const row of rows) {
                insert.run(
                    ...columns.map((column) => {
                        const value = row[column];
                        if (value === null || value === undefined) {
                            return null;
                        }
                        return typeof value === "boolean" ? Number(value) : value;
                    }),
                );
            }
            reportDatabase.exec("COMMIT;");
        } catch (error) {
            reportDatabase.exec("ROLLBACK;");
            throw error;
        }

        reviewedDatasets[dataset] = reportDatabase.prepare(datasetQueries[dataset]).all();
    }
} finally {
    reportDatabase.close();
}

const artifact = {
    surface: "report",
    manifest: {
        version: 1,
        surface: "report",
        title,
        description:
            "Technical validation of HCC1395 biological-ID specificity and the evidence gates required for a multi-locus clone/subclone fusion tree.",
        generatedAt,
        blocks,
        charts,
        tables,
        sources,
    },
    snapshot: {
        version: 1,
        status: "ready",
        generatedAt,
        datasets: reviewedDatasets,
    },
    sources,
};

mkdirSync(dirname(outputPath), { recursive: true });
writeFileSync(outputPath, `${JSON.stringify(artifact, null, 2)}\n`, "utf8");

const outputSummary = {
    artifact_path: outputPath,
    sqlite_path: reportDatabasePath,
    all_strict_validations_pass: strict.all_validations_pass,
    blocks: blocks.length,
    charts: charts.length,
    tables: tables.length,
    sources: sources.length,
    datasets: Object.keys(artifact.snapshot.datasets).length,
    executed_queries: Object.keys(datasetQueries).length,
    snapshot_rows: Object.values(artifact.snapshot.datasets).reduce(
        (sum, rows) => sum + rows.length,
        0,
    ),
    headline: headline[0],
};

process.stdout.write(`${JSON.stringify(outputSummary, null, 2)}\n`);
