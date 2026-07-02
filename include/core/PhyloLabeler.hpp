#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "TreeStructure.hpp"

namespace InterSubMod {

/**
 * @brief Configuration for phylo-v4.1 hierarchical cluster labeling.
 *
 * Pinned defaults match the Python reference (phylo_v4.py) for cross-language
 * statistical equivalence. The seed is fixed at compile time for Hard-Gate
 * reproducibility (no std::random_device).
 */
struct PhyloConfig {
    int min_sz = 3;            ///< Minimum group / outlier size (MIN_SZ)
    double sep_min = 1.3;      ///< between/within ratio floor (SEP_MIN, safe-plateau convention)
    int rnull = 40;            ///< Null permutation replicates (RNULL)
    int min_common = 3;        ///< Minimum common CpG sites for BERNOULLI (C_min)
    double null_pct = 95.0;    ///< Null percentile gate (coarse=95; fine=90 in B2)
    uint64_t seed = 20260622;  ///< Pinned base seed (Hard-Gate deterministic)
};

/**
 * @brief Per-locus phylo labeling result (Phase B1: coarse only).
 */
struct PhyloResult {
    std::vector<std::string> labels;  ///< Per-read hierarchical label ("1","1-1",...) or "outlier"
    int n_groups = 0;                 ///< Number of distinct non-outlier terminal labels
    int n_outlier = 0;                ///< Number of reads labeled "outlier"
    bool valid = false;               ///< False if too few reads
};

/**
 * @brief phylo-v4.1 hierarchical cluster labeler (C++ port of phylo_v4.py v4_label).
 *
 * Phase B1 (coarse): descend-quarantine to a balanced split node, gate the split by
 * a per-subgroup column-permutation BERNOULLI null (null95), and emit recursive
 * hierarchical labels (larger child = lower number). Replaces the silhouette-based
 * TreeCutter::find_optimal_clusters path (quantified 75-91% pure-noise FP).
 *
 * Phases B2/B3 (modal instability, fine layer, other group, hidden-het, alignment,
 * RNG verification harness) extend this class.
 */
class PhyloLabeler {
public:
    explicit PhyloLabeler(const PhyloConfig& cfg = PhyloConfig());

    /**
     * @brief Label reads into a hierarchical cluster partition.
     *
     * @param dist Symmetric NxN read-read BERNOULLI distance (already peeled; -1 = invalid pair).
     * @param meth NxC methylation matrix aligned with dist rows (NaN = missing CpG).
     * @return Coarse labeling result.
     */
    PhyloResult label(const Eigen::MatrixXd& dist, const Eigen::MatrixXd& meth);

    /**
     * @brief Recompute BERNOULLI distance from a methylation matrix.
     *
     * Matches DistanceMatrix::calculate_bernoulli (weight=2|p-0.5|, delta=p(1-q)+(1-p)q,
     * Dist=sum(w*delta)/sum(w)); invalid pairs (common<min_sz or zero weight) set to -1.
     * Exposed for the per-subgroup null permutation.
     */
    Eigen::MatrixXd bernoulli_dist(const Eigen::MatrixXd& meth) const;

private:
    PhyloConfig cfg_;

    /// between(S1,S2)/within(S1∪S2) mean BERNOULLI ratio (sibling separation); <0 dropped.
    double between_within(const Eigen::MatrixXd& s, const std::vector<int>& S1, const std::vector<int>& S2) const;

    /// 95-percentile of per-subgroup column-permutation null root-split ratio (sub-tree of `leaves`).
    double subgroup_null_pct(const Eigen::MatrixXd& meth, const std::vector<int>& leaves, uint64_t seed) const;

    /// Build a UPGMA tree from a distance matrix and return desc/ch maps keyed by node id.
    void build_desc_ch(const Eigen::MatrixXd& dist, std::vector<std::vector<int>>& desc,
                       std::vector<std::pair<int, int>>& ch, int& root, int n) const;
};

}  // namespace InterSubMod
