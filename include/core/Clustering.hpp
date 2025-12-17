#pragma once

#include <map>
#include <vector>

#include "Types.hpp"

namespace InterSubMod {

/**
 * @brief Statistics for a single cluster of reads.
 */
struct ClusterStats {
    int cluster_id;
    int size;

    // Haplotype counts
    int count_hp1;
    int count_hp2;
    int count_hp_unknown;

    // Sample source counts
    int count_tumor;
    int count_normal;

    // Somatic allele counts
    int count_alt;
    int count_ref;

    // Statistical association (p-values)
    double p_value_hp;       ///< Association with Haplotype (Fisher's Exact)
    double p_value_somatic;  ///< Association with Somatic Allele (Fisher's Exact)
};

/**
 * @brief Result of clustering analysis for a single Region.
 */
class ClusteringResult {
public:
    int region_id;
    int num_clusters;
    std::vector<int> labels;  // labels[i] corresponds to read_ids[i] in distance matrix
    std::vector<double> silhouette_scores;
    std::map<int, ClusterStats> stats;

    // Placeholder for PhyloNode* phylo_root;
};

}  // namespace InterSubMod
