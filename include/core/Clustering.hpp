#pragma once

#include <vector>
#include <map>
#include "Types.hpp"

namespace InterSubMod {

struct ClusterStats {
    int cluster_id;
    int size;
    int count_hp1;
    int count_hp2;
    int count_hp_unknown;
    int count_tumor;
    int count_normal;
    int count_alt;
    int count_ref;
    double p_value_hp;
    double p_value_somatic;
};

class ClusteringResult {
public:
    int region_id;
    int num_clusters;
    std::vector<int> labels; // labels[i] for read_ids[i]
    std::vector<double> silhouette_scores;
    std::map<int, ClusterStats> stats;
    
    // Placeholder for PhyloNode* phylo_root;
};

} // namespace InterSubMod

