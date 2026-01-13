/**
 * @file SignificanceAnalyzer.cpp
 * @brief Implementation of main significance analyzer
 */

#include "core/SignificanceAnalyzer.hpp"

#include <cmath>
#include <numeric>

namespace InterSubMod {

SignificanceAnalyzer::SignificanceAnalyzer(const SignificanceConfig& config) : config_(config) {
    init_analyzers();
}

void SignificanceAnalyzer::init_analyzers() {
    if (config_.seed != 0) {
        config_.global_config.fisher_config.seed = config_.seed;
        config_.local_config.fisher_config.seed = config_.seed + 1;
        config_.structure_config.seed = config_.seed + 2;
        config_.bootstrap_config.seed = config_.seed + 3;
        config_.label_config.seed = config_.seed + 4;  // NEW: LabelTest seed
    }

    global_test_ = std::make_unique<GlobalTest>(config_.global_config);
    local_test_ = std::make_unique<LocalTest>(config_.local_config);
    structure_test_ = std::make_unique<StructureTest>(config_.structure_config);
    bootstrap_ = std::make_unique<Bootstrap>(config_.bootstrap_config);
    label_test_ = std::make_unique<LabelTest>(config_.label_config);  // NEW: Initialize LabelTest
}

void SignificanceAnalyzer::set_seed(uint64_t seed) {
    config_.seed = seed;
    init_analyzers();
}

SignificanceResult SignificanceAnalyzer::analyze(const std::vector<int>& cluster_labels,
                                                 const std::vector<FullLabel>& full_labels,
                                                 const Eigen::MatrixXd& dist_matrix,
                                                 const Eigen::MatrixXd& methylation_matrix, int region_id,
                                                 const std::string& anchor_key) {
    SignificanceResult result;

    // Set metadata
    result.run_id = config_.run_id;
    result.vcf_id = config_.vcf_id;
    result.bam_id = config_.bam_id;
    result.region_id = region_id;
    result.anchor_key = anchor_key;
    result.seed = config_.seed;

    // Basic info
    result.n_reads = static_cast<int>(cluster_labels.size());
    result.n_cpg = methylation_matrix.cols();

    if (result.n_reads < 2) {
        result.valid_flag = false;
        result.invalid_reason = "insufficient_reads";
        return result;
    }

    // Count labels
    int n_alt = 0, n_ref = 0;
    for (const auto& label : full_labels) {
        if (label.allele == AltSupport::ALT)
            ++n_alt;
        else if (label.allele == AltSupport::REF)
            ++n_ref;
    }

    // Phase 1: Global Tests
    if (config_.enable_global_test) {
        auto [alt_result, hp_result, sample_result] = global_test_->test_all(cluster_labels, full_labels);

        result.global_alt = alt_result;
        result.global_hp = hp_result;
        result.global_sample = sample_result;

        // Determine passed_gate based on ALT/REF or HP (any significant signal)
        result.passed_gate = alt_result.passed_gate || hp_result.passed_gate;
    }

    // Phase 2: Local Tests (only if passed gate)
    if (config_.enable_local_test && result.passed_gate) {
        result.local_result = local_test_->run(cluster_labels, full_labels);
    }

    // Phase 3: Structure Tests (only if passed gate)
    if (result.passed_gate) {
        // Filter reads for complete distance matrix
        std::vector<int> valid_read_indices;
        bool can_run_permanova = structure_test_->filter_reads_for_complete_matrix(dist_matrix, valid_read_indices);

        result.n_reads_valid = static_cast<int>(valid_read_indices.size());

        if (config_.enable_permanova && can_run_permanova) {
            // Build filtered distance matrix
            int n_valid = static_cast<int>(valid_read_indices.size());
            Eigen::MatrixXd filtered_dist(n_valid, n_valid);

            for (int i = 0; i < n_valid; ++i) {
                for (int j = 0; j < n_valid; ++j) {
                    filtered_dist(i, j) = dist_matrix(valid_read_indices[i], valid_read_indices[j]);
                }
            }

            // Get filtered labels
            std::vector<int> filtered_labels(n_valid);
            for (int i = 0; i < n_valid; ++i) {
                filtered_labels[i] = cluster_labels[valid_read_indices[i]];
            }

            result.permanova = structure_test_->run_permanova(filtered_dist, filtered_labels);

            if (config_.enable_dispersion) {
                result.dispersion = structure_test_->check_dispersion(filtered_dist, filtered_labels);
                result.dispersion_warning = result.dispersion.warning;
            }
        } else {
            result.permanova.valid = false;
            result.permanova.invalid_reason = "insufficient_overlap";
        }

        // Bootstrap (expensive, disabled by default)
        if (config_.enable_bootstrap) {
            // Simple re-clustering function based on hierarchical clustering
            // In practice, this would use the same clustering method as the original
            auto cluster_func = [](const Eigen::MatrixXd& dist) -> std::vector<int> {
                int n = dist.rows();
                std::vector<int> labels(n, 0);
                // Simple 2-cluster split based on average distance
                double total_mean = 0.0;
                int count = 0;
                for (int i = 0; i < n; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        if (!std::isnan(dist(i, j))) {
                            total_mean += dist(i, j);
                            ++count;
                        }
                    }
                }
                total_mean /= count;

                for (int i = 0; i < n; ++i) {
                    double sum = 0.0;
                    int cnt = 0;
                    for (int j = 0; j < n; ++j) {
                        if (i != j && !std::isnan(dist(i, j))) {
                            sum += dist(i, j);
                            ++cnt;
                        }
                    }
                    labels[i] = (cnt > 0 && sum / cnt > total_mean) ? 1 : 0;
                }
                return labels;
            };

            result.bootstrap = bootstrap_->run(methylation_matrix, cluster_labels, cluster_func);
        }
    }

    // Phase 4: Label-First Verification (NEW)
    // Run label-based delta tests to verify if HP/Allele labels explain distance structure
    if (dist_matrix.rows() >= config_.label_config.min_total_reads) {
        LabelTestResult label_result = label_test_->test_all(dist_matrix, full_labels);

        result.label_hp = label_result.hp_result;
        result.label_allele = label_result.allele_result;
        result.dominant_label_dimension = label_result.dominant_dimension;
        result.label_significant = label_result.any_significant;
    }

    // Phase 5: Classification (Bidirectional Verification)
    // Determine verification_class based on Label-First and Cluster-First results
    bool cluster_significant = result.passed_gate && 
                               (result.global_alt.fisher_ffh.p_value <= 0.05 || 
                                result.global_hp.fisher_ffh.p_value <= 0.05);
    bool label_sig = result.label_significant;

    if (label_sig && cluster_significant) {
        result.verification_class = "Strong";  // Strong association: Both paths agree
    } else if (!label_sig && cluster_significant) {
        result.verification_class = "Subclone";  // Subclone exists: Cluster structure without label correlation
    } else if (label_sig && !cluster_significant) {
        result.verification_class = "Weak";  // Weak association: Label signal but weak clustering
    } else {
        result.verification_class = "Noise";  // No association/noise: No significant signal
    }

    // Compute heuristic score
    result.heuristic_score = compute_heuristic_score(result);

    // Mark as ambiguous if dispersion warning or PERMANOVA invalid
    result.score_ambiguous = result.dispersion_warning || (config_.enable_permanova && !result.permanova.valid);

    result.valid_flag = true;
    return result;
}

SignificanceResult SignificanceAnalyzer::analyze_simple(const std::vector<int>& cluster_labels,
                                                        const std::vector<FullLabel>& full_labels, int region_id,
                                                        const std::string& anchor_key) {
    // Create dummy matrices
    Eigen::MatrixXd empty_dist(0, 0);
    Eigen::MatrixXd empty_meth(0, 0);

    // Temporarily disable structure tests
    bool orig_permanova = config_.enable_permanova;
    bool orig_dispersion = config_.enable_dispersion;
    bool orig_bootstrap = config_.enable_bootstrap;

    config_.enable_permanova = false;
    config_.enable_dispersion = false;
    config_.enable_bootstrap = false;

    SignificanceResult result = analyze(cluster_labels, full_labels, empty_dist, empty_meth, region_id, anchor_key);

    // Restore settings
    config_.enable_permanova = orig_permanova;
    config_.enable_dispersion = orig_dispersion;
    config_.enable_bootstrap = orig_bootstrap;

    return result;
}

double SignificanceAnalyzer::compute_heuristic_score(const SignificanceResult& result) {
    double score = 0.0;
    double p_val_alt = (result.global_alt.valid) ? result.global_alt.fisher_ffh.p_value : 1.0;
    double p_val_hp = (result.global_hp.valid) ? result.global_hp.fisher_ffh.p_value : 1.0;

    // Take the best (lowest) p-value from ALT or HP
    double best_p = std::min(p_val_alt, p_val_hp);

    // Primary: -log10(best_p)
    if (best_p > 0) {
        score = -std::log10(best_p);
    } else {
        score = 20.0;  // Cap for p=0
    }

    // Bonus for Cramér's V (best of ALT or HP)
    double best_v = 0.0;
    if (result.global_alt.cramers_v_reliable) best_v = std::max(best_v, result.global_alt.cramers_v);
    if (result.global_hp.cramers_v_reliable) best_v = std::max(best_v, result.global_hp.cramers_v);

    score += best_v * 2.0;

    // Penalty for dispersion warning
    if (result.dispersion_warning) {
        score *= 0.5;
    }

    // Penalty if PERMANOVA not significant
    if (result.permanova.valid && result.permanova.p_value > 0.05) {
        score *= 0.7;
    }

    return score;
}

}  // namespace InterSubMod
