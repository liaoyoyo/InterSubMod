#pragma once

#include <memory>
#include <string>
#include <vector>

#include "core/BamReader.hpp"
#include "core/Config.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/HierarchicalClustering.hpp"
#include "core/MatrixBuilder.hpp"
#include "core/MethylationMatrix.hpp"
#include "core/MethylationParser.hpp"
#include "core/ReadParser.hpp"
#include "core/SignificanceAnalyzer.hpp"
#include "core/SomaticSnv.hpp"
#include "core/TreeStructure.hpp"
#include "io/RegionWriter.hpp"
#include "io/TreeWriter.hpp"
#include "utils/FastaReader.hpp"

namespace InterSubMod {

/**
 * @brief Statistics result for processing a single SNV region
 */
struct RegionResult {
    int region_id;
    int snv_id;
    int num_reads;
    int num_cpgs;
    int num_forward_reads;   ///< Forward strand reads
    int num_reverse_reads;   ///< Reverse strand reads
    int num_filtered_reads;  ///< Reads filtered out (debug mode)
    double elapsed_ms;
    double peak_memory_mb;
    bool success;
    std::string error_message;

    // Distance matrix statistics
    int num_valid_pairs;         ///< Number of valid distance pairs
    int num_invalid_pairs;       ///< Number of invalid pairs (insufficient overlap)
    double avg_common_coverage;  ///< Average common CpG coverage per pair
    double pairwise_mean_distance;    ///< Mean pairwise distance across valid pairs
    double pairwise_median_distance;  ///< Median pairwise distance across valid pairs

    // Significance analysis results (Cluster-First)
    bool significance_computed;  ///< Whether significance analysis was run
    double global_p_value;       ///< Global association p-value (FFH)
    double cramers_v;            ///< Effect size (Cramér's V)
    double local_best_p_value;   ///< Best cluster's local p-value
    int local_best_cluster;      ///< Best cluster ID
    double heuristic_score;      ///< Combined heuristic score [0-1]
    bool passed_gating;          ///< Whether passed gating (global_p <= 0.1)
    double cluster_permanova_f;      ///< Cluster-based PERMANOVA pseudo-F
    double cluster_permanova_p;      ///< Cluster-based PERMANOVA p-value
    bool cluster_permanova_valid;    ///< Whether cluster-based PERMANOVA is valid
    double cluster_dispersion_p;     ///< Cluster-based dispersion p-value
    bool cluster_dispersion_warning; ///< Whether cluster-based dispersion warns

    // Label-First verification results
    bool label_test_computed;    ///< Whether label-first test was run
    double label_delta;          ///< Delta = d_between - d_within (deprecated, use Stage 1)
    double label_p_value;        ///< Label permutation test p-value (deprecated, use Stage 1)
    bool label_significant;      ///< Whether label test is significant (p <= 0.05)
    std::string dominant_label;  ///< Which label dimension is most significant ("hp", "allele", "none")
    double label_hp_permanova_f;      ///< Label-based HP PERMANOVA pseudo-F
    double label_hp_permanova_p;      ///< Label-based HP PERMANOVA p-value
    bool label_hp_permanova_valid;    ///< Whether label-based HP PERMANOVA is valid
    double label_hp_dispersion_p;     ///< Label-based HP dispersion p-value
    bool label_hp_dispersion_warning; ///< Whether label-based HP dispersion warns
    double label_allele_permanova_f;      ///< Label-based allele PERMANOVA pseudo-F
    double label_allele_permanova_p;      ///< Label-based allele PERMANOVA p-value
    bool label_allele_permanova_valid;    ///< Whether label-based allele PERMANOVA is valid
    double label_allele_dispersion_p;     ///< Label-based allele dispersion p-value
    bool label_allele_dispersion_warning; ///< Whether label-based allele dispersion warns

    // Multi-Stage HP Verification results (NEW)
    // Stage 1: HP Family Merged Test
    double hp_merged_delta;      ///< Delta for (HP1+HP1-1) vs (HP2+HP2-1)
    double hp_merged_p;          ///< P-value for merged HP test
    bool hp_merged_sig;          ///< Significant at alpha
    int hp_merged_n_hp1;         ///< HP1-family count (HP1 + HP1-1)
    int hp_merged_n_hp2;         ///< HP2-family count (HP2 + HP2-1)

    // Stage 2: HP Fine-Grained Test
    double hp_fine_f;            ///< Pseudo-F statistic
    double hp_fine_p;            ///< P-value for fine-grained test
    bool hp_fine_sig;            ///< Significant at alpha
    int hp_fine_n_groups;        ///< Number of valid groups (max 4)

    // Stage 3: Allele Test
    double allele_delta;         ///< Delta for ALT vs REF
    double allele_p;             ///< P-value for allele test
    bool allele_sig;             ///< Significant at alpha

    // Stage 4: Unassigned Affinity Test
    double unassigned_affinity;  ///< Affinity score for HP3/HP0 reads
    double unassigned_affinity_p;///< P-value for affinity test
    std::string unassigned_dir;  ///< Affinity direction ("HP1", "HP2", "None")
    int unassigned_n_hp3;        ///< Number of HP3 reads
    int unassigned_n_hp0;        ///< Number of HP0/unphased reads

    // Cluster stability results (NEW)
    double cluster_stability;    ///< Stability score from subsampling [0-1]
    bool has_outlier_cluster;    ///< True if any cluster has < 3 reads after retry
    int n_clusters;              ///< Number of clusters found

    // Bidirectional verification classification (NEW)
    std::string verification_class;  ///< "Strong", "Subclone", "Weak", or "Noise"

    // Multi-Layer Validation Quality Metrics (NEW - Phase 5)
    double hp_ratio;                 ///< HP1/(HP1+HP2), range [0,1]
    bool potential_loh;              ///< True if HP ratio < 0.1 or > 0.9
    double coverage_multiple;        ///< NumReads/75.0 (expected coverage)
    std::string coverage_category;   ///< "Normal", "Low", "High", "CNV_Loss", "CNV_Gain", "High_Copy"
    std::string loh_subtype;         ///< "None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"
    float quality_score;             ///< Composite quality score [0-100]
    std::string quality_tier;        ///< "High" (>=70), "Medium" (40-69), "Low" (<40)

    RegionResult()
        : region_id(-1),
          snv_id(-1),
          num_reads(0),
          num_cpgs(0),
          num_forward_reads(0),
          num_reverse_reads(0),
          num_filtered_reads(0),
          elapsed_ms(0.0),
          peak_memory_mb(0.0),
          success(false),
          num_valid_pairs(0),
          num_invalid_pairs(0),
          avg_common_coverage(0.0),
          pairwise_mean_distance(0.0),
          pairwise_median_distance(0.0),
          significance_computed(false),
          global_p_value(1.0),
          cramers_v(0.0),
          local_best_p_value(1.0),
          local_best_cluster(-1),
          heuristic_score(0.0),
          passed_gating(false),
          cluster_permanova_f(0.0),
          cluster_permanova_p(1.0),
          cluster_permanova_valid(false),
          cluster_dispersion_p(1.0),
          cluster_dispersion_warning(false),
          label_test_computed(false),
          label_delta(0.0),
          label_p_value(1.0),
          label_significant(false),
          dominant_label("none"),
          label_hp_permanova_f(0.0),
          label_hp_permanova_p(1.0),
          label_hp_permanova_valid(false),
          label_hp_dispersion_p(1.0),
          label_hp_dispersion_warning(false),
          label_allele_permanova_f(0.0),
          label_allele_permanova_p(1.0),
          label_allele_permanova_valid(false),
          label_allele_dispersion_p(1.0),
          label_allele_dispersion_warning(false),
          hp_merged_delta(0.0),
          hp_merged_p(1.0),
          hp_merged_sig(false),
          hp_merged_n_hp1(0),
          hp_merged_n_hp2(0),
          hp_fine_f(0.0),
          hp_fine_p(1.0),
          hp_fine_sig(false),
          hp_fine_n_groups(0),
          allele_delta(0.0),
          allele_p(1.0),
          allele_sig(false),
          unassigned_affinity(0.0),
          unassigned_affinity_p(1.0),
          unassigned_dir("None"),
          unassigned_n_hp3(0),
          unassigned_n_hp0(0),
          cluster_stability(0.0),
          has_outlier_cluster(false),
          n_clusters(0),
          verification_class("Noise"),
          hp_ratio(0.5),
          potential_loh(false),
          coverage_multiple(1.0),
          coverage_category("Normal"),
          loh_subtype("None"),
          quality_score(50.0f),
          quality_tier("Medium") {
    }
};

/**
 * @brief Core class for parallel processing of multiple SNV regions
 *
 * This class is responsible for:
 * 1. Loading SNV table
 * 2. Defining region for each SNV (e.g., +/-2000bp)
 * 3. Parallel processing of multiple regions using OpenMP
 * 4. Managing thread-local resources (BamReader, FastaReader)
 * 5. Collecting and reporting processing results for each region
 * 6. Recording filtered reads in debug mode
 *
 * Thread-safety:
 * - Each thread maintains its own BAM/FASTA readers
 * - MatrixBuilder and RegionWriter are used within critical sections
 * - Result collection is mutex-protected
 */
class RegionProcessor {
public:
    /**
     * @brief Construct RegionProcessor (simplified version, for backward compatibility)
     */
    RegionProcessor(const std::string& tumor_bam_path, const std::string& normal_bam_path,
                    const std::string& ref_fasta_path, const std::string& output_dir, int num_threads = 4,
                    int32_t window_size = 2000);

    /**
     * @brief Construct RegionProcessor (full version, using Config)
     *
     * @param config Configuration object containing all parameters
     */
    explicit RegionProcessor(const Config& config);

    /**
     * @brief Load SNV table (TSV format)
     *
     * Format: chr  pos  ref  alt  qual
     * Example: chr17  7578000  C  T  100.0
     *
     * @param snv_table_path SNV table file path
     * @return Number of SNVs successfully loaded
     */
    int load_snvs(const std::string& snv_table_path);

    /**
     * @brief Load SNVs from VCF file
     *
     * @param vcf_path Path to VCF file
     * @return Number of SNVs loaded
     */
    int load_snvs_from_vcf(const std::string& vcf_path);

    /**
     * @brief Process all SNVs (parallelized)
     *
     * @param max_snvs Maximum number of SNVs to process (0 = all)
     * @return Vector of processing results
     */
    std::vector<RegionResult> process_all_regions(int max_snvs = 0);

    /**
     * @brief Process a single region (called by OpenMP worker thread)
     *
     * @param snv SNV information
     * @param region_id Region ID
     * @param bam_reader Thread-local BAM reader
     * @param fasta_reader Thread-local FASTA reader
     * @return RegionResult
     */
    RegionResult process_single_region(const SomaticSnv& snv, int region_id, BamReader& bam_reader,
                                       FastaReader& fasta_reader);

    /**
     * @brief Get loaded SNVs list
     */
    const std::vector<SomaticSnv>& get_snvs() const {
        return snvs_;
    }

    /**
     * @brief Output processing summary report
     */
    void print_summary(const std::vector<RegionResult>& results) const;

private:
    /**
     * @brief Write significance summary CSV and statistics report
     */
    void write_significance_summary(const std::vector<RegionResult>& results) const;

    // ========== Refactored Helper Methods ==========

    /**
     * @brief Process fetched reads and build methylation matrix
     *
     * @param reads Vector of BAM records (already fetched)
     * @param snv SNV information
     * @param ref_seq Reference sequence
     * @param region_start Region start coordinate (for coordinate mapping)
     * @param matrix_builder Matrix builder to populate
     * @param filtered_reads Vector to collect filtered reads (debug mode)
     * @param result RegionResult to update with strand counts
     */
    void process_reads(const std::vector<bam1_t*>& reads, const SomaticSnv& snv, const std::string& ref_seq,
                       int32_t region_start, MatrixBuilder& matrix_builder,
                       std::vector<FilteredReadInfo>& filtered_reads, RegionResult& result);

    /**
     * @brief Build MethylationMatrix from MatrixBuilder for distance calculation
     *
     * @param matrix_builder Source matrix builder
     * @param region_id Region identifier
     * @return Constructed MethylationMatrix
     */
    MethylationMatrix build_methylation_matrix(const MatrixBuilder& matrix_builder, int region_id);

    /**
     * @brief Compute distance matrices and optionally write to disk
     *
     * @param meth_mat MethylationMatrix input
     * @param read_list Read information list
     * @param region_dir Region output directory
     * @param metric Distance metric to use
     * @param result RegionResult to update with statistics
     * @return Computed distance matrix
     */
    DistanceMatrix compute_and_write_distance_matrix(const MethylationMatrix& meth_mat,
                                                      const std::vector<ReadInfo>& read_list,
                                                      const std::string& region_dir, DistanceMetricType metric,
                                                      RegionResult& result);

    /**
     * @brief Perform hierarchical clustering and significance analysis
     *
     * @param all_dist Distance matrix for all reads
     * @param read_list Read information list
     * @param meth_mat MethylationMatrix for analysis
     * @param clustering_dir Directory for clustering output
     * @param chr_name Chromosome name
     * @param snv SNV information
     * @param region_id Region identifier
     * @param result RegionResult to update with analysis results
     */
    void perform_clustering_and_significance(const DistanceMatrix& all_dist, const std::vector<ReadInfo>& read_list,
                                              const MethylationMatrix& meth_mat, const std::string& clustering_dir,
                                              const std::string& chr_name, const SomaticSnv& snv, int region_id,
                                              RegionResult& result);

    /**
     * @brief Write strand-specific clustering trees
     *
     * @param forward_dist Forward strand distance matrix
     * @param reverse_dist Reverse strand distance matrix
     * @param read_list Read information list
     * @param clustering_dir Directory for clustering output
     */
    void write_strand_specific_trees(const DistanceMatrix& forward_dist, const DistanceMatrix& reverse_dist,
                                      const std::vector<ReadInfo>& read_list, const std::string& clustering_dir);

    std::string tumor_bam_path_;
    std::string normal_bam_path_;
    std::string ref_fasta_path_;
    std::string output_dir_;
    std::string debug_output_dir_;
    std::string vcf_filename_;  ///< VCF filename (without path and extension)
    int num_threads_;
    int32_t window_size_;

    // Configuration
    LogLevel log_level_;
    bool output_filtered_reads_;
    bool no_filter_output_;
    ReadFilterConfig filter_config_;
    bool use_full_read_span_;  ///< If true, dynamically expand window to cover full span of reads

    // Distance matrix configuration
    bool compute_distance_matrix_;
    bool output_distance_matrix_;
    bool output_strand_distance_matrices_;
    DistanceConfig distance_config_;
    std::vector<DistanceMetricType> distance_metrics_;

    // Hierarchical clustering configuration
    bool compute_clustering_;
    bool output_tree_files_;
    bool output_linkage_matrix_;
    LinkageMethod linkage_method_;
    int clustering_min_reads_;

    std::vector<SomaticSnv> snvs_;
    ChromIndex chrom_index_;  // Manage chromosome name to ID mapping

    // Thread-local resources are created in process_single_region
};

}  // namespace InterSubMod
