#pragma once

#include <memory>
#include <string>
#include <vector>

#include "core/BamReader.hpp"
#include "core/Config.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/HierarchicalClustering.hpp"
#include "core/MatrixBuilder.hpp"
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

    // Significance analysis results (Cluster-First)
    bool significance_computed;  ///< Whether significance analysis was run
    double global_p_value;       ///< Global association p-value (FFH)
    double cramers_v;            ///< Effect size (Cramér's V)
    double local_best_p_value;   ///< Best cluster's local p-value
    int local_best_cluster;      ///< Best cluster ID
    double heuristic_score;      ///< Combined heuristic score [0-1]
    bool passed_gating;          ///< Whether passed gating (global_p <= 0.1)

    // Label-First verification results (NEW)
    bool label_test_computed;    ///< Whether label-first test was run
    double label_delta;          ///< Delta = d_between - d_within
    double label_p_value;        ///< Label permutation test p-value
    bool label_significant;      ///< Whether label test is significant (p <= 0.05)
    std::string dominant_label;  ///< Which label dimension is most significant ("hp", "allele", "none")

    // Cluster stability results (NEW)
    double cluster_stability;    ///< Stability score from subsampling [0-1]
    bool has_outlier_cluster;    ///< True if any cluster has < 3 reads after retry
    int n_clusters;              ///< Number of clusters found

    // Bidirectional verification classification (NEW)
    std::string verification_class;  ///< "Strong", "Subclone", "Weak", or "Noise"

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
          significance_computed(false),
          global_p_value(1.0),
          cramers_v(0.0),
          local_best_p_value(1.0),
          local_best_cluster(-1),
          heuristic_score(0.0),
          passed_gating(false),
          label_test_computed(false),
          label_delta(0.0),
          label_p_value(1.0),
          label_significant(false),
          dominant_label("none"),
          cluster_stability(0.0),
          has_outlier_cluster(false),
          n_clusters(0),
          verification_class("Noise") {
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
