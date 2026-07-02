#include "core/RegionProcessor.hpp"

#include <omp.h>
#include <sys/stat.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <random>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <tuple>
#include <unordered_set>

#include "core/HierarchicalClustering.hpp"
#include "core/PhyloLabeler.hpp"
#include "core/MethylationMatrix.hpp"
#include "core/NormalBaseline.hpp"
#include "core/SignificanceAnalyzer.hpp"
#include "core/StructureTest.hpp"
#include "io/TreeWriter.hpp"
#include "utils/Logger.hpp"

namespace InterSubMod {

// ============================================================================
// Quality Assessment Helper Functions (Multi-Layer Validation)
// ============================================================================

namespace {

/**
 * @brief Compute HP ratio from HP family counts
 */
double compute_hp_ratio(int hp1_family_n, int hp2_family_n) {
    // Add small epsilon to avoid division by zero
    double total = static_cast<double>(hp1_family_n + hp2_family_n) + 0.002;
    return (static_cast<double>(hp1_family_n) + 0.001) / total;
}

/**
 * @brief Determine if HP ratio indicates potential LOH
 */
bool is_potential_loh(double hp_ratio) {
    return (hp_ratio < 0.1) || (hp_ratio > 0.9);
}

}  // anonymous namespace

/**
 * @brief Compute coverage multiple relative to expected diploid coverage.
 *
 * The expected_coverage defaults to 75.0 for initial per-region computation.
 * After all regions are processed, this is recalculated with the sample-specific
 * diploid coverage estimated via KDE mode of the NumReads distribution.
 *
 * Exposed at namespace scope (not anonymous) to allow unit testing.
 */
double compute_coverage_multiple(int num_reads, double expected_coverage) {
    return static_cast<double>(num_reads) / expected_coverage;
}

/**
 * @brief Estimate diploid coverage from the distribution of per-region read counts.
 *
 * Uses a histogram-based mode estimation (equivalent to KDE mode with fixed bandwidth).
 * The mode of the NumReads distribution corresponds to diploid (CN=2) regions,
 * which are typically the most abundant genomic state.
 *
 * Validated against SEQC2 ground truth: estimate error = 7.5% (HCC1395).
 *
 * Falls back to 75.0 if insufficient data (< 100 regions).
 */
DiploidEstimate estimate_diploid_coverage(const std::vector<RegionResult>& results) {
    // Collect valid num_reads
    std::vector<double> nr_values;
    nr_values.reserve(results.size());
    for (const auto& r : results) {
        if (r.success && r.num_reads > 5) {
            nr_values.push_back(static_cast<double>(r.num_reads));
        }
    }

    if (nr_values.size() < 100) {
        LOG_WARN("KDE diploid estimate falling back to 75.0x: insufficient valid regions (" +
                 std::to_string(nr_values.size()) + " < 100). CovM values will NOT be per-sample calibrated.");
        return {75.0, true};  // Fallback: insufficient data
    }

    // Find 80th percentile as upper search bound (exclude extreme gains)
    std::vector<double> sorted_nr = nr_values;
    std::sort(sorted_nr.begin(), sorted_nr.end());
    double p80 = sorted_nr[static_cast<size_t>(sorted_nr.size() * 0.8)];
    double lower = 10.0;
    double upper = std::max(p80, 50.0);

    // Histogram-based mode with Gaussian smoothing (bin_size=2, sigma=3 bins)
    int bin_size = 2;
    int n_bins = static_cast<int>((upper - lower) / bin_size) + 1;
    if (n_bins < 10) {
        LOG_WARN("KDE diploid estimate falling back to 75.0x: histogram range too narrow (n_bins=" +
                 std::to_string(n_bins) + " < 10, lower=" + std::to_string(static_cast<int>(lower)) +
                 ", upper=" + std::to_string(static_cast<int>(upper)) + ").");
        return {75.0, true};
    }

    std::vector<double> hist(n_bins, 0.0);
    for (double v : nr_values) {
        if (v >= lower && v < upper) {
            int idx = static_cast<int>((v - lower) / bin_size);
            if (idx >= 0 && idx < n_bins) {
                hist[idx] += 1.0;
            }
        }
    }

    // Gaussian smoothing (sigma = 3 bins, window = 2*3*sigma+1 = 19)
    int sigma_bins = 3;
    int half_win = 3 * sigma_bins;
    std::vector<double> kernel(2 * half_win + 1);
    double kernel_sum = 0.0;
    for (int i = -half_win; i <= half_win; i++) {
        kernel[i + half_win] = std::exp(-0.5 * (i * i) / (sigma_bins * sigma_bins));
        kernel_sum += kernel[i + half_win];
    }
    for (auto& k : kernel) k /= kernel_sum;

    std::vector<double> smoothed(n_bins, 0.0);
    for (int i = 0; i < n_bins; i++) {
        for (int j = -half_win; j <= half_win; j++) {
            int src = i + j;
            if (src >= 0 && src < n_bins) {
                smoothed[i] += hist[src] * kernel[j + half_win];
            }
        }
    }

    // Find the peak (mode)
    int peak_idx = 0;
    double peak_val = 0.0;
    for (int i = 0; i < n_bins; i++) {
        if (smoothed[i] > peak_val) {
            peak_val = smoothed[i];
            peak_idx = i;
        }
    }

    double diploid_est = lower + (peak_idx + 0.5) * bin_size;

    // Sanity check: diploid estimate should be reasonable (10-300)
    if (diploid_est < 10.0 || diploid_est > 300.0) {
        LOG_WARN("KDE diploid estimate falling back to 75.0x: mode estimate " +
                 std::to_string(static_cast<int>(diploid_est)) + "x out of sanity range [10, 300].");
        return {75.0, true};  // Fallback
    }

    return {diploid_est, false};
}

namespace {

/**
 * @brief Determine coverage category based on coverage multiple
 */
std::string determine_coverage_category(double coverage_multiple) {
    if (coverage_multiple < 0.5) {
        return "CNV_Loss";
    } else if (coverage_multiple < 0.8) {
        return "Low";
    } else if (coverage_multiple <= 1.2) {
        return "Normal";
    } else if (coverage_multiple <= 1.5) {
        return "Elevated";
    } else if (coverage_multiple <= 2.0) {
        return "CNV_Gain";
    } else {
        return "High_Copy";
    }
}

/**
 * @brief Determine LOH subtype based on verification class and LOH status
 */
std::string determine_loh_subtype(bool potential_loh, const std::string& verification_class) {
    if (!potential_loh) {
        return "None";
    }

    if (verification_class == "Noise") {
        return "LOH_Noise";
    } else if (verification_class == "Weak") {
        return "LOH_Weak";
    } else if (verification_class == "Strong") {
        return "LOH_Strong";
    } else if (verification_class == "Subclone") {
        return "LOH_Subclone";
    }

    return "None";
}

/**
 * @brief Mode-aware quality score weights.
 *
 * Paired and tumor-only modes use different weights because TO phasing artifacts
 * inflate LOH-like rates and the verification bonus rewards FP > TP in TO mode
 * (validated by AUC analysis 2026-03-31).
 */
struct QualityScoreWeights {
    float read_penalty_severe = 20.0f;   // num_reads < 30
    float read_penalty_moderate = 10.0f; // num_reads < 50
    float cpg_penalty_severe = 15.0f;    // num_cpgs < 20
    float cpg_penalty_moderate = 10.0f;  // num_cpgs < 30
    float cov_penalty_cnv_loss = 30.0f;  // coverage_multiple < 0.5
    float cov_penalty_low = 15.0f;       // coverage_multiple < 0.8
    float cov_penalty_high_copy = 20.0f; // coverage_multiple > 2.0
    float loh_penalty = 25.0f;           // potential LOH
    float verify_bonus = 10.0f;          // HP and Allele tests both significant
    float effect_bonus_strong = 15.0f;   // cramers_v >= 0.3
    float effect_bonus_moderate = 5.0f;  // cramers_v >= 0.2
};

QualityScoreWeights get_paired_weights() {
    return QualityScoreWeights{};  // all defaults
}

QualityScoreWeights get_tumor_only_weights() {
    QualityScoreWeights w{};
    w.loh_penalty = 0.0f;    // TO LOH-like is dominated by phasing artifacts
    w.verify_bonus = 0.0f;   // TO verify bonus rewards FP more than TP
    return w;
}

/**
 * @brief Compute comprehensive quality score using mode-aware weights.
 */
float compute_quality_score(int num_reads, int num_cpgs, double coverage_multiple,
                            bool potential_loh, bool hp_sig, bool allele_sig, double cramers_v,
                            const QualityScoreWeights& weights) {
    float score = 100.0f;

    // Read count penalty
    if (num_reads < 30) {
        score -= weights.read_penalty_severe;
    } else if (num_reads < 50) {
        score -= weights.read_penalty_moderate;
    }

    // CpG count penalty
    if (num_cpgs < 20) {
        score -= weights.cpg_penalty_severe;
    } else if (num_cpgs < 30) {
        score -= weights.cpg_penalty_moderate;
    }

    // Coverage anomaly penalty
    if (coverage_multiple < 0.5) {
        score -= weights.cov_penalty_cnv_loss;
    } else if (coverage_multiple < 0.8) {
        score -= weights.cov_penalty_low;
    } else if (coverage_multiple > 2.0) {
        score -= weights.cov_penalty_high_copy;
    }

    // LOH penalty
    if (potential_loh) {
        score -= weights.loh_penalty;
    }

    // Verification consistency bonus
    if (hp_sig && allele_sig) {
        score += weights.verify_bonus;
    }

    // Effect size bonus
    if (cramers_v >= 0.3) {
        score += weights.effect_bonus_strong;
    } else if (cramers_v >= 0.2) {
        score += weights.effect_bonus_moderate;
    }

    // Clamp to valid range
    if (score < 0.0f) score = 0.0f;
    if (score > 100.0f) score = 100.0f;

    return score;
}

/**
 * @brief Determine quality tier from quality score
 */
std::string determine_quality_tier(float quality_score) {
    if (quality_score >= 70.0f) {
        return "High";
    } else if (quality_score >= 40.0f) {
        return "Medium";
    } else {
        return "Low";
    }
}

/**
 * @brief Summarize valid off-diagonal distances from a distance matrix.
 */
std::pair<double, double> summarize_pairwise_distances(const Eigen::MatrixXd& dist_matrix) {
    std::vector<double> valid_distances;
    valid_distances.reserve(static_cast<size_t>(dist_matrix.rows()) * static_cast<size_t>(dist_matrix.rows()) / 2);

    for (int i = 0; i < dist_matrix.rows(); ++i) {
        for (int j = i + 1; j < dist_matrix.cols(); ++j) {
            double dist = dist_matrix(i, j);
            if (!std::isnan(dist)) {
                valid_distances.push_back(dist);
            }
        }
    }

    if (valid_distances.empty()) {
        return {0.0, 0.0};
    }

    double mean = std::accumulate(valid_distances.begin(), valid_distances.end(), 0.0) /
                  static_cast<double>(valid_distances.size());

    std::sort(valid_distances.begin(), valid_distances.end());
    size_t mid = valid_distances.size() / 2;
    double median = 0.0;
    if (valid_distances.size() % 2 == 0) {
        median = (valid_distances[mid - 1] + valid_distances[mid]) / 2.0;
    } else {
        median = valid_distances[mid];
    }

    return {mean, median};
}

}  // anonymous namespace

RegionProcessor::RegionProcessor(const std::string& tumor_bam_path, const std::string& normal_bam_path,
                                 const std::string& ref_fasta_path, const std::string& output_dir, int num_threads,
                                 int32_t window_size)
    : tumor_bam_path_(tumor_bam_path),
      normal_bam_path_(normal_bam_path),
      ref_fasta_path_(ref_fasta_path),
      output_dir_(output_dir),
      debug_output_dir_(output_dir + "/debug"),
      num_threads_(num_threads),
      window_size_(window_size),
      log_level_(LogLevel::LOG_INFO),
      output_filtered_reads_(false),
      no_filter_output_(false),
      compute_clustering_(true),
      output_tree_files_(true),
      output_linkage_matrix_(true),
      linkage_method_(LinkageMethod::UPGMA),
      clustering_min_reads_(10) {
    // Set OpenMP threads
    omp_set_num_threads(num_threads_);

    use_full_read_span_ = false;
    expected_coverage_ = 0.0;

    LOG_INFO("RegionProcessor initialized with " + std::to_string(num_threads_) + " threads, window_size=±" +
             std::to_string(window_size_) + "bp");
}

RegionProcessor::RegionProcessor(const Config& config)
    : tumor_bam_path_(config.tumor_bam_path),
      normal_bam_path_(config.normal_bam_path),
      ref_fasta_path_(config.reference_fasta_path),
      output_dir_(config.output_dir),
      debug_output_dir_(config.get_debug_output_dir()),
      num_threads_(config.threads),
      window_size_(config.window_size_bp),
      log_level_(config.log_level),
      output_filtered_reads_(config.output_filtered_reads),
      no_filter_output_(config.no_filter_output),
      compute_distance_matrix_(config.compute_distance_matrix),
      output_distance_matrix_(config.output_distance_matrix),
      output_strand_distance_matrices_(config.output_strand_distance_matrices),
      distance_metrics_(config.distance_metrics),
      compute_clustering_(config.compute_clustering),
      output_tree_files_(config.output_tree_files),
      output_linkage_matrix_(config.output_linkage_matrix),
      clustering_min_reads_(config.clustering_min_reads),
      binary_methyl_high_(config.binary_methyl_high),
      binary_methyl_low_(config.binary_methyl_low) {
    // Extract VCF filename from config (remove path and extension)
    std::filesystem::path vcf_path(config.somatic_vcf_path);
    std::string filename = vcf_path.stem().string();
    // If extension is .vcf.gz, remove .vcf as well
    if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".vcf") {
        filename = filename.substr(0, filename.size() - 4);
    }
    vcf_filename_ = filename;

    // Parse linkage method from string
    linkage_method_ = HierarchicalClustering::string_to_method(config.linkage_method);

    // Set full read span mode
    use_full_read_span_ = config.use_full_read_span;

    // Coverage normalization
    expected_coverage_ = config.expected_coverage;

    // Set filter configuration
    filter_config_.min_mapq = config.min_mapq;
    filter_config_.min_read_length = config.min_read_length;
    filter_config_.min_base_quality = config.min_base_quality;
    filter_config_.require_mm_ml = true;
    filter_config_.germline_hp_only = config.germline_hp_only;

    // Set distance matrix configuration
    if (!distance_metrics_.empty()) {
        distance_config_.metric = distance_metrics_[0];  // Default to first
    } else {
        distance_config_.metric = DistanceMetricType::NHD;
        distance_metrics_.push_back(DistanceMetricType::NHD);
    }
    distance_config_.min_common_coverage = config.min_common_coverage;
    distance_config_.nan_strategy = config.nan_distance_strategy;
    distance_config_.max_distance_value = config.max_distance_value;
    distance_config_.use_binary_matrix = config.distance_use_binary;
    distance_config_.pearson_center = config.distance_pearson_center;
    distance_config_.jaccard_include_unmeth = config.distance_jaccard_include_unmeth;
    distance_config_.num_threads = 1;  // Single thread for per-region computation

    // Set OpenMP threads
    omp_set_num_threads(num_threads_);

    // Create debug directory if needed
    if (output_filtered_reads_) {
        mkdir(debug_output_dir_.c_str(), 0755);
    }

    std::stringstream ss;
    ss << "RegionProcessor initialized:\n"
       << "  Threads: " << num_threads_ << "\n"
       << "  Window size: ±" << window_size_ << " bp\n"
       << "  Log level: " << static_cast<int>(log_level_) << "\n"
       << "  Distance metrics: ";
    for (size_t i = 0; i < distance_metrics_.size(); ++i) {
        ss << (i > 0 ? ", " : "") << DistanceCalculator::metric_to_string(distance_metrics_[i]);
    }
    ss << "\n";
    ss << "  Min common coverage (C_min): " << distance_config_.min_common_coverage << "\n";
    if (output_filtered_reads_) {
        ss << "  Debug output: " << debug_output_dir_ << "\n";
    }
    if (no_filter_output_) {
        ss << "  Mode: No-filter (outputting all reads)\n";
    }
    if (compute_distance_matrix_) {
        ss << "  Distance matrix: enabled\n";
        if (output_strand_distance_matrices_) {
            ss << "  Strand-specific matrices: enabled\n";
        }
    }
    if (compute_clustering_) {
        ss << "  Hierarchical clustering: enabled\n"
           << "  Linkage method: " << HierarchicalClustering::method_to_string(linkage_method_) << "\n"
           << "  Clustering min reads: " << clustering_min_reads_;
    }

    LOG_INFO(ss.str());

    // Load LOH BED file if provided (Phase C)
    if (!config.loh_bed_path.empty()) {
        int n_regions = loh_annotator_.load(config.loh_bed_path);
        if (n_regions < 0) {
            LOG_WARN("Failed to load LOH BED file: " + config.loh_bed_path);
        } else {
            LOG_INFO("Loaded " + std::to_string(n_regions) + " LOH regions from BED file");
        }
    }
}

int RegionProcessor::load_snvs(const std::string& snv_table_path) {
    std::ifstream ifs(snv_table_path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open SNV table: " << snv_table_path << std::endl;
        return 0;
    }

    snvs_.clear();
    std::string line;
    int line_num = 0;

    // Skip header if present
    std::getline(ifs, line);
    line_num++;

    // Check if first line is header
    bool has_header = (line.find("chr") != std::string::npos || line.find("pos") != std::string::npos);

    if (!has_header) {
        // First line is data, parse it
        std::istringstream iss(line);
        std::string chr_str;
        uint32_t pos;
        char ref, alt;
        float qual = 0.0f;

        if (iss >> chr_str >> pos >> ref >> alt) {
            iss >> qual;  // Optional

            SomaticSnv snv;
            snv.snv_id = 0;
            snv.chr_id = chrom_index_.get_or_create_id(chr_str);
            snv.pos = pos;
            snv.ref_base = ref;
            snv.alt_base = alt;
            snv.qual = qual;

            snvs_.push_back(snv);
        }
    }

    // Parse remaining lines
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string chr_str;
        uint32_t pos;
        char ref, alt;
        float qual = 0.0f;

        if (iss >> chr_str >> pos >> ref >> alt) {
            iss >> qual;  // Optional

            SomaticSnv snv;
            snv.snv_id = snvs_.size();
            snv.chr_id = chrom_index_.get_or_create_id(chr_str);
            snv.pos = pos;
            snv.ref_base = ref;
            snv.alt_base = alt;
            snv.qual = qual;

            snvs_.push_back(snv);
        } else {
            std::cerr << "Failed to parse SNV at line " << line_num << ": " << line << std::endl;
        }
    }

    ifs.close();
    LOG_INFO("Loaded " + std::to_string(snvs_.size()) + " SNVs from " + snv_table_path);
    return snvs_.size();
}

int RegionProcessor::load_snvs_from_vcf(const std::string& vcf_path) {
    SomaticSnvTable table;
    if (table.load_from_vcf(vcf_path, chrom_index_)) {
        snvs_ = table.all();
        LOG_INFO("Loaded " + std::to_string(snvs_.size()) + " SNVs from VCF: " + vcf_path);
        return snvs_.size();
    }

    LOG_ERROR("Failed to load SNVs from VCF: " + vcf_path);
    return 0;
}

std::vector<RegionResult> RegionProcessor::process_all_regions(int max_snvs) {
    int num_to_process = (max_snvs > 0 && max_snvs < static_cast<int>(snvs_.size())) ? max_snvs : snvs_.size();

    std::cout << "Processing " << num_to_process << " regions with " << num_threads_ << " threads..." << std::endl;

    std::vector<RegionResult> results(num_to_process);

    auto t_start = std::chrono::high_resolution_clock::now();

    LOG_INFO("Processing " + std::to_string(num_to_process) + " regions with " + 
             std::to_string(num_threads_) + " threads...");

    // Progress tracking variables (atomic for thread safety)
    std::atomic<int> processed_count{0};
    std::string current_chr = "";
    std::mutex chr_mutex;
    int last_progress_pct = -1;

// OpenMP parallel loop
// OpenMP parallel region to manage thread-local resources
#pragma omp parallel
    {
        // Initialize thread-local readers once per thread
        // This avoids opening/closing the BAM file and loading the index for every region
        BamReader tumor_reader(tumor_bam_path_);
        FastaReader ref_reader(ref_fasta_path_);

        // Optional: Normal BAM reader if path is provided
        std::unique_ptr<BamReader> normal_reader;
        if (!normal_bam_path_.empty()) {
            normal_reader = std::make_unique<BamReader>(normal_bam_path_);
        }

#pragma omp for schedule(dynamic)
        for (int i = 0; i < num_to_process; i++) {
            const auto& snv = snvs_[i];
            std::string chr_name = chrom_index_.get_name(snv.chr_id);

            // Process the region using the thread-local readers
            results[i] = process_single_region(snv, i, tumor_reader, ref_reader,
                                                    normal_reader.get());

            // Update progress counter
            int count = ++processed_count;

            // Calculate progress percentage
            int progress_pct = (count * 100) / num_to_process;

            // Log per-region completion at DEBUG level
            if (results[i].success) {
                std::stringstream ss;
                ss << "Region " << i << " (" << chr_name << ":" << snv.pos << ") completed: " << results[i].num_reads
                   << " reads, " << std::fixed << std::setprecision(1) << results[i].elapsed_ms << " ms";
                LOG_DEBUG(ss.str());
            } else {
                std::stringstream ss;
                ss << "Region " << i << " (" << chr_name << ":" << snv.pos << ") failed: " << results[i].error_message;
                LOG_ERROR(ss.str());
            }

            // Progress reporting at INFO level (every 5% or chromosome change)
            bool chr_changed = false;
            {
                std::lock_guard<std::mutex> lock(chr_mutex);
                if (current_chr != chr_name) {
                    current_chr = chr_name;
                    chr_changed = true;
                }
            }

            // Report progress at 5% intervals or when chromosome changes
            if (chr_changed || (progress_pct >= last_progress_pct + 5 && progress_pct != last_progress_pct)) {
                last_progress_pct = progress_pct;

                // Calculate ETA
                auto now = std::chrono::high_resolution_clock::now();
                double wall_elapsed_ms = std::chrono::duration<double, std::milli>(now - t_start).count();
                double avg_ms_per_region = wall_elapsed_ms / count;
                int remaining = num_to_process - count;
                double eta_seconds = (avg_ms_per_region * remaining) / 1000.0;

                std::stringstream ss;
                ss << "Progress: " << count << "/" << num_to_process 
                   << " (" << progress_pct << "%) | " << chr_name
                   << " | Avg: " << std::fixed << std::setprecision(0) << avg_ms_per_region << " ms/region"
                   << " | ETA: ";
                
                if (eta_seconds < 60) {
                    ss << std::fixed << std::setprecision(0) << eta_seconds << "s";
                } else if (eta_seconds < 3600) {
                    ss << std::fixed << std::setprecision(1) << (eta_seconds / 60.0) << " min";
                } else {
                    ss << std::fixed << std::setprecision(1) << (eta_seconds / 3600.0) << " hr";
                }
                
                LOG_INFO(ss.str());
            }
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_elapsed = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // Final summary
    std::stringstream final_ss;
    final_ss << "Completed " << num_to_process << " regions in " 
             << std::fixed << std::setprecision(1) << (total_elapsed / 1000.0) << "s"
             << " (avg: " << std::setprecision(1) << (total_elapsed / num_to_process) << " ms/region)";
    LOG_INFO(final_ss.str());

    // Determine diploid coverage: user-specified or auto-estimated
    double diploid_coverage;
    std::string diploid_source;
    if (expected_coverage_ > 0.0) {
        diploid_coverage = expected_coverage_;
        diploid_source = "user_specified";
        LOG_INFO("[CovM] Pass 2 diploid_coverage=" + std::to_string(diploid_coverage) +
                 "x source=user_specified (--expected-coverage)");
    } else {
        DiploidEstimate est = estimate_diploid_coverage(results);
        diploid_coverage = est.value;
        diploid_source = est.used_fallback ? "auto_kde_fallback_default" : "auto_kde";
        LOG_INFO("[CovM] Pass 2 diploid_coverage=" + std::to_string(diploid_coverage) +
                 "x source=" + diploid_source + " (KDE mode of NumReads distribution)");
    }

    auto qs_weights = normal_bam_path_.empty() ? get_tumor_only_weights() : get_paired_weights();
    for (auto& r : results) {
        if (r.success) {
            r.coverage_multiple = compute_coverage_multiple(r.num_reads, diploid_coverage);
            r.diploid_coverage_used = diploid_coverage;
            r.coverage_category = determine_coverage_category(r.coverage_multiple);
            r.quality_score = compute_quality_score(r.num_reads, r.num_cpgs, r.coverage_multiple,
                                                    r.potential_loh, r.hp_merged_sig, r.allele_sig,
                                                    r.cramers_v, qs_weights);
            r.quality_tier = determine_quality_tier(r.quality_score);
        }
    }

    // Phase D: Cross-region subclone analysis (when normal BAM is available)
    if (!normal_bam_path_.empty()) {
        SubcloneAnalyzer subclone_analyzer;
        SubcloneResult subclone_result = subclone_analyzer.analyze(results);
        if (subclone_result.valid) {
            // Write back subclone assignments to results
            // Assignments map 1:1 with valid (success + min reads/cpgs + significance) results
            int profile_idx = 0;
            for (auto& r : results) {
                if (!r.success || r.num_reads < 10 || r.num_cpgs < 3 || !r.significance_computed) continue;
                if (profile_idx < static_cast<int>(subclone_result.region_assignments.size())) {
                    r.subclone_id = subclone_result.region_assignments[profile_idx];
                }
                ++profile_idx;
            }

            // Write reports
            SubcloneAnalyzer::write_report(subclone_result, output_dir_ + "/subclone_structure.txt");
            LOG_INFO("Subclone analysis complete: " + std::to_string(subclone_result.n_subclones)
                     + " groups from " + std::to_string(subclone_result.total_regions_analyzed) + " regions");
        } else {
            LOG_INFO("Subclone analysis skipped: insufficient valid regions ("
                     + std::to_string(subclone_result.total_regions_analyzed) + ")");
        }
    }

    // Write significance summary and print stats
    write_significance_summary(results);
    print_summary(results);

    return results;
}

RegionResult RegionProcessor::process_single_region(const SomaticSnv& snv, int region_id, BamReader& bam_reader,
                                                    FastaReader& fasta_reader, BamReader* normal_reader) {
    RegionResult result;
    result.region_id = region_id;
    result.snv_id = snv.snv_id;

    auto t_start = std::chrono::high_resolution_clock::now();

    try {
        // Compute initial region window using RegionBounds helper
        std::string chr_name = chrom_index_.get_name(snv.chr_id);
        int32_t chr_length = fasta_reader.get_chr_length(chr_name);

        RegionBounds bounds = compute_region_bounds(chr_name,
                                                    static_cast<int32_t>(snv.pos),
                                                    static_cast<int32_t>(window_size_),
                                                    chr_length);
        int32_t region_start = bounds.start;
        int32_t region_end   = bounds.end;

        // Fetch reads from BAM.
        // #3 (2026-06-13) tumor SNV-point fetch: only reads crossing the somatic SNV are
        // useful for tumor (ALT/REF haplotyping); reads that merely overlap the ±window but
        // do not cross the SNV are filtered out downstream anyway (SNV_NOT_COVERED, ~40% at
        // ±5000). Fetching at the SNV point avoids fetching + alt-support-checking them.
        // Methylation is still read over the full ref_seq window (region_start..region_end)
        // since covering reads are long (ONT) and span the window. The alt_support filter is
        // kept as a safety net (e.g. SNV at a deletion within a covering read).
        // NOTE: normal reads (below) MUST keep the window fetch — they contribute the
        // methylation baseline even when not covering the SNV.
        int32_t snv_pos = static_cast<int32_t>(snv.pos);
        auto reads = bam_reader.fetch_reads(chr_name, snv_pos, snv_pos + 1);

        // Expand window to cover full span of all reads (if enabled)
        if (use_full_read_span_ && !reads.empty()) {
            int32_t min_start = region_start;
            int32_t max_end = region_end;

            for (auto* b : reads) {
                int32_t r_start = b->core.pos;
                int32_t r_end = bam_endpos(b);
                if (r_start < min_start) min_start = r_start;
                if (r_end > max_end) max_end = r_end;
            }

            region_start = std::max(1, min_start - 10);
            region_end = (chr_length > 0) ? std::min(chr_length, max_end + 10) : (max_end + 10);
        }

        // Fetch reference sequence
        std::string ref_seq = fasta_reader.fetch_sequence(chr_name, region_start, region_end);
        if (ref_seq.empty()) {
            throw std::runtime_error("Failed to fetch reference sequence");
        }

        // Aggregate reads using ReadAggregator (filter → parse → matrix)
        ReadAggregateConfig agg_config;
        agg_config.filter_config         = filter_config_;
        agg_config.no_filter_output      = no_filter_output_;
        agg_config.collect_filtered_reads = output_filtered_reads_;

        auto agg = ReadAggregator(agg_config).aggregate(reads, snv, ref_seq, region_start,
                                                          /*is_tumor=*/true);

        result.num_forward_reads = agg.num_forward_reads;
        result.num_reverse_reads = agg.num_reverse_reads;

        // Fetch and merge normal reads (if normal BAM is provided)
        // Normal reads: use same region bounds, skip alt-support filter. They KEEP their real
        // germline HP tag from the normal BAM (NOT forced to "0") — the germline-baseline
        // subtraction (germline ASM Δβ + tumor-vs-normal HP residual) depends on this; empirically
        // Normal_HP_Valid=true in ~96% of regions. Do NOT "fix" normal HP to "0".
        std::vector<bam1_t*> normal_reads;
        if (normal_reader != nullptr) {
            normal_reads = normal_reader->fetch_reads(chr_name, region_start, region_end);
            if (!normal_reads.empty()) {
                auto normal_agg = ReadAggregator(agg_config).aggregate(
                    normal_reads, snv, ref_seq, region_start, /*is_tumor=*/false);
                result.num_normal_reads = normal_agg.matrix_builder.num_reads();
                // Merge normal reads into tumor matrix before finalization
                agg.matrix_builder.merge(normal_agg.matrix_builder);
            }
        }
        result.num_tumor_reads = agg.matrix_builder.num_reads() - result.num_normal_reads;

        // Finalize matrix construction (now contains both tumor + normal reads)
        agg.matrix_builder.finalize();
        result.num_reads          = agg.matrix_builder.num_reads();
        result.num_cpgs           = agg.matrix_builder.num_cpgs();
        result.num_filtered_reads = agg.filtered_reads.size();

        // Write output to disk
        RegionWriter writer(output_dir_, debug_output_dir_, true, vcf_filename_);
        writer.write_region(snv, chr_name, region_id, region_start, region_end,
                            agg.matrix_builder.get_reads(),
                            agg.matrix_builder.get_cpg_positions(),
                            agg.matrix_builder.get_matrix(), 0.0, 0.0);

        // Write filtered reads in debug mode
        if (output_filtered_reads_ && !agg.filtered_reads.empty()) {
            std::string region_dir = writer.get_region_dir(chr_name, snv.pos, region_start, region_end);
            writer.write_filtered_reads(region_dir, chr_name, agg.filtered_reads);
        }

        // Distance matrix, clustering, and significance analysis
        if (compute_distance_matrix_ && result.num_reads >= 2 && result.num_cpgs >= 1) {
            MethylationMatrix meth_mat = build_methylation_matrix(agg.matrix_builder, region_id);
            const auto& read_list = agg.matrix_builder.get_reads();
            std::string region_dir = writer.get_region_dir(chr_name, snv.pos, region_start, region_end);

            // Audit: count somatic HP tags (HP:i:11/21/33) from the raw pre-demotion value.
            // Always runs regardless of the germline_hp_only flag so users can diagnose
            // self-phasing prevalence even when the downstream HP features already treat
            // these reads as unphased.
            for (const auto& r : read_list) {
                if (r.hp_tag_raw == "1-1") {
                    ++result.n_hp_somatic_11;
                } else if (r.hp_tag_raw == "2-1") {
                    ++result.n_hp_somatic_21;
                } else if (r.hp_tag_raw == "3") {
                    ++result.n_hp_somatic_33;
                }
            }

            // Compute distance matrix (using first metric for main analysis)
            DistanceMatrix all_dist =
                compute_and_write_distance_matrix(meth_mat, read_list, region_dir, distance_metrics_[0], result);

            // Clustering and significance analysis
            if (compute_clustering_ && result.num_reads >= clustering_min_reads_) {
                std::string clustering_dir = region_dir + "/clustering";
                perform_clustering_and_significance(all_dist, read_list, meth_mat, clustering_dir, chr_name, snv,
                                                    region_id, result);

                // Strand-specific trees
                if (output_strand_distance_matrices_ && output_tree_files_) {
                    DistanceConfig config = distance_config_;
                    config.metric = distance_metrics_[0];
                    DistanceCalculator dist_calc(config);
                    auto [forward_dist, reverse_dist] = dist_calc.compute_strand_specific(meth_mat, read_list);
                    write_strand_specific_trees(forward_dist, reverse_dist, read_list, clustering_dir);
                }
            }

            // Additional metrics (write distance matrices only, no clustering)
            for (size_t i = 1; i < distance_metrics_.size(); ++i) {
                RegionResult dummy_result;
                compute_and_write_distance_matrix(meth_mat, read_list, region_dir, distance_metrics_[i], dummy_result);
            }

            // Normal baseline + Somatic HP ASM (Phase B)
            // Compute normal baseline and HP ASM separately for tumor vs normal.
            // Key insight: subtracting per-column constant from distance is identity
            // (d(x-μ, y-μ) = d(x,y)), so residual matrix approach doesn't work.
            // Instead: compute HP delta on tumor-only and normal-only subsets,
            // then somatic HP ASM = tumor HP delta - normal HP delta.
            if (result.num_normal_reads >= 5) {
                // Build is_tumor mask from read list
                std::vector<bool> is_tumor_mask;
                is_tumor_mask.reserve(read_list.size());
                for (const auto& r : read_list) {
                    is_tumor_mask.push_back(r.is_tumor);
                }

                NormalBaseline baseline = build_normal_baseline(meth_mat.raw_matrix, is_tumor_mask);
                if (baseline.valid) {
                    result.normal_baseline_mean = baseline.overall_mean;
                    result.normal_baseline_coverage = baseline.mean_coverage;
                }

                // Build HP merged labels and separate tumor/normal index sets
                std::vector<int> tumor_indices, normal_indices;
                std::vector<int> hp_merged_labels(read_list.size(), -1);
                for (size_t i = 0; i < read_list.size(); ++i) {
                    const std::string& hp = read_list[i].hp_tag;
                    if (hp == "1" || hp == "HP1" || hp == "1-1") {
                        hp_merged_labels[i] = 0;
                    } else if (hp == "2" || hp == "HP2" || hp == "2-1") {
                        hp_merged_labels[i] = 1;
                    }
                    if (read_list[i].is_tumor) {
                        tumor_indices.push_back(static_cast<int>(i));
                    } else {
                        normal_indices.push_back(static_cast<int>(i));
                    }
                }

                // Extract tumor-only and normal-only distance submatrices
                // Then run HP merged test on each separately
                auto extract_and_test_hp = [&](const std::vector<int>& indices) -> LabelDeltaResult {
                    int n = static_cast<int>(indices.size());
                    if (n < 6) {
                        LabelDeltaResult r;
                        r.valid = false;
                        r.invalid_reason = "insufficient_reads";
                        return r;
                    }
                    Eigen::MatrixXd sub_dist(n, n);
                    std::vector<int> sub_labels(n, -1);
                    for (int i = 0; i < n; ++i) {
                        sub_labels[i] = hp_merged_labels[indices[i]];
                        for (int j = 0; j < n; ++j) {
                            sub_dist(i, j) = all_dist.dist_matrix(indices[i], indices[j]);
                        }
                    }
                    LabelTestConfig lt_config;
                    lt_config.seed = 42;
                    LabelTest lt(lt_config);
                    return lt.test_binary_groups(sub_dist, sub_labels);
                };

                // Count HP reads per sample for diagnostics
                for (size_t i = 0; i < read_list.size(); ++i) {
                    int hp_label = hp_merged_labels[i];
                    if (read_list[i].is_tumor) {
                        if (hp_label == 0) ++result.tumor_hp1_count;
                        else if (hp_label == 1) ++result.tumor_hp2_count;
                    } else {
                        if (hp_label == 0) ++result.normal_hp1_count;
                        else if (hp_label == 1) ++result.normal_hp2_count;
                    }
                }

                LabelDeltaResult tumor_hp = extract_and_test_hp(tumor_indices);
                LabelDeltaResult normal_hp = extract_and_test_hp(normal_indices);

                // Record diagnostic validity
                result.tumor_hp_valid = tumor_hp.valid;
                result.normal_hp_valid = normal_hp.valid;
                if (tumor_hp.valid) result.tumor_hp_delta = tumor_hp.delta;
                if (normal_hp.valid) result.normal_hp_delta = normal_hp.delta;

                // Somatic HP ASM computation
                // Key insight: at somatic SNV sites, tumor reads are almost all in
                // one HP (somatic haplotagging), so tumor-only HP test is usually
                // invalid. This is EXPECTED biology, not a bug.
                if (tumor_hp.valid && normal_hp.valid) {
                    // Both valid — somatic HP ASM = tumor delta - normal delta
                    result.hp_residual_delta = tumor_hp.delta - normal_hp.delta;
                    result.hp_residual_p = tumor_hp.p_value;
                    result.hp_residual_sig = tumor_hp.significant;
                } else if (tumor_hp.valid) {
                    // Normal HP test invalid (rare) — tumor HP is somatic
                    result.hp_residual_delta = tumor_hp.delta;
                    result.hp_residual_p = tumor_hp.p_value;
                    result.hp_residual_sig = tumor_hp.significant;
                } else if (normal_hp.valid) {
                    // Tumor HP test invalid (common: single HP dominance)
                    // Normal HP delta serves as germline baseline reference.
                    // The full-matrix HP test (hp_merged_delta) already captures
                    // the combined signal; here we record normal as negative
                    // to indicate "germline HP ASM present, tumor too skewed to test"
                    result.hp_residual_delta = -normal_hp.delta;
                    result.hp_residual_p = normal_hp.p_value;
                    result.hp_residual_sig = false;  // not a somatic signal
                }

                // Signed HP mean methylation delta: mean(HP1_meth) - mean(HP2_meth)
                // Captures direction of ASM unlike distance-based delta
                auto compute_signed_hp_delta = [&](const std::vector<int>& indices) -> double {
                    double hp1_sum = 0.0, hp2_sum = 0.0;
                    int hp1_count = 0, hp2_count = 0;
                    for (int idx : indices) {
                        int hp = hp_merged_labels[idx];
                        if (hp != 0 && hp != 1) continue;
                        for (int j = 0; j < meth_mat.raw_matrix.cols(); ++j) {
                            double val = meth_mat.raw_matrix(idx, j);
                            if (std::isnan(val)) continue;
                            if (hp == 0) {
                                hp1_sum += val;
                                ++hp1_count;
                            } else {
                                hp2_sum += val;
                                ++hp2_count;
                            }
                        }
                    }
                    if (hp1_count < 1 || hp2_count < 1) return NAN;
                    return (hp1_sum / hp1_count) - (hp2_sum / hp2_count);
                };

                result.tumor_hp_signed_delta = compute_signed_hp_delta(tumor_indices);
                result.normal_hp_signed_delta = compute_signed_hp_delta(normal_indices);

                // Combined: all reads (tumor + normal)
                std::vector<int> all_indices;
                all_indices.reserve(tumor_indices.size() + normal_indices.size());
                all_indices.insert(all_indices.end(), tumor_indices.begin(), tumor_indices.end());
                all_indices.insert(all_indices.end(), normal_indices.begin(), normal_indices.end());
                result.combined_hp_signed_delta = compute_signed_hp_delta(all_indices);

                // Signed residual: somatic directional ASM change
                if (!std::isnan(result.tumor_hp_signed_delta) && !std::isnan(result.normal_hp_signed_delta)) {
                    result.hp_signed_residual = result.tumor_hp_signed_delta - result.normal_hp_signed_delta;
                }

                // [Δβ module] read-level somatic residual Δβ + permutation test (#3 hp_residual 修法)
                // min_group=3: a (sample×HP) or (germline/carrier) cell must have >=3 reads for sig (tiny-group guard).
                const int kDbetaMinGroup = 3;
                compute_somatic_residual_dbeta_test(meth_mat.raw_matrix, hp_merged_labels, is_tumor_mask, 999, 42,
                                                    kDbetaMinGroup, result.somatic_residual_dbeta,
                                                    result.somatic_residual_dbeta_p, result.somatic_residual_dbeta_sig);

                // [Δβ module stage 2] germline ASM Δβ (normal) + fine same-hap subclone Δβ (tumor germline vs carrier).
                // Fine HP labels: 0=HP1 germline, 1=HP1-1 carrier, 2=HP2 germline, 3=HP2-1 carrier, -1=unlabeled.
                std::vector<int> hp_fine_labels(read_list.size(), -1);
                for (size_t i = 0; i < read_list.size(); ++i) {
                    const std::string& hp = read_list[i].hp_tag;
                    if (hp == "1" || hp == "HP1") {
                        hp_fine_labels[i] = 0;
                    } else if (hp == "1-1") {
                        hp_fine_labels[i] = 1;
                    } else if (hp == "2" || hp == "HP2") {
                        hp_fine_labels[i] = 2;
                    } else if (hp == "2-1") {
                        hp_fine_labels[i] = 3;
                    }
                }

                // Germline ASM Δβ (goal 1): normal reads HP1-family vs HP2-family (normal has no carriers → germline-only).
                std::vector<int> germ_grp(read_list.size(), -1);
                for (size_t i = 0; i < read_list.size(); ++i) {
                    if (!read_list[i].is_tumor && hp_merged_labels[i] >= 0) germ_grp[i] = hp_merged_labels[i];
                }
                int germ_n0 = 0, germ_n1 = 0;  // normal HP1/HP2 counts also in Normal_HP1/HP2; throwaway here
                compute_group_dbeta_test(meth_mat.raw_matrix, germ_grp, 999, 42, kDbetaMinGroup,
                                         result.germline_asm_dbeta, result.germline_asm_dbeta_p,
                                         result.germline_asm_dbeta_sig, germ_n0, germ_n1);

                // Fine subclone HP1 Δβ (goal 3): tumor HP1 germline vs HP1-1 somatic-carrier (same-haplotype subclone).
                std::vector<int> sub1_grp(read_list.size(), -1);
                for (size_t i = 0; i < read_list.size(); ++i) {
                    if (!read_list[i].is_tumor) continue;
                    if (hp_fine_labels[i] == 0) {
                        sub1_grp[i] = 0;
                    } else if (hp_fine_labels[i] == 1) {
                        sub1_grp[i] = 1;
                    }
                }
                compute_group_dbeta_test(meth_mat.raw_matrix, sub1_grp, 999, 42, kDbetaMinGroup,
                                         result.subclone_dbeta_hp1, result.subclone_dbeta_hp1_p,
                                         result.subclone_dbeta_hp1_sig, result.subclone_dbeta_hp1_n_germ,
                                         result.subclone_dbeta_hp1_n_carrier);

                // Fine subclone HP2 Δβ (goal 3): tumor HP2 germline vs HP2-1 somatic-carrier (same-haplotype subclone).
                std::vector<int> sub2_grp(read_list.size(), -1);
                for (size_t i = 0; i < read_list.size(); ++i) {
                    if (!read_list[i].is_tumor) continue;
                    if (hp_fine_labels[i] == 2) {
                        sub2_grp[i] = 0;
                    } else if (hp_fine_labels[i] == 3) {
                        sub2_grp[i] = 1;
                    }
                }
                compute_group_dbeta_test(meth_mat.raw_matrix, sub2_grp, 999, 42, kDbetaMinGroup,
                                         result.subclone_dbeta_hp2, result.subclone_dbeta_hp2_p,
                                         result.subclone_dbeta_hp2_sig, result.subclone_dbeta_hp2_n_germ,
                                         result.subclone_dbeta_hp2_n_carrier);

                // [Δβ component A] alt-axis subclone: within tumor HP-family, split by the read's OWN
                // alt_support (ALT=0 / REF=1) — phasing-independent somatic contrast. HP-family membership
                // uses the HP tag (1/1-1 = HP1-family), the group split uses alt_support.
                auto build_alt_grp = [&](bool hp1_family) {
                    std::vector<int> grp(read_list.size(), -1);
                    for (size_t i = 0; i < read_list.size(); ++i) {
                        if (!read_list[i].is_tumor) continue;
                        const std::string& h = read_list[i].hp_tag;
                        const bool in_family = hp1_family ? (h == "1" || h == "HP1" || h == "1-1")
                                                          : (h == "2" || h == "HP2" || h == "2-1");
                        if (!in_family) continue;
                        if (read_list[i].alt_support == AltSupport::ALT) {
                            grp[i] = 0;  // carries somatic allele
                        } else if (read_list[i].alt_support == AltSupport::REF) {
                            grp[i] = 1;
                        }
                    }
                    return grp;
                };
                std::vector<int> alt1_grp = build_alt_grp(true);
                compute_group_dbeta_test(meth_mat.raw_matrix, alt1_grp, 999, 42, kDbetaMinGroup,
                                         result.alt_subclone_hp1_dbeta, result.alt_subclone_hp1_p,
                                         result.alt_subclone_hp1_sig, result.alt_subclone_hp1_n_alt,
                                         result.alt_subclone_hp1_n_ref);
                std::vector<int> alt2_grp = build_alt_grp(false);
                compute_group_dbeta_test(meth_mat.raw_matrix, alt2_grp, 999, 42, kDbetaMinGroup,
                                         result.alt_subclone_hp2_dbeta, result.alt_subclone_hp2_p,
                                         result.alt_subclone_hp2_sig, result.alt_subclone_hp2_n_alt,
                                         result.alt_subclone_hp2_n_ref);

                // [Δβ component C] subclone per-CpG localization: per CpG, tumor germline vs carrier mean β;
                // a CpG is a "driver" if |Δβ_cpg|>0.2 with both groups >= min_group reads AT that CpG.
                auto count_driver_cpgs = [&](int germ_fine, int carrier_fine, int& out_n, int& out_tested) {
                    out_n = 0;
                    out_tested = 0;
                    const Eigen::MatrixXd& raw = meth_mat.raw_matrix;
                    for (int c = 0; c < raw.cols(); ++c) {
                        double gs = 0.0, cs = 0.0;
                        int gc = 0, cc = 0;
                        for (size_t i = 0; i < read_list.size(); ++i) {
                            if (!read_list[i].is_tumor) continue;
                            const double v = raw(static_cast<int>(i), c);
                            if (std::isnan(v)) continue;
                            if (hp_fine_labels[i] == germ_fine) {
                                gs += v;
                                ++gc;
                            } else if (hp_fine_labels[i] == carrier_fine) {
                                cs += v;
                                ++cc;
                            }
                        }
                        if (gc >= kDbetaMinGroup && cc >= kDbetaMinGroup) {
                            ++out_tested;
                            if (std::abs(gs / gc - cs / cc) > 0.2) ++out_n;
                        }
                    }
                };
                count_driver_cpgs(0, 1, result.subclone_hp1_driver_cpg_n, result.subclone_hp1_driver_cpg_tested);
                count_driver_cpgs(2, 3, result.subclone_hp2_driver_cpg_n, result.subclone_hp2_driver_cpg_tested);

                // [Δβ component B] full label-combination Δβ table (sample×HP-family×alt pairwise + BH-FDR).
                compute_combo_dbeta(meth_mat.raw_matrix, read_list, kDbetaMinGroup, 42,
                                    result.combo_dbeta_n_tested, result.combo_dbeta_n_sig,
                                    result.combo_dbeta_sig_pairs);
            }

            // [Stage② within-HP substructure] does a single germline-HP (tumor reads) split into >=2 clean
            // methylation clusters? Pattern-based (per-CpG distance from all_dist), independent of normal.
            {
                std::vector<int> hp1_idx, hp2_idx;
                for (size_t i = 0; i < read_list.size(); ++i) {
                    if (!read_list[i].is_tumor) continue;
                    const std::string& h = read_list[i].hp_tag;
                    if (h == "1" || h == "HP1" || h == "1-1") {
                        hp1_idx.push_back(static_cast<int>(i));
                    } else if (h == "2" || h == "HP2" || h == "2-1") {
                        hp2_idx.push_back(static_cast<int>(i));
                    }
                }
                double sil1 = 0.0, sil2 = 0.0, mf1 = 1.0, mf2 = 1.0;
                compute_within_hp_substructure(all_dist.dist_matrix, hp1_idx, result.within_hp1_ngroups, sil1, mf1);
                compute_within_hp_substructure(all_dist.dist_matrix, hp2_idx, result.within_hp2_ngroups, sil2, mf2);
                // [D] graded output: best within-HP silhouette + its smallest-cluster fraction (exposes weak/minor
                // splits the binary clean_multigroup gate drops; downstream/human can re-threshold).
                if (sil1 >= sil2) {
                    result.within_hp_best_sil = sil1;
                    result.within_hp_min_frac = mf1;
                } else {
                    result.within_hp_best_sil = sil2;
                    result.within_hp_min_frac = mf2;
                }
                // [Stage② level axis] distance clustering above catches PATTERN substructure only; add the LEVEL
                // axis (mean-β bimodality) which distance misses (chr1: distance 1.2% vs level 26.1%).
                constexpr int kWithinHpMinGroup = 3;
                result.within_hp_level_bimodal =
                    within_hp_level_clean(meth_mat.raw_matrix, hp1_idx, kWithinHpMinGroup) ||
                    within_hp_level_clean(meth_mat.raw_matrix, hp2_idx, kWithinHpMinGroup);
                result.within_hp_clean_multigroup = (result.within_hp1_ngroups >= 2 || result.within_hp2_ngroups >= 2 ||
                                                     result.within_hp_level_bimodal);
                // [C'] a-priori germline-vs-carrier within-HP PERMANOVA — statistical "mark/confirm": does the
                // within-HP distance separate germline(tag 1/2) vs somatic-carrier(tag 1-1/2-1) reads? a-priori
                // labels => no double-dip (unlike unsupervised tumor_intrinsic which fires ~96% on noise).
                std::vector<int> gc1, gc2;
                gc1.reserve(hp1_idx.size());
                gc2.reserve(hp2_idx.size());
                for (int idx : hp1_idx) gc1.push_back(read_list[idx].hp_tag == "1-1" ? 1 : 0);
                for (int idx : hp2_idx) gc2.push_back(read_list[idx].hp_tag == "2-1" ? 1 : 0);
                double scp1 = 1.0, scf1 = 0.0, scp2 = 1.0, scf2 = 0.0;
                bool scv1 = false, scv2 = false, scd1 = false, scd2 = false;
                compute_within_hp_subclone_permanova(all_dist.dist_matrix, hp1_idx, gc1, scp1, scf1, scv1, scd1);
                compute_within_hp_subclone_permanova(all_dist.dist_matrix, hp2_idx, gc2, scp2, scf2, scv2, scd2);
                const bool sig1 = scv1 && scp1 < 0.05 && !scd1;
                const bool sig2 = scv2 && scp2 < 0.05 && !scd2;
                result.within_hp_subclone_sig = sig1 || sig2;
                // report the more-significant valid HP
                if (scv1 && (!scv2 || scp1 <= scp2)) {
                    result.within_hp_subclone_permanova_p = scp1;
                    result.within_hp_subclone_permanova_f = scf1;
                    result.within_hp_subclone_valid = scv1;
                    result.within_hp_subclone_dispersion_warn = scd1;
                } else if (scv2) {
                    result.within_hp_subclone_permanova_p = scp2;
                    result.within_hp_subclone_permanova_f = scf2;
                    result.within_hp_subclone_valid = scv2;
                    result.within_hp_subclone_dispersion_warn = scd2;
                }
            }

            // [tumor-only structure axis] does structure exist in TUMOR reads alone? The all-pool main path cannot
            // self-distinguish a tumor-vs-normal split (sampleOrphan confound); this re-clusters TUMOR reads only and
            // gates with PERMANOVA permutation null. tumor_intrinsic = tumor alone has clean (non-dispersion) location.
            {
                std::vector<int> tumor_idx;
                tumor_idx.reserve(read_list.size());
                for (size_t i = 0; i < read_list.size(); ++i) {
                    if (read_list[i].is_tumor) tumor_idx.push_back(static_cast<int>(i));
                }
                compute_tumor_only_cluster_structure(all_dist.dist_matrix, tumor_idx, result.tumor_only_cluster_k,
                                                     result.tumor_only_silhouette, result.tumor_only_permanova_f,
                                                     result.tumor_only_permanova_p, result.tumor_only_permanova_valid,
                                                     result.tumor_only_dispersion_warn);
                result.tumor_intrinsic = result.tumor_only_permanova_valid && result.tumor_only_permanova_p < 0.05 &&
                                         !result.tumor_only_dispersion_warn;
            }

            // [Strength score] EQUAL-weighted mean of 5 sub-scores. Weight choice is empirically ranking-robust
            // (Spearman rho=0.998 equal-vs-judgment, top-1000 94.6% overlap); the saturated sig term (median 1.0,
            // non-discriminating) is dropped. germline demoted, somatic + tumor-only added (subclone-aligned). Each
            // sub-score is OUTPUT as its own column so strength is observable/verifiable/re-weightable downstream.
            {
                auto clamp01 = [](double v) { return std::max(0.0, std::min(1.0, v)); };
                result.strength_struct = clamp01((result.hp_auc_all - 0.5) / 0.5);  // -1 default -> 0
                // [fix] strength_tumor now gated by C' (within-HP germline-vs-carrier a-priori PERMANOVA): double-dip
                // REMOVED, ~2.4x enriched on legacy VC (NOT a sensitive detector, NOT "0% noise"). strength_tumor>0
                // rate: noise 7.2% vs structure 17.2% (all-VC); 12.7% vs 28.4% (within C'-valid coverage). Replaces
                // tumor_intrinsic (unsupervised cluster DERIVED from the same distance matrix it tests -> double-dip:
                // fired ~91% on noise). a-priori labels = read hp_tag (1 vs 1-1), not distance-derived. Magnitude =
                // within_hp_best_sil. tumor_only emit cols kept as raw diagnostics (observability), excluded from verdict.
                result.strength_tumor =
                    result.within_hp_subclone_sig ? clamp01(result.within_hp_best_sil / 0.5) : 0.0;
                result.strength_somatic = std::isnan(result.somatic_residual_dbeta)
                                              ? 0.0
                                              : clamp01(std::abs(result.somatic_residual_dbeta) / 0.3);
                result.strength_assoc = clamp01(result.cramers_v);
                result.strength_germline = std::isnan(result.germline_asm_dbeta)
                                               ? 0.0
                                               : clamp01(std::abs(result.germline_asm_dbeta) / 0.3);
                result.strength_score = (result.strength_struct + result.strength_tumor + result.strength_somatic +
                                         result.strength_assoc + result.strength_germline) /
                                        5.0;
                result.strength_grade = result.strength_score >= 0.50   ? "A"
                                        : result.strength_score >= 0.35 ? "B"
                                        : result.strength_score >= 0.22 ? "C"
                                        : result.strength_score >= 0.12 ? "D"
                                                                        : "E";
            }

            // Per-CpG ASM and epiallele metrics (Phase E)
            // Compute after distance matrix and HP analysis are done.
            // Uses raw_matrix + binary_matrix + HP labels from read_list.
            {
                PerCpgAsmResult pcr = compute_per_cpg_asm(meth_mat, read_list, 5);
                result.per_cpg_asm_valid = pcr.valid;
                if (pcr.valid) {
                    result.per_cpg_fisher_n_sig = pcr.fisher_n_sig;
                    result.per_cpg_fisher_frac_sig = pcr.fisher_frac_sig;
                    result.per_cpg_fisher_n_tested = pcr.fisher_n_tested;
                    result.per_cpg_fisher_max_neg_log_fdr = pcr.fisher_max_neg_log_fdr;
                    result.per_cpg_nme_hp1 = pcr.nme_hp1;
                    result.per_cpg_nme_hp2 = pcr.nme_hp2;
                    result.per_cpg_entropy_imbalance = pcr.entropy_imbalance;
                    result.per_cpg_epipoly_hp1 = pcr.epipoly_hp1;
                    result.per_cpg_epipoly_hp2 = pcr.epipoly_hp2;
                    result.per_cpg_epipoly_delta = pcr.epipoly_delta;
                }
            }
        }

        // Cleanup BAM records
        for (auto* r : reads) {
            bam_destroy1(r);
        }
        for (auto* r : normal_reads) {
            bam_destroy1(r);
        }

        // [Stage④ 新 VerificationClass] override with evidence-based reclassification (user philosophy:
        // retain if ANY structure). Computed HERE because Δβ/HP-AUC/Epipoly are only final after the HP +
        // Per-CpG blocks (the legacy 2x2 in SignificanceAnalyzer ran earlier, before Δβ). Original kept in
        // verification_class_legacy. within-HP multigroup + SEQC2-LOH deferred to Stage②③.
        {
            result.verification_class_legacy = result.verification_class;
            const bool dbeta_sig = result.germline_asm_dbeta_sig || result.subclone_dbeta_hp1_sig ||
                                   result.subclone_dbeta_hp2_sig || result.somatic_residual_dbeta_sig ||
                                   result.alt_subclone_hp1_sig || result.alt_subclone_hp2_sig;
            const bool hp_auc_struct = result.hp_auc_all >= 0.7;
            const std::string& lvc = result.verification_class_legacy;
            const bool cluster_match = (lvc == "Strong" || lvc == "Subclone");  // GlobalTest matched a label
            const bool loh_struct = result.potential_loh && (hp_auc_struct || dbeta_sig);
            const bool within_hp = result.within_hp_clean_multigroup;  // Stage② within-HP clean substructure

            // [D1] label-PERMANOVA significance split by PERMDISP (D2) into clean location-shift vs
            // dispersion-confounded. clean_location_permanova adds a valid structural axis; dispersion-only
            // significance gets its own DispersionStructure tier (structure exists but ambiguous attribution).
            const bool hp_perm_sig = result.label_hp_permanova_valid && result.label_hp_permanova_p < 0.05;
            const bool al_perm_sig = result.label_allele_permanova_valid && result.label_allele_permanova_p < 0.05;
            const bool clean_location_permanova = (hp_perm_sig && !result.label_hp_dispersion_warning) ||
                                                  (al_perm_sig && !result.label_allele_dispersion_warning);
            const bool dispersion_structure = (hp_perm_sig && result.label_hp_dispersion_warning) ||
                                              (al_perm_sig && result.label_allele_dispersion_warning);

            if (dbeta_sig || hp_auc_struct || cluster_match || loh_struct || within_hp || clean_location_permanova) {
                if (cluster_match) {
                    result.verification_class = "Strong";
                } else if (loh_struct) {
                    result.verification_class = "LOH-Structure";
                } else if (within_hp) {
                    result.verification_class = "MultiGroupNoLabel";  // S4/S6 within-HP clean multigroup
                } else if (dbeta_sig) {
                    result.verification_class = "LabelShift";  // Δβ mean-shift signal (dominant rescue)
                } else if (clean_location_permanova) {
                    result.verification_class = "PermanovaLocation";  // [D1] clean (PERMDISP-passed) PERMANOVA
                } else {
                    result.verification_class = "StructureNoLabel";  // HP-AUC structure, no label match
                }
            } else if (dispersion_structure) {
                result.verification_class = "DispersionStructure";  // [D1] dispersion-type structure (ambiguous)
            } else {
                const double ep = result.per_cpg_epipoly_hp1;
                const double pw = result.pairwise_mean_distance;
                if ((pw < 0.15) || (!std::isnan(ep) && ep < 0.2)) {
                    result.verification_class = "Noise_Uniform";  // 甲基均勻, 無可分
                } else if ((!std::isnan(ep) && ep > 0.5) || (pw > 0.35)) {
                    result.verification_class = "Noise_Chaotic";  // 高熵混亂無乾淨群
                } else {
                    result.verification_class = "Noise_Uncorrelated";
                }
            }
        }

        result.success = true;
    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = e.what();
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    return result;
}

void RegionProcessor::write_significance_summary(const std::vector<RegionResult>& results) const {
    if (results.empty()) return;

    // 1. Write CSV Summary
    std::string csv_path = output_dir_ + "/significance_summary.csv";
    std::ofstream csv_file(csv_path);
    if (!csv_file.is_open()) {
        LOG_ERROR("Failed to open significance summary CSV for writing: " + csv_path);
        return;
    }

    // Header (updated with Multi-Stage HP verification and Quality Assessment columns)
    csv_file << "RegionID,Chr,Pos,Ref,Alt,NumReads,NumCpGs,NReadsValid,HP_AUC_Normal,HP_AUC_Tumor,HP_AUC_All,GlobalP,CramersV,"
                "GlobalP_HPFamily,CramersV_HPFamily,GlobalP_HPFine,CramersV_HPFine,HPFine_NGroups_CF,"
                "HeuristicScore,StrengthScore,StrengthGrade,PassedGating,"
                "PairwiseMeanDist,PairwiseMedianDist,"
                "ClusterPermanovaF,ClusterPermanovaP,ClusterPermanovaValid,ClusterDispersionP,ClusterDispersionWarn,"
                // Stage 1: HP Merged
                "HPMergedDelta,HPMergedP,HPMergedSig,HP1FamilyN,HP2FamilyN,"
                // Stage 2: HP Fine-Grained
                "HPFineF,HPFineP,HPFineSig,HPFineNGroups,"
                "HPFineN_HP1,HPFineN_HP1S,HPFineN_HP2,HPFineN_HP2S,"
                "HPFineD_HP1_HP1S,HPFineD_HP1_HP2,HPFineD_HP1_HP2S,"
                "HPFineD_HP1S_HP2,HPFineD_HP1S_HP2S,HPFineD_HP2_HP2S,"
                // Stage 3: Allele
                "AlleleDelta,AlleleP,AlleleSig,"
                "LabelHPPermanovaF,LabelHPPermanovaP,LabelHPPermanovaValid,LabelHPDispersionP,LabelHPDispersionWarn,"
                "LabelAllelePermanovaF,LabelAllelePermanovaP,LabelAllelePermanovaValid,LabelAlleleDispersionP,LabelAlleleDispersionWarn,"
                // Stage 4: Unassigned Affinity
                "UnassignedAffinity,UnassignedAffinityP,UnassignedDir,NHP3,NHP0,"
                // Summary
                "DominantLabel,Stability,VerificationClass,VerificationClass_Legacy,"
                "WithinHP1_NGroups,WithinHP2_NGroups,WithinHP_LevelBimodal,WithinHP_CleanMultigroup,"
                "WithinHP_BestSil,WithinHP_MinFrac,"
                "WithinHP_SubclonePermanovaP,WithinHP_SubclonePermanovaF,WithinHP_SubcloneValid,WithinHP_SubcloneDispWarn,WithinHP_SubcloneSig,"
                // Multi-Layer Validation Quality Assessment (NEW)
                "HP_Ratio,Potential_LOH,Coverage_Multiple,Diploid_Coverage_Used,Coverage_Category,LOH_Subtype,"
                "Quality_Score,Quality_Tier,"
                // Phase B: Normal Baseline + Sample ASM + Residual HP
                "NTumorReads,NNormalReads,"
                "SampleASM_Delta,SampleASM_P,SampleASM_Sig,SampleASM_NTumor,SampleASM_NNormal,"
                "NormalBaseline_Mean,NormalBaseline_Coverage,"
                "HP_Residual_Delta,HP_Residual_P,HP_Residual_Sig,"
                "Tumor_HP_Delta,Tumor_HP_Valid,Normal_HP_Delta,Normal_HP_Valid,"
                "Tumor_HP1,Tumor_HP2,Normal_HP1,Normal_HP2,"
                // Signed HP methylation delta (direction-aware)
                "Tumor_HP_Signed_Delta,Normal_HP_Signed_Delta,HP_Signed_Residual,Combined_HP_Signed_Delta,"
                "SomaticResidualDbeta,SomaticResidualDbeta_P,SomaticResidualDbeta_Sig,"
                "GermlineAsmDbeta,GermlineAsmDbeta_P,GermlineAsmDbeta_Sig,"
                "SubcloneDbeta_HP1,SubcloneDbeta_HP1_P,SubcloneDbeta_HP1_Sig,SubcloneDbeta_HP1_NGerm,"
                "SubcloneDbeta_HP1_NCarrier,"
                "SubcloneDbeta_HP2,SubcloneDbeta_HP2_P,SubcloneDbeta_HP2_Sig,SubcloneDbeta_HP2_NGerm,"
                "SubcloneDbeta_HP2_NCarrier,"
                // [Δβ component A] alt-axis subclone (within HP-family, ALT vs REF by read's own allele)
                "AltSubcloneDbeta_HP1,AltSubcloneDbeta_HP1_P,AltSubcloneDbeta_HP1_Sig,"
                "AltSubcloneDbeta_HP1_NAlt,AltSubcloneDbeta_HP1_NRef,"
                "AltSubcloneDbeta_HP2,AltSubcloneDbeta_HP2_P,AltSubcloneDbeta_HP2_Sig,"
                "AltSubcloneDbeta_HP2_NAlt,AltSubcloneDbeta_HP2_NRef,"
                // [Δβ component C] subclone per-CpG driver localization
                "SubcloneHP1_DriverCpG_N,SubcloneHP1_DriverCpG_Tested,"
                "SubcloneHP2_DriverCpG_N,SubcloneHP2_DriverCpG_Tested,"
                // [Δβ component B] full label-combination Δβ (pairwise + BH-FDR)
                "ComboDbeta_NTested,ComboDbeta_NSig,ComboDbeta_SigPairs,"
                // Phase C: LOH BED Annotation
                "LOH_Bed_Overlap,LOH_Source,LOH_Bed_Annotation,"
                // Phase D: Subclone Assignment
                "Subclone_ID,"
                // Phase E: Per-CpG ASM + Epiallele Metrics
                "PerCpgASM_Valid,"
                "Fisher_N_Sig,Fisher_Frac_Sig,Fisher_N_Tested,Fisher_MaxNegLogFDR,"
                "NME_HP1,NME_HP2,Entropy_Imbalance,"
                "Epipoly_HP1,Epipoly_HP2,Epipoly_Delta,"
                // Original columns
                "LocalBestCluster,LocalBestP,Significant,SuggestFilter,"
                // Self-phasing audit (NEW): raw somatic HP tag counts (independent of --germline-hp-only flag)
                "NHP_Somatic11,NHP_Somatic21,NHP_Somatic33,"
                // tumor-only structure axis (clustering + PERMANOVA on TUMOR reads only) + StrengthScore components
                "TumorOnlyClusterK,TumorOnlySilhouette,TumorOnlyPermanovaF,TumorOnlyPermanovaP,"
                "TumorOnlyPermanovaValid,TumorOnlyDispersionWarn,TumorIntrinsic,"
                "StrengthStruct,StrengthTumor,StrengthSomatic,StrengthAssoc,StrengthGermline\n";

    // Statistics conatiners
    struct ChrStats {
        int total = 0;
        int significant = 0;
    };
    std::map<std::string, ChrStats> chr_stats;
    int total_significant = 0;
    int total_analyzed = 0;  // successfully processed and significance computed

    for (const auto& r : results) {
        if (!r.success || !r.significance_computed) continue;

        total_analyzed++;
        const auto& snv = snvs_[r.region_id];
        std::string chr_name = chrom_index_.get_name(snv.chr_id);

        // Determine significance (Binary classification)
        // [D1] Unified structural verdict: Significant = VerificationClass is a clean location-structure
        // class. Replaces the obsolete F1-filter gate (passed_gating && p && CramersV && depth) — the
        // methylation F1-filter direction is concluded DEAD, so there is one unified verdict and Significant
        // is its boolean projection. DispersionStructure (dispersion-type, ambiguous) and Noise_* are NOT
        // location-significant. Per-axis raw numbers (PERMANOVA/dispersion/Δβ/CramersV/HP-AUC/Strength) stay.
        const std::string& vc = r.verification_class;
        bool is_significant = (vc == "Strong" || vc == "LOH-Structure" || vc == "MultiGroupNoLabel" ||
                               vc == "LabelShift" || vc == "PermanovaLocation" || vc == "StructureNoLabel");
        
        // Suggest filtering sites with LabelDelta > 0.3 (F1 optimization)
        bool suggest_filter = (r.label_delta > 0.3);

        if (is_significant) {
            total_significant++;
            chr_stats[chr_name].significant++;
        }
        chr_stats[chr_name].total++;

        csv_file << r.region_id << "," << chr_name << "," << snv.pos << "," << snv.ref_base << "," << snv.alt_base
                 << "," << r.num_reads << "," << r.num_cpgs << "," << r.num_reads_valid << "," << std::fixed
                 << std::setprecision(4) << r.hp_auc_normal << "," << r.hp_auc_tumor << "," << r.hp_auc_all << ","
                 << std::scientific << std::setprecision(6) << r.global_p_value << "," << std::fixed
                 << std::setprecision(4) << r.cramers_v << ","
                 // Multi-layer HP (Cluster-First)
                 << std::scientific << std::setprecision(6) << r.global_hp_family_p << ","
                 << std::fixed << std::setprecision(4) << r.global_hp_family_v << ","
                 << std::scientific << std::setprecision(6) << r.global_hp_fine_p << ","
                 << std::fixed << std::setprecision(4) << r.global_hp_fine_v << ","
                 << r.global_hp_fine_n_groups << ","
                 << std::fixed << std::setprecision(4) << r.heuristic_score << ","
                 << std::fixed << std::setprecision(4) << r.strength_score << "," << r.strength_grade << ","
                 << (r.passed_gating ? "true" : "false") << ","
                 << std::fixed << std::setprecision(4) << r.pairwise_mean_distance << ","
                 << std::fixed << std::setprecision(4) << r.pairwise_median_distance << ","
                 << std::fixed << std::setprecision(4) << r.cluster_permanova_f << ","
                 << std::scientific << std::setprecision(6) << r.cluster_permanova_p << ","
                 << (r.cluster_permanova_valid ? "true" : "false") << ","
                 << std::scientific << std::setprecision(6) << r.cluster_dispersion_p << ","
                 << (r.cluster_dispersion_warning ? "true" : "false") << ","
                 // Stage 1: HP Merged
                 << std::fixed << std::setprecision(4) << r.hp_merged_delta << ","
                 << std::scientific << std::setprecision(6) << r.hp_merged_p << ","
                 << (r.hp_merged_sig ? "true" : "false") << ","
                 << r.hp_merged_n_hp1 << "," << r.hp_merged_n_hp2 << ","
                 // Stage 2: HP Fine-Grained
                 << std::fixed << std::setprecision(4) << r.hp_fine_f << ","
                 << std::scientific << std::setprecision(6) << r.hp_fine_p << ","
                 << (r.hp_fine_sig ? "true" : "false") << ","
                 << r.hp_fine_n_groups << ","
                 // Stage 2 cont: Fine group counts + pairwise distances
                 << r.hp_fine_n_hp1 << "," << r.hp_fine_n_hp1s << ","
                 << r.hp_fine_n_hp2 << "," << r.hp_fine_n_hp2s << ","
                 << std::fixed << std::setprecision(4)
                 << r.hp_fine_d_hp1_hp1s << "," << r.hp_fine_d_hp1_hp2 << ","
                 << r.hp_fine_d_hp1_hp2s << "," << r.hp_fine_d_hp1s_hp2 << ","
                 << r.hp_fine_d_hp1s_hp2s << "," << r.hp_fine_d_hp2_hp2s << ","
                 // Stage 3: Allele
                 << std::fixed << std::setprecision(4) << r.allele_delta << ","
                 << std::scientific << std::setprecision(6) << r.allele_p << ","
                 << (r.allele_sig ? "true" : "false") << ","
                 << std::fixed << std::setprecision(4) << r.label_hp_permanova_f << ","
                 << std::scientific << std::setprecision(6) << r.label_hp_permanova_p << ","
                 << (r.label_hp_permanova_valid ? "true" : "false") << ","
                 << std::scientific << std::setprecision(6) << r.label_hp_dispersion_p << ","
                 << (r.label_hp_dispersion_warning ? "true" : "false") << ","
                 << std::fixed << std::setprecision(4) << r.label_allele_permanova_f << ","
                 << std::scientific << std::setprecision(6) << r.label_allele_permanova_p << ","
                 << (r.label_allele_permanova_valid ? "true" : "false") << ","
                 << std::scientific << std::setprecision(6) << r.label_allele_dispersion_p << ","
                 << (r.label_allele_dispersion_warning ? "true" : "false") << ","
                 // Stage 4: Unassigned Affinity
                 << std::fixed << std::setprecision(4) << r.unassigned_affinity << ","
                 << std::scientific << std::setprecision(6) << r.unassigned_affinity_p << ","
                 << r.unassigned_dir << ","
                 << r.unassigned_n_hp3 << "," << r.unassigned_n_hp0 << ","
                 // Summary
                 << r.dominant_label << ","
                 << std::fixed << std::setprecision(4) << r.cluster_stability << ","
                 << r.verification_class << "," << r.verification_class_legacy << ","
                 << r.within_hp1_ngroups << "," << r.within_hp2_ngroups << ","
                 << (r.within_hp_level_bimodal ? "true" : "false") << ","
                 << (r.within_hp_clean_multigroup ? "true" : "false") << ","
                 << std::fixed << std::setprecision(4) << r.within_hp_best_sil << ","
                 << r.within_hp_min_frac << ","
                 << std::scientific << std::setprecision(6) << r.within_hp_subclone_permanova_p << ","
                 << std::fixed << std::setprecision(4) << r.within_hp_subclone_permanova_f << ","
                 << (r.within_hp_subclone_valid ? "true" : "false") << ","
                 << (r.within_hp_subclone_dispersion_warn ? "true" : "false") << ","
                 << (r.within_hp_subclone_sig ? "true" : "false") << ","
                 // Multi-Layer Validation Quality Assessment (NEW)
                 << std::fixed << std::setprecision(4) << r.hp_ratio << ","
                 << (r.potential_loh ? "true" : "false") << ","
                 << std::fixed << std::setprecision(3) << r.coverage_multiple << ","
                 << std::fixed << std::setprecision(2) << r.diploid_coverage_used << ","
                 << r.coverage_category << ","
                 << r.loh_subtype << ","
                 << std::fixed << std::setprecision(1) << r.quality_score << ","
                 << r.quality_tier << ","
                 // Phase B: Normal Baseline + Sample ASM + Residual HP
                 << r.num_tumor_reads << "," << r.num_normal_reads << ","
                 << std::fixed << std::setprecision(6) << r.sample_asm_delta << ","
                 << std::scientific << std::setprecision(6) << r.sample_asm_p << ","
                 << (r.sample_asm_sig ? "true" : "false") << ","
                 << r.sample_asm_n_tumor << "," << r.sample_asm_n_normal << ","
                 << std::fixed << std::setprecision(4) << r.normal_baseline_mean << ","
                 << std::fixed << std::setprecision(2) << r.normal_baseline_coverage << ","
                 << std::fixed << std::setprecision(6) << r.hp_residual_delta << ","
                 << std::scientific << std::setprecision(6) << r.hp_residual_p << ","
                 << (r.hp_residual_sig ? "true" : "false") << ","
                 // HP subset diagnostics
                 << (std::isnan(r.tumor_hp_delta) ? "NA" : std::to_string(r.tumor_hp_delta)) << ","
                 << (r.tumor_hp_valid ? "true" : "false") << ","
                 << (std::isnan(r.normal_hp_delta) ? "NA" : std::to_string(r.normal_hp_delta)) << ","
                 << (r.normal_hp_valid ? "true" : "false") << ","
                 << r.tumor_hp1_count << "," << r.tumor_hp2_count << ","
                 << r.normal_hp1_count << "," << r.normal_hp2_count << ","
                 // Signed HP methylation delta (direction-aware)
                 << (std::isnan(r.tumor_hp_signed_delta) ? "NA" : std::to_string(r.tumor_hp_signed_delta)) << ","
                 << (std::isnan(r.normal_hp_signed_delta) ? "NA" : std::to_string(r.normal_hp_signed_delta)) << ","
                 << (std::isnan(r.hp_signed_residual) ? "NA" : std::to_string(r.hp_signed_residual)) << ","
                 << (std::isnan(r.combined_hp_signed_delta) ? "NA" : std::to_string(r.combined_hp_signed_delta)) << ","
                 // [Δβ module] somatic residual Δβ + permutation test
                 << (std::isnan(r.somatic_residual_dbeta) ? "NA" : std::to_string(r.somatic_residual_dbeta)) << ","
                 << r.somatic_residual_dbeta_p << "," << (r.somatic_residual_dbeta_sig ? "true" : "false") << ","
                 // [Δβ module stage 2] germline ASM Δβ + fine same-hap subclone Δβ (HP1/HP2)
                 << (std::isnan(r.germline_asm_dbeta) ? "NA" : std::to_string(r.germline_asm_dbeta)) << ","
                 << r.germline_asm_dbeta_p << "," << (r.germline_asm_dbeta_sig ? "true" : "false") << ","
                 << (std::isnan(r.subclone_dbeta_hp1) ? "NA" : std::to_string(r.subclone_dbeta_hp1)) << ","
                 << r.subclone_dbeta_hp1_p << "," << (r.subclone_dbeta_hp1_sig ? "true" : "false") << ","
                 << r.subclone_dbeta_hp1_n_germ << "," << r.subclone_dbeta_hp1_n_carrier << ","
                 << (std::isnan(r.subclone_dbeta_hp2) ? "NA" : std::to_string(r.subclone_dbeta_hp2)) << ","
                 << r.subclone_dbeta_hp2_p << "," << (r.subclone_dbeta_hp2_sig ? "true" : "false") << ","
                 << r.subclone_dbeta_hp2_n_germ << "," << r.subclone_dbeta_hp2_n_carrier << ","
                 // [Δβ component A] alt-axis subclone
                 << (std::isnan(r.alt_subclone_hp1_dbeta) ? "NA" : std::to_string(r.alt_subclone_hp1_dbeta)) << ","
                 << r.alt_subclone_hp1_p << "," << (r.alt_subclone_hp1_sig ? "true" : "false") << ","
                 << r.alt_subclone_hp1_n_alt << "," << r.alt_subclone_hp1_n_ref << ","
                 << (std::isnan(r.alt_subclone_hp2_dbeta) ? "NA" : std::to_string(r.alt_subclone_hp2_dbeta)) << ","
                 << r.alt_subclone_hp2_p << "," << (r.alt_subclone_hp2_sig ? "true" : "false") << ","
                 << r.alt_subclone_hp2_n_alt << "," << r.alt_subclone_hp2_n_ref << ","
                 // [Δβ component C] subclone per-CpG driver localization
                 << r.subclone_hp1_driver_cpg_n << "," << r.subclone_hp1_driver_cpg_tested << ","
                 << r.subclone_hp2_driver_cpg_n << "," << r.subclone_hp2_driver_cpg_tested << ","
                 // [Δβ component B] full label-combination Δβ (SigPairs quoted: contains ; ~ =)
                 << r.combo_dbeta_n_tested << "," << r.combo_dbeta_n_sig << ","
                 << "\"" << r.combo_dbeta_sig_pairs << "\","
                 // Phase C: LOH BED Annotation
                 << (r.loh_bed_overlap ? "true" : "false") << ","
                 << r.loh_source << ","
                 << "\"" << r.loh_bed_annotation << "\","
                 // Phase D: Subclone Assignment
                 << r.subclone_id << ","
                 // Phase E: Per-CpG ASM + Epiallele Metrics
                 << (r.per_cpg_asm_valid ? "true" : "false") << ","
                 << r.per_cpg_fisher_n_sig << ","
                 << std::fixed << std::setprecision(4) << r.per_cpg_fisher_frac_sig << ","
                 << r.per_cpg_fisher_n_tested << ","
                 << std::fixed << std::setprecision(4) << r.per_cpg_fisher_max_neg_log_fdr << ","
                 << (std::isnan(r.per_cpg_nme_hp1) ? "NA" : std::to_string(r.per_cpg_nme_hp1)) << ","
                 << (std::isnan(r.per_cpg_nme_hp2) ? "NA" : std::to_string(r.per_cpg_nme_hp2)) << ","
                 << (std::isnan(r.per_cpg_entropy_imbalance) ? "NA" : std::to_string(r.per_cpg_entropy_imbalance)) << ","
                 << (std::isnan(r.per_cpg_epipoly_hp1) ? "NA" : std::to_string(r.per_cpg_epipoly_hp1)) << ","
                 << (std::isnan(r.per_cpg_epipoly_hp2) ? "NA" : std::to_string(r.per_cpg_epipoly_hp2)) << ","
                 << (std::isnan(r.per_cpg_epipoly_delta) ? "NA" : std::to_string(r.per_cpg_epipoly_delta)) << ","
                 // Original local test columns
                 << r.local_best_cluster << "," << std::scientific
                 << std::setprecision(6) << r.local_best_p_value << ","
                 << (is_significant ? "true" : "false") << ","
                 // SuggestFilter column for F1 optimization
                 << (suggest_filter ? "true" : "false") << ","
                 // Self-phasing audit: raw somatic HP tag counts (independent of --germline-hp-only flag)
                 << r.n_hp_somatic_11 << "," << r.n_hp_somatic_21 << "," << r.n_hp_somatic_33 << ","
                 // tumor-only structure axis + StrengthScore components
                 << r.tumor_only_cluster_k << "," << std::fixed << std::setprecision(4) << r.tumor_only_silhouette
                 << "," << std::setprecision(4) << r.tumor_only_permanova_f << "," << std::scientific
                 << std::setprecision(6) << r.tumor_only_permanova_p << ","
                 << (r.tumor_only_permanova_valid ? "true" : "false") << ","
                 << (r.tumor_only_dispersion_warn ? "true" : "false") << ","
                 << (r.tumor_intrinsic ? "true" : "false") << "," << std::fixed << std::setprecision(4)
                 << r.strength_struct << "," << r.strength_tumor << "," << r.strength_somatic << ","
                 << r.strength_assoc << "," << r.strength_germline << "\n";
    }
    csv_file.close();
    LOG_INFO("Significance summary written to: " + csv_path);

    // 2. Write Statistics Report
    std::string stats_path = output_dir_ + "/significance_statistics.txt";
    std::ofstream stats_file(stats_path);
    if (!stats_file.is_open()) {
        LOG_ERROR("Failed to open significance statistics file: " + stats_path);
        return;
    }

    stats_file << "=== Significance Analysis Statistics ===\n\n";
    stats_file << "Total Regions Processed: " << results.size() << "\n";
    stats_file << "Regions Analyzed (Success): " << total_analyzed << "\n";

    double sig_rate = total_analyzed > 0 ? (100.0 * total_significant / total_analyzed) : 0.0;
    stats_file << "Total Significant: " << total_significant << " (" << std::fixed << std::setprecision(2) << sig_rate
               << "%)\n";
    stats_file << "Significance Threshold: Passed Gating AND Global P-value <= 0.05\n\n";

    stats_file << "--- Per-Chromosome Breakdown ---\n";
    stats_file << std::left << std::setw(10) << "Chr" << std::setw(10) << "Total" << std::setw(15) << "Significant"
               << std::setw(10) << "% Sig" << "\n";

    for (const auto& [chr, stats] : chr_stats) {
        double chr_rate = stats.total > 0 ? (100.0 * stats.significant / stats.total) : 0.0;
        stats_file << std::left << std::setw(10) << chr << std::setw(10) << stats.total << std::setw(15)
                   << stats.significant << std::setw(10) << std::fixed << std::setprecision(2) << chr_rate << "%\n";
    }
    stats_file.close();
    LOG_INFO("Significance statistics summary written to: " + stats_path);
}

void RegionProcessor::print_summary(const std::vector<RegionResult>& results) const {
    int success_count = 0;
    int total_reads = 0;
    int total_cpgs = 0;
    int total_forward = 0;
    int total_reverse = 0;
    int total_filtered = 0;
    double total_time = 0.0;

    // Distance matrix statistics
    int total_valid_pairs = 0;
    int total_invalid_pairs = 0;
    double total_avg_coverage = 0.0;
    int regions_with_distance = 0;

    for (const auto& r : results) {
        if (r.success) {
            success_count++;
            total_reads += r.num_reads;
            total_cpgs += r.num_cpgs;
            total_forward += r.num_forward_reads;
            total_reverse += r.num_reverse_reads;
            total_filtered += r.num_filtered_reads;
            total_time += r.elapsed_ms;

            // Distance matrix stats
            if (r.num_valid_pairs > 0 || r.num_invalid_pairs > 0) {
                total_valid_pairs += r.num_valid_pairs;
                total_invalid_pairs += r.num_invalid_pairs;
                total_avg_coverage += r.avg_common_coverage;
                regions_with_distance++;
            }
        }
    }

    std::stringstream ss;
    ss << "\n=== Processing Summary ===\n"
       << "Total regions: " << results.size() << "\n"
       << "Successful: " << success_count << "\n"
       << "Failed: " << (results.size() - success_count) << "\n"
       << "Total reads processed: " << total_reads << "\n"
       << "  Forward strand (+): " << total_forward << "\n"
       << "  Reverse strand (-): " << total_reverse << "\n";
    if (output_filtered_reads_) {
        ss << "  Filtered out: " << total_filtered << "\n";
    }
    ss << "Total CpG sites found: " << total_cpgs << "\n"
       << "Total processing time: " << total_time << " ms\n"
       << "Average time per region: " << (results.size() > 0 ? (total_time / results.size()) : 0) << " ms\n"
       << "Average reads per region: " << (success_count > 0 ? (total_reads / static_cast<double>(success_count)) : 0)
       << "\n"
       << "Average CpGs per region: " << (success_count > 0 ? (total_cpgs / static_cast<double>(success_count)) : 0)
       << "\n";

    // Distance matrix summary
    if (compute_distance_matrix_ && regions_with_distance > 0) {
        ss << "\n=== Distance Matrix Summary (First Metric) ===\n"
           << "Metric: " << DistanceCalculator::metric_to_string(distance_metrics_[0]) << "\n"
           << "Min common coverage (C_min): " << distance_config_.min_common_coverage << "\n"
           << "Regions with distance matrices: " << regions_with_distance << "\n"
           << "Total valid read pairs: " << total_valid_pairs << "\n"
           << "Total invalid pairs (insufficient overlap): " << total_invalid_pairs << "\n";
        if (total_valid_pairs + total_invalid_pairs > 0) {
            double valid_ratio = 100.0 * total_valid_pairs / (total_valid_pairs + total_invalid_pairs);
            ss << "Valid pair ratio: " << std::fixed << std::setprecision(1) << valid_ratio << "%\n";
        }
        ss << "Average common CpG coverage: " << std::fixed << std::setprecision(2)
           << (total_avg_coverage / regions_with_distance) << "\n";
    }

    LOG_INFO(ss.str());
}

MethylationMatrix RegionProcessor::build_methylation_matrix(const MatrixBuilder& matrix_builder, int region_id) {
    MethylationMatrix meth_mat;
    meth_mat.region_id = region_id;

    const auto& read_list = matrix_builder.get_reads();
    const auto& raw_matrix = matrix_builder.get_matrix();
    const auto& cpg_positions = matrix_builder.get_cpg_positions();

    // Set up read IDs
    meth_mat.read_ids.resize(read_list.size());
    for (size_t i = 0; i < read_list.size(); ++i) {
        meth_mat.read_ids[i] = read_list[i].read_id;
    }

    // Set up CpG IDs
    meth_mat.cpg_ids.resize(cpg_positions.size());
    for (size_t i = 0; i < cpg_positions.size(); ++i) {
        meth_mat.cpg_ids[i] = static_cast<int>(i);
    }

    // Convert matrix to Eigen format
    int n_reads = static_cast<int>(raw_matrix.size());
    int n_cpgs = n_reads > 0 ? static_cast<int>(raw_matrix[0].size()) : 0;

    meth_mat.raw_matrix = Eigen::MatrixXd(n_reads, n_cpgs);
    meth_mat.binary_matrix = Eigen::MatrixXi(n_reads, n_cpgs);

    for (int i = 0; i < n_reads; ++i) {
        for (int j = 0; j < n_cpgs; ++j) {
            double val = raw_matrix[i][j];
            if (val < 0) {  // -1.0 indicates no coverage
                meth_mat.raw_matrix(i, j) = NAN;
                meth_mat.binary_matrix(i, j) = -1;
            } else {
                meth_mat.raw_matrix(i, j) = val;
                // Binary threshold (configurable via --methyl-high/--methyl-low)
                if (val >= binary_methyl_high_) {
                    meth_mat.binary_matrix(i, j) = 1;
                } else if (val <= binary_methyl_low_) {
                    meth_mat.binary_matrix(i, j) = 0;
                } else {
                    meth_mat.binary_matrix(i, j) = -1;  // Ambiguous
                }
            }
        }
    }

    return meth_mat;
}

DistanceMatrix RegionProcessor::compute_and_write_distance_matrix(const MethylationMatrix& meth_mat,
                                                                   const std::vector<ReadInfo>& read_list,
                                                                   const std::string& region_dir,
                                                                   DistanceMetricType metric, RegionResult& result) {
    DistanceConfig config = distance_config_;
    config.metric = metric;
    DistanceCalculator dist_calc(config);

    DistanceMatrix all_dist = dist_calc.compute(meth_mat, read_list);

    // Update result statistics
    result.num_valid_pairs = all_dist.num_valid_pairs;
    result.num_invalid_pairs = all_dist.num_invalid_pairs;
    result.avg_common_coverage = all_dist.avg_common_coverage;
    std::tie(result.pairwise_mean_distance, result.pairwise_median_distance) =
        summarize_pairwise_distances(all_dist.dist_matrix);

    // Write distance matrix if enabled
    if (output_distance_matrix_) {
        DistanceMatrix forward_dist, reverse_dist;
        if (output_strand_distance_matrices_) {
            std::tie(forward_dist, reverse_dist) = dist_calc.compute_strand_specific(meth_mat, read_list);
        }

        RegionWriter writer(output_dir_, debug_output_dir_, true, vcf_filename_);
        writer.write_distance_matrices(region_dir, all_dist, forward_dist, reverse_dist, metric,
                                        output_strand_distance_matrices_);
    }

    if (log_level_ >= LogLevel::LOG_TRACE) {
        std::stringstream ss;
        ss << "  Distance matrix (" << DistanceCalculator::metric_to_string(metric) << "): " << all_dist.size() << "x"
           << all_dist.size() << ", valid pairs: " << all_dist.num_valid_pairs
           << ", avg coverage: " << std::fixed << std::setprecision(1) << all_dist.avg_common_coverage;
        LOG_TRACE(ss.str());
    }

    return all_dist;
}

void RegionProcessor::extract_complete_submatrix_indices(const Eigen::MatrixXd& dist, std::vector<int>& out_indices) {
    // Greedy peel-off: drop the read with the fewest non-NaN partners until the
    // retained sub-matrix is NaN-free. Same algorithm as
    // StructureTest::filter_reads_for_complete_matrix, so the clustering subset
    // matches the subset PERMANOVA re-derives downstream.
    out_indices.clear();
    const int n = static_cast<int>(dist.rows());
    if (n == 0) return;

    std::vector<bool> keep(n, true);
    int n_keep = n;

    while (true) {
        bool has_nan = false;
        for (int i = 0; i < n && !has_nan; ++i) {
            if (!keep[i]) continue;
            for (int j = i + 1; j < n; ++j) {
                if (!keep[j]) continue;
                if (std::isnan(dist(i, j))) {
                    has_nan = true;
                    break;
                }
            }
        }
        if (!has_nan) break;
        if (n_keep <= 2) {  // cannot form a clusterable complete sub-matrix
            out_indices.clear();
            return;
        }

        int worst = -1;
        int min_degree = n_keep;
        for (int i = 0; i < n; ++i) {
            if (!keep[i]) continue;
            int degree = 0;
            for (int j = 0; j < n; ++j) {
                if (i == j || !keep[j]) continue;
                if (!std::isnan(dist(i, j))) ++degree;
            }
            if (degree < min_degree) {
                min_degree = degree;
                worst = i;
            }
        }
        if (worst < 0) break;
        keep[worst] = false;
        --n_keep;
    }

    for (int i = 0; i < n; ++i) {
        if (keep[i]) out_indices.push_back(i);
    }
}

double RegionProcessor::compute_hp_auc(const Eigen::MatrixXd& dist, const std::vector<int>& hp_fam,
                                       const std::vector<int>& idx) {
    // Collect same-HP and diff-HP pairwise distances (skip NaN = SKIP invalid pairs).
    std::vector<double> same, diff;
    for (size_t a = 0; a < idx.size(); ++a) {
        for (size_t b = a + 1; b < idx.size(); ++b) {
            const int i = idx[a], j = idx[b];
            const double d = dist(i, j);
            if (std::isnan(d)) continue;
            (hp_fam[i] == hp_fam[j] ? same : diff).push_back(d);
        }
    }
    if (same.empty() || diff.empty()) return -1.0;
    // Rank-based AUC = mean over diff of [#same<d + 0.5*#same==d] / n_same.
    std::sort(same.begin(), same.end());
    double wins = 0.0;
    for (const double d : diff) {
        const auto lo = std::lower_bound(same.begin(), same.end(), d);
        const auto hi = std::upper_bound(same.begin(), same.end(), d);
        wins += static_cast<double>(lo - same.begin()) + 0.5 * static_cast<double>(hi - lo);
    }
    return wins / (static_cast<double>(diff.size()) * static_cast<double>(same.size()));
}

void RegionProcessor::compute_somatic_residual_dbeta_test(const Eigen::MatrixXd& raw,
                                                          const std::vector<int>& hp_fam,
                                                          const std::vector<bool>& is_tumor, int n_perm,
                                                          uint64_t seed, int min_group, double& out_dbeta,
                                                          double& out_p, bool& out_sig) {
    out_dbeta = std::numeric_limits<double>::quiet_NaN();
    out_p = 1.0;
    out_sig = false;
    const int n = static_cast<int>(raw.rows());

    // Per-read mean β over valid (non-NaN) CpG. read-level so permutation is O(reads)/shuffle.
    std::vector<double> read_mean(n, std::numeric_limits<double>::quiet_NaN());
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        int c = 0;
        for (int j = 0; j < raw.cols(); ++j) {
            const double v = raw(i, j);
            if (!std::isnan(v)) {
                s += v;
                ++c;
            }
        }
        if (c > 0) read_mean[i] = s / c;
    }

    // Collect HP-labeled reads (0=HP1f, 1=HP2f) with valid mean + sample(T/N).
    std::vector<double> mval;
    std::vector<int> hp;
    std::vector<char> tum;
    for (int i = 0; i < n; ++i) {
        if (hp_fam[i] < 0 || std::isnan(read_mean[i])) continue;
        mval.push_back(read_mean[i]);
        hp.push_back(hp_fam[i]);
        tum.push_back(is_tumor[i] ? 1 : 0);
    }

    // residual = (tumor: mean β HP1 - mean β HP2) - (normal: mean β HP1 - mean β HP2)
    auto residual = [&](const std::vector<char>& sample) -> double {
        double t1 = 0, t2 = 0, nn1 = 0, nn2 = 0;
        int ct1 = 0, ct2 = 0, cn1 = 0, cn2 = 0;
        for (size_t k = 0; k < mval.size(); ++k) {
            if (sample[k]) {
                if (hp[k] == 0) {
                    t1 += mval[k];
                    ++ct1;
                } else {
                    t2 += mval[k];
                    ++ct2;
                }
            } else {
                if (hp[k] == 0) {
                    nn1 += mval[k];
                    ++cn1;
                } else {
                    nn2 += mval[k];
                    ++cn2;
                }
            }
        }
        if (ct1 < 1 || ct2 < 1 || cn1 < 1 || cn2 < 1) return std::numeric_limits<double>::quiet_NaN();
        return (t1 / ct1 - t2 / ct2) - (nn1 / cn1 - nn2 / cn2);
    };

    const double obs = residual(tum);
    if (std::isnan(obs)) return;
    out_dbeta = obs;

    // Observed counts of the 4 (sample×HP) cells, for the tiny-group guard on sig.
    int g_t1 = 0, g_t2 = 0, g_n1 = 0, g_n2 = 0;
    for (size_t k = 0; k < mval.size(); ++k) {
        if (tum[k]) {
            (hp[k] == 0 ? g_t1 : g_t2)++;
        } else {
            (hp[k] == 0 ? g_n1 : g_n2)++;
        }
    }

    // Permutation: shuffle sample(T/N) label among HP-labeled reads (HP fixed). Null = tumor/normal share HP-ASM.
    std::mt19937_64 rng(seed);
    std::vector<char> perm = tum;
    int n_extreme = 1;
    for (int p = 0; p < n_perm; ++p) {
        std::shuffle(perm.begin(), perm.end(), rng);
        const double r = residual(perm);
        if (!std::isnan(r) && std::abs(r) >= std::abs(obs)) ++n_extreme;
    }
    out_p = static_cast<double>(n_extreme) / static_cast<double>(n_perm + 1);
    // Tiny-group guard: every (sample×HP) cell must have >= min_group reads.
    out_sig = (out_p <= 0.05) && (std::min(std::min(g_t1, g_t2), std::min(g_n1, g_n2)) >= min_group);
}

void RegionProcessor::compute_group_dbeta_test(const Eigen::MatrixXd& raw, const std::vector<int>& group, int n_perm,
                                               uint64_t seed, int min_group, double& out_dbeta, double& out_p,
                                               bool& out_sig, int& out_n0, int& out_n1) {
    out_dbeta = std::numeric_limits<double>::quiet_NaN();
    out_p = 1.0;
    out_sig = false;
    out_n0 = 0;
    out_n1 = 0;
    const int n = static_cast<int>(raw.rows());

    // Collect group-labeled reads (0/1) with valid per-read mean β. read-level → permutation O(reads)/shuffle.
    std::vector<double> mval;
    std::vector<char> grp;
    for (int i = 0; i < n; ++i) {
        if (group[i] != 0 && group[i] != 1) continue;
        double s = 0.0;
        int c = 0;
        for (int j = 0; j < raw.cols(); ++j) {
            const double v = raw(i, j);
            if (!std::isnan(v)) {
                s += v;
                ++c;
            }
        }
        if (c == 0) continue;
        mval.push_back(s / c);
        grp.push_back(static_cast<char>(group[i]));
        (group[i] == 0 ? out_n0 : out_n1)++;
    }

    // Δβ = mean β(group 0) − mean β(group 1).
    auto dbeta = [&](const std::vector<char>& g) -> double {
        double a = 0, b = 0;
        int ca = 0, cb = 0;
        for (size_t k = 0; k < mval.size(); ++k) {
            if (g[k] == 0) {
                a += mval[k];
                ++ca;
            } else {
                b += mval[k];
                ++cb;
            }
        }
        if (ca < 1 || cb < 1) return std::numeric_limits<double>::quiet_NaN();
        return a / ca - b / cb;
    };

    const double obs = dbeta(grp);
    if (std::isnan(obs)) return;
    out_dbeta = obs;

    // Permutation: shuffle the 0/1 group label among labeled reads. Null = label exchangeable.
    std::mt19937_64 rng(seed);
    std::vector<char> perm = grp;
    int n_extreme = 1;
    for (int p = 0; p < n_perm; ++p) {
        std::shuffle(perm.begin(), perm.end(), rng);
        const double r = dbeta(perm);
        if (!std::isnan(r) && std::abs(r) >= std::abs(obs)) ++n_extreme;
    }
    out_p = static_cast<double>(n_extreme) / static_cast<double>(n_perm + 1);
    // Tiny-group guard: require both groups >= min_group so a 1-read group cannot drive a "significant" subclone.
    out_sig = (out_p <= 0.05) && (std::min(out_n0, out_n1) >= min_group);
}

void RegionProcessor::compute_combo_dbeta(const Eigen::MatrixXd& raw, const std::vector<ReadInfo>& read_list,
                                          int min_group, uint64_t seed, int& out_n_tested, int& out_n_sig,
                                          std::string& out_sig_pairs) {
    out_n_tested = 0;
    out_n_sig = 0;
    out_sig_pairs.clear();

    // Build groups keyed by sample(N/T) × HP-family(HP1/HP2) × alt(REF/ALT); keep groups with >= min_group reads.
    struct Group {
        std::string name;
        std::vector<int> idx;
    };
    std::vector<Group> groups;
    const char* sample_name[2] = {"N", "T"};
    const char* fam_name[2] = {"HP1", "HP2"};
    const char* alt_name[2] = {"REF", "ALT"};
    for (int s = 0; s < 2; ++s) {
        for (int f = 0; f < 2; ++f) {
            for (int a = 0; a < 2; ++a) {
                Group g;
                g.name = std::string(sample_name[s]) + "." + fam_name[f] + "." + alt_name[a];
                for (size_t i = 0; i < read_list.size(); ++i) {
                    if ((s == 1) != read_list[i].is_tumor) continue;
                    const std::string& h = read_list[i].hp_tag;
                    const bool in_family = (f == 0) ? (h == "1" || h == "HP1" || h == "1-1")
                                                    : (h == "2" || h == "HP2" || h == "2-1");
                    if (!in_family) continue;
                    const AltSupport as = read_list[i].alt_support;
                    if ((a == 1 && as != AltSupport::ALT) || (a == 0 && as != AltSupport::REF)) continue;
                    g.idx.push_back(static_cast<int>(i));
                }
                if (static_cast<int>(g.idx.size()) >= min_group) groups.push_back(std::move(g));
            }
        }
    }

    // All pairwise read-level Δβ (reuse compute_group_dbeta_test).
    const int n = static_cast<int>(raw.rows());
    struct Pair {
        std::string g1, g2;
        double dbeta, p;
    };
    std::vector<Pair> pairs;
    for (size_t i = 0; i < groups.size(); ++i) {
        for (size_t j = i + 1; j < groups.size(); ++j) {
            std::vector<int> grp(n, -1);
            for (int k : groups[i].idx) grp[k] = 0;
            for (int k : groups[j].idx) grp[k] = 1;
            double db = 0, p = 1;
            bool sig = false;
            int n0 = 0, n1 = 0;
            compute_group_dbeta_test(raw, grp, 999, seed, min_group, db, p, sig, n0, n1);
            if (std::isnan(db)) continue;
            pairs.push_back({groups[i].name, groups[j].name, db, p});
        }
    }
    out_n_tested = static_cast<int>(pairs.size());
    if (pairs.empty()) return;

    // BH-FDR within region.
    std::vector<int> order(pairs.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return pairs[a].p < pairs[b].p; });
    const int m = static_cast<int>(pairs.size());
    std::vector<double> q(pairs.size());
    double running_min = 1.0;
    for (int rank = m; rank >= 1; --rank) {
        const int idx = order[rank - 1];
        running_min = std::min(running_min, pairs[idx].p * m / rank);
        q[idx] = std::min(running_min, 1.0);
    }

    std::string out;
    for (size_t i = 0; i < pairs.size(); ++i) {
        if (q[i] <= 0.05) {
            ++out_n_sig;
            if (!out.empty()) out += ";";
            out += pairs[i].g1 + "~" + pairs[i].g2 + "=" + std::to_string(pairs[i].dbeta) + "(" +
                   std::to_string(q[i]) + ")";
        }
    }
    out_sig_pairs = out;
}

void RegionProcessor::compute_within_hp_substructure(const Eigen::MatrixXd& all_dist, const std::vector<int>& hp_idx,
                                                     int& out_ngroups, double& out_silhouette,
                                                     double& out_min_frac) const {
    out_ngroups = 1;
    out_silhouette = 0.0;
    out_min_frac = 1.0;  // [D] default = single group (no split)
    const int m = static_cast<int>(hp_idx.size());
    if (m < 8) return;  // need enough reads for a meaningful within-HP split

    // Extract the HP-family sub-distance-matrix from the already-computed all_dist (no recompute).
    Eigen::MatrixXd sub(m, m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) sub(i, j) = all_dist(hp_idx[i], hp_idx[j]);
    }
    // NaN peel-off to a complete submatrix (mirrors the main clustering path).
    std::vector<int> valid;
    extract_complete_submatrix_indices(sub, valid);
    const int n = static_cast<int>(valid.size());
    if (n < 8) return;
    Eigen::MatrixXd work(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) work(i, j) = sub(valid[i], valid[j]);
    }

    std::vector<std::string> names(n);
    for (int i = 0; i < n; ++i) names[i] = std::to_string(i);
    HierarchicalClustering clusterer(linkage_method_);
    Tree tree = clusterer.build_tree(work, names);
    if (tree.empty()) return;
    auto kl = TreeCutter::find_optimal_clusters(tree, work, 2, std::min(6, n / 2));
    std::vector<int> labels = kl.second;
    if (static_cast<int>(labels.size()) != n) return;

    // find_optimal_clusters returns best_k (kl.first) but its labels are NOT guaranteed to be the
    // contiguous range [0, best_k): cut_by_num_clusters cuts to a height yielding >= k clusters, and
    // tied merge heights (rife under MAX_DIST sentinel distances) make it overshoot to more clusters
    // than requested. Indexing a best_k-sized array by a raw label is then an out-of-bounds write
    // (heap corruption -> delayed SIGSEGV). Derive the real cluster count and remap to [0, K).
    int max_label = -1;
    for (int l : labels) max_label = std::max(max_label, l);
    if (max_label < 0) return;
    std::vector<int> remap(max_label + 1, -1);
    int k = 0;
    for (int l : labels) {
        if (remap[l] < 0) remap[l] = k++;
    }
    if (k < 2) return;
    for (int& l : labels) l = remap[l];

    // Mean silhouette (1 read = 1 obs): a = mean dist to own cluster, b = min mean dist to another cluster.
    double sil_sum = 0.0;
    int sil_n = 0;
    for (int i = 0; i < n; ++i) {
        double a = 0.0;
        int ca = 0;
        std::vector<double> other_sum(k, 0.0);
        std::vector<int> other_cnt(k, 0);
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (labels[j] == labels[i]) {
                a += work(i, j);
                ++ca;
            } else {
                other_sum[labels[j]] += work(i, j);
                ++other_cnt[labels[j]];
            }
        }
        if (ca == 0) continue;
        a /= ca;
        double b = std::numeric_limits<double>::max();
        for (int c = 0; c < k; ++c) {
            if (c != labels[i] && other_cnt[c] > 0) b = std::min(b, other_sum[c] / other_cnt[c]);
        }
        if (b == std::numeric_limits<double>::max()) continue;
        const double denom = std::max(a, b);
        if (denom > 0) {
            sil_sum += (b - a) / denom;
            ++sil_n;
        }
    }
    out_silhouette = sil_n > 0 ? sil_sum / sil_n : 0.0;

    // Clean multigroup = good silhouette + minimum cluster size. [B] balance relaxed from max(3, n/5) to >=3:
    // the 20% relative balance gate systematically blocked clean BUT MINOR (low-CCF) subclones — within_hp_rederive
    // audit (N=113): 26.5% were sil>=0.5 clean splits rejected only by balance. silhouette>=0.5 still guards cleanness.
    std::vector<int> sizes(k, 0);
    for (int l : labels) ++sizes[l];
    const int min_sz = *std::min_element(sizes.begin(), sizes.end());
    out_min_frac = static_cast<double>(min_sz) / n;  // [D] expose smallest-cluster fraction (minor-subclone indicator)
    if (out_silhouette >= 0.5 && min_sz >= 3) out_ngroups = k;  // [B] relax balance: min_sz>=3 (was max(3, n/5))
}

// [Stage② C'] within-HP a-priori subclone PERMANOVA: does the within-HP methylation DISTANCE separate germline-tag
// vs somatic-carrier-tag reads? Labels are A-PRIORI (read hp_tag), NOT derived from the distance — so this is a VALID
// hypothesis test (no double-dip, unlike the unsupervised tumor_only axis). Statistical "mark/confirm" for within-HP
// subclone structure: reuses StructureTest::run_permanova (location, permutation null) + check_dispersion (PERMDISP).
void RegionProcessor::compute_within_hp_subclone_permanova(const Eigen::MatrixXd& all_dist,
                                                           const std::vector<int>& hp_idx,
                                                           const std::vector<int>& gc_labels, double& out_p,
                                                           double& out_f, bool& out_valid, bool& out_disp_warn) const {
    out_p = 1.0;
    out_f = 0.0;
    out_valid = false;
    out_disp_warn = false;
    const int m = static_cast<int>(hp_idx.size());
    if (m < 8) return;
    Eigen::MatrixXd sub(m, m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) sub(i, j) = all_dist(hp_idx[i], hp_idx[j]);
    }
    // NaN peel-off to a complete submatrix (mirrors the within-HP / tumor-only paths).
    std::vector<int> valid;
    extract_complete_submatrix_indices(sub, valid);
    const int n = static_cast<int>(valid.size());
    if (n < 8) return;
    Eigen::MatrixXd work(n, n);
    std::vector<int> labels(n);
    int n0 = 0, n1 = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) work(i, j) = sub(valid[i], valid[j]);
        labels[i] = gc_labels[valid[i]];  // a-priori germline(0)/carrier(1) — independent of distance => no double-dip
        if (labels[i] == 0) {
            ++n0;
        } else if (labels[i] == 1) {
            ++n1;
        }
    }
    if (n0 < 3 || n1 < 3) return;  // need both germline AND carrier present with enough reads for a valid test
    StructureTestConfig st_config;
    st_config.n_permutations = 999;
    st_config.seed = 42;
    st_config.dispersion_alpha = 0.05;
    StructureTest st(st_config);
    PermanovaResult perm = st.run_permanova(work, labels);
    out_f = perm.pseudo_f;
    out_p = perm.p_value;
    out_valid = perm.valid;
    DispersionResult disp = st.check_dispersion(work, labels);
    out_disp_warn = disp.warning;
}

void RegionProcessor::compute_tumor_only_cluster_structure(const Eigen::MatrixXd& all_dist,
                                                           const std::vector<int>& tumor_idx, int& out_k,
                                                           double& out_silhouette, double& out_permanova_f,
                                                           double& out_permanova_p, bool& out_permanova_valid,
                                                           bool& out_dispersion_warn) const {
    out_k = 1;
    out_silhouette = 0.0;
    out_permanova_f = 0.0;
    out_permanova_p = 1.0;
    out_permanova_valid = false;
    out_dispersion_warn = false;
    const int m = static_cast<int>(tumor_idx.size());
    if (m < 8) return;  // need enough tumor reads for a stable 2-group split + PERMANOVA (Nmin=5)

    // Extract the TUMOR-only sub-distance-matrix from the already-computed all_dist (no recompute).
    Eigen::MatrixXd sub(m, m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) sub(i, j) = all_dist(tumor_idx[i], tumor_idx[j]);
    }
    // NaN peel-off to a complete submatrix (mirrors the main clustering + within-HP path).
    std::vector<int> valid;
    extract_complete_submatrix_indices(sub, valid);
    const int n = static_cast<int>(valid.size());
    if (n < 8) return;
    Eigen::MatrixXd work(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) work(i, j) = sub(valid[i], valid[j]);
    }

    // Unsupervised clustering on tumor reads (UPGMA + silhouette-optimal k), same as the main path.
    std::vector<std::string> names(n);
    for (int i = 0; i < n; ++i) names[i] = std::to_string(i);
    HierarchicalClustering clusterer(linkage_method_);
    Tree tree = clusterer.build_tree(work, names);
    if (tree.empty()) return;
    auto kl = TreeCutter::find_optimal_clusters(tree, work, 2, std::min(6, n / 2));
    std::vector<int> labels = kl.second;
    if (static_cast<int>(labels.size()) != n) return;
    // Remap labels to contiguous [0, K) (same OOB/segfault guard as within-HP; see compute_within_hp_substructure).
    int max_label = -1;
    for (int l : labels) max_label = std::max(max_label, l);
    if (max_label < 0) return;
    std::vector<int> remap(max_label + 1, -1);
    int k = 0;
    for (int l : labels) {
        if (remap[l] < 0) remap[l] = k++;
    }
    if (k < 2) return;
    for (int& l : labels) l = remap[l];
    out_k = k;

    // Mean silhouette (same definition as compute_within_hp_substructure).
    double sil_sum = 0.0;
    int sil_n = 0;
    for (int i = 0; i < n; ++i) {
        double a = 0.0;
        int ca = 0;
        std::vector<double> other_sum(k, 0.0);
        std::vector<int> other_cnt(k, 0);
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (labels[j] == labels[i]) {
                a += work(i, j);
                ++ca;
            } else {
                other_sum[labels[j]] += work(i, j);
                ++other_cnt[labels[j]];
            }
        }
        if (ca == 0) continue;
        a /= ca;
        double b = std::numeric_limits<double>::max();
        for (int c = 0; c < k; ++c) {
            if (c != labels[i] && other_cnt[c] > 0) b = std::min(b, other_sum[c] / other_cnt[c]);
        }
        if (b == std::numeric_limits<double>::max()) continue;
        const double denom = std::max(a, b);
        if (denom > 0) {
            sil_sum += (b - a) / denom;
            ++sil_n;
        }
    }
    out_silhouette = sil_n > 0 ? sil_sum / sil_n : 0.0;

    // Null-gated structure test: silhouette over-detects (random reads still split into compact groups), so
    // gate the tumor-only clusters with PERMANOVA (permutation null, location) + PERMDISP (analytic-F dispersion).
    StructureTestConfig st_config;
    st_config.n_permutations = 999;
    st_config.seed = 42;
    st_config.dispersion_alpha = 0.05;
    StructureTest st(st_config);
    PermanovaResult perm = st.run_permanova(work, labels);
    out_permanova_f = perm.pseudo_f;
    out_permanova_p = perm.p_value;
    out_permanova_valid = perm.valid;
    DispersionResult disp = st.check_dispersion(work, labels);
    out_dispersion_warn = disp.warning;
}

bool RegionProcessor::within_hp_level_clean(const Eigen::MatrixXd& raw, const std::vector<int>& hp_idx,
                                            int min_group) {
    // Per-read mean β over valid CpG (level, not pattern).
    std::vector<double> means;
    for (int idx : hp_idx) {
        double s = 0.0;
        int c = 0;
        for (int j = 0; j < raw.cols(); ++j) {
            const double v = raw(idx, j);
            if (!std::isnan(v)) {
                s += v;
                ++c;
            }
        }
        if (c > 0) means.push_back(s / c);
    }
    const int n = static_cast<int>(means.size());
    if (n < 2 * min_group) return false;
    std::sort(means.begin(), means.end());
    double total = 0.0, mean_all = 0.0;
    for (double m : means) mean_all += m;
    mean_all /= n;
    for (double m : means) total += (m - mean_all) * (m - mean_all);
    if (total <= 0) return false;
    const int min_bal = std::max(min_group, n / 5);  // balanced: both groups >= max(min_group, 20%)
    double best_within = total;
    double best_gap = 0.0;
    for (int s = min_bal; s <= n - min_bal; ++s) {
        double ma = 0.0, mb = 0.0;
        for (int i = 0; i < s; ++i) ma += means[i];
        for (int i = s; i < n; ++i) mb += means[i];
        ma /= s;
        mb /= (n - s);
        double wv = 0.0;
        for (int i = 0; i < s; ++i) wv += (means[i] - ma) * (means[i] - ma);
        for (int i = s; i < n; ++i) wv += (means[i] - mb) * (means[i] - mb);
        if (wv < best_within) {
            best_within = wv;
            best_gap = mb - ma;
        }
    }
    return best_gap > 0.15 && (1.0 - best_within / total) > 0.5;
}

// phylo-v4.1 native labeling (replaces silhouette for the subclone-oriented tumor-only partition).
// Filters tumor + germline-HP-tagged reads, peels to a complete sub-matrix, runs PhyloLabeler with
// modal instability (K=10 null95) + fine (null90) + "other" residual group, writes phylo_groups.tsv.
// Mirrors phylo_v4.py analyze_v4; the binary emits labels so Python only reads + plots.
void RegionProcessor::compute_phylo_groups(const DistanceMatrix& all_dist, const std::vector<ReadInfo>& read_list,
                                           const MethylationMatrix& meth_mat,
                                           const std::string& clustering_dir) const {
    static const std::set<std::string> kHpLabmap = {"1", "HP1", "1-1", "2", "HP2", "2-1"};
    const int MINSZ = 3;
    // 1. Filter: tumor + HP-tagged (matches Python is_tumor + hp in LABMAP).
    std::vector<int> sel;
    for (size_t i = 0; i < read_list.size(); ++i) {
        if (read_list[i].is_tumor && kHpLabmap.count(read_list[i].hp_tag)) sel.push_back(static_cast<int>(i));
    }
    if (static_cast<int>(sel.size()) < 2 * MINSZ) return;
    // 2. Extract tumor sub-distance + peel to a complete sub-matrix.
    const int ns = static_cast<int>(sel.size());
    Eigen::MatrixXd sub(ns, ns);
    for (int i = 0; i < ns; ++i)
        for (int j = 0; j < ns; ++j) sub(i, j) = all_dist.dist_matrix(sel[i], sel[j]);
    std::vector<int> valid;
    extract_complete_submatrix_indices(sub, valid);
    const int n = static_cast<int>(valid.size());
    if (n < 2 * MINSZ) return;
    Eigen::MatrixXd dist(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) dist(i, j) = sub(valid[i], valid[j]);
    Eigen::MatrixXd meth(n, meth_mat.raw_matrix.cols());
    for (int i = 0; i < n; ++i) meth.row(i) = meth_mat.raw_matrix.row(sel[valid[i]]);

    // 3. modal coarse (K=10 null95) — robust verdict + modal_frac + unstable (matches analyze_v4).
    const int K = 10;
    std::vector<PhyloResult> coarse_runs;
    std::map<int, int> ng_count;
    for (int k = 0; k < K; ++k) {
        PhyloConfig cfg;
        cfg.null_pct = 95.0;
        cfg.seed = 20260622ull + static_cast<uint64_t>(k) * 101ull;
        PhyloLabeler lab(cfg);
        PhyloResult r = lab.label(dist, meth);
        coarse_runs.push_back(r);
        ng_count[r.n_groups]++;
    }
    int modal_ng = 0, modal_cnt = 0;
    for (auto& kv : ng_count)
        if (kv.second > modal_cnt) {
            modal_cnt = kv.second;
            modal_ng = kv.first;
        }
    double modal_frac = static_cast<double>(modal_cnt) / K;
    bool unstable = modal_frac < 0.7;
    int ng_min = 1 << 30, ng_max = 0;
    for (auto& kv : ng_count) {
        ng_min = std::min(ng_min, kv.first);
        ng_max = std::max(ng_max, kv.first);
    }
    // representative coarse labeling = first run that produced modal_ng
    PhyloResult coarse = coarse_runs.front();
    for (auto& r : coarse_runs)
        if (r.n_groups == modal_ng) {
            coarse = r;
            break;
        }
    // 4. fine (null90 candidate, low confidence).
    PhyloConfig fcfg;
    fcfg.null_pct = 90.0;
    fcfg.seed = 20260622ull;
    PhyloResult fine = PhyloLabeler(fcfg).label(dist, meth);

    // 5. "other" residual group: coarse outliers >= MINSZ -> "other".
    std::vector<std::string> clab = coarse.labels;
    int n_out = static_cast<int>(std::count(clab.begin(), clab.end(), std::string("outlier")));
    if (n_out >= MINSZ)
        for (auto& l : clab)
            if (l == "outlier") l = "other";
    int n_other = static_cast<int>(std::count(clab.begin(), clab.end(), std::string("other")));
    bool hidden_het = n > 0 && (static_cast<double>(n_other) / n > 0.30);

    // 6. write phylo_groups.tsv (read_id/coarse/fine/is_other/is_outlier) + summary.
    std::ofstream tsv(clustering_dir + "/phylo_groups.tsv");
    if (tsv.is_open()) {
        tsv << "read_id\tread_name\thp\talt_support\tcoarse_label\tfine_label\tis_other\tis_outlier\n";
        for (int i = 0; i < n; ++i) {
            const ReadInfo& r = read_list[sel[valid[i]]];
            const std::string& cl = clab[i];
            const std::string& fl = fine.valid ? fine.labels[i] : std::string("NA");
            std::string al = (r.alt_support == AltSupport::REF)   ? "REF"
                             : (r.alt_support == AltSupport::ALT) ? "ALT"
                                                                  : "UNKNOWN";
            tsv << r.read_id << "\t" << r.read_name << "\t" << r.hp_tag << "\t" << al << "\t" << cl << "\t" << fl
                << "\t" << (cl == "other" ? 1 : 0) << "\t" << (cl == "outlier" ? 1 : 0) << "\n";
        }
        tsv.close();
    }
    std::ofstream js(clustering_dir + "/phylo_groups_summary.json");
    if (js.is_open()) {
        js << "{\"n\":" << n << ",\"coarse_ng\":" << modal_ng << ",\"modal_frac\":" << modal_frac
           << ",\"fine_ng\":" << (fine.valid ? fine.n_groups : 0) << ",\"n_other\":" << n_other
           << ",\"unstable\":" << (unstable ? "true" : "false") << ",\"ng_range\":[" << ng_min << "," << ng_max
           << "],\"hidden_het\":" << (hidden_het ? "true" : "false") << ",\"seed\":20260622,\"rnull\":40}\n";
        js.close();
    }
}

void RegionProcessor::perform_clustering_and_significance(const DistanceMatrix& all_dist,
                                                           const std::vector<ReadInfo>& read_list,
                                                           const MethylationMatrix& meth_mat,
                                                           const std::string& clustering_dir,
                                                           const std::string& chr_name, const SomaticSnv& snv,
                                                           int region_id, RegionResult& result) {
    std::filesystem::create_directories(clustering_dir);

    // phylo-v4.1 native cluster labeling (binary emits labels; Python only reads+plots).
    compute_phylo_groups(all_dist, read_list, meth_mat, clustering_dir);

    // [P1a] HP-AUC: does the read-read distance recover the germline-HP label?
    // Ground truth = HP tag; computed on the full matrix (NaN pairs skipped under SKIP),
    // per subset (normal-only cleanest germline truth; tumor often single-HP -> -1).
    {
        std::vector<int> hp_fam(read_list.size(), -1);
        for (size_t i = 0; i < read_list.size(); ++i) {
            const std::string& h = read_list[i].hp_tag;
            if (!h.empty() && h[0] == '1')
                hp_fam[i] = 0;
            else if (!h.empty() && h[0] == '2')
                hp_fam[i] = 1;
        }
        std::vector<int> idx_n, idx_t, idx_a;
        for (size_t i = 0; i < read_list.size(); ++i) {
            if (hp_fam[i] < 0) continue;
            idx_a.push_back(static_cast<int>(i));
            (read_list[i].is_tumor ? idx_t : idx_n).push_back(static_cast<int>(i));
        }
        result.hp_auc_normal = compute_hp_auc(all_dist.dist_matrix, hp_fam, idx_n);
        result.hp_auc_tumor = compute_hp_auc(all_dist.dist_matrix, hp_fam, idx_t);
        result.hp_auc_all = compute_hp_auc(all_dist.dist_matrix, hp_fam, idx_a);
    }

    // [SKIP-fix] Filter reads to a complete (NaN-free) sub-matrix before clustering.
    // Hierarchical clustering / tree-cutting cannot consume NaN distances: a NaN pair is
    // never selected as the UPGMA minimum, producing a degenerate tree whose traversal
    // stack-overflows (SIGSEGV). This mirrors the downstream PERMANOVA filter, so
    // clustering and significance operate on the SAME valid read subset. For MAX_DIST
    // (no NaN) the subset == all reads => identical to prior behavior (regression-safe).
    // See fix/clustering-nan-skip-segfault (2026-06-14).
    const Eigen::MatrixXd& full_dist = all_dist.dist_matrix;
    const int n_total = static_cast<int>(read_list.size());
    std::vector<int> valid_idx;
    extract_complete_submatrix_indices(full_dist, valid_idx);
    const int n_valid = static_cast<int>(valid_idx.size());
    if (n_valid < 2) return;  // not enough complete-overlap reads to cluster
    const bool is_subset = (n_valid != n_total);
    result.num_reads_valid = n_valid;  // [gap1] reads used for clustering after NaN-pair drop (== num_reads under MAX_DIST)

    // Working distance matrix + read names in valid-read space (zero-copy when no drop).
    std::vector<std::string> read_names;
    read_names.reserve(n_valid);
    Eigen::MatrixXd sub_dist;
    if (is_subset) {
        sub_dist.resize(n_valid, n_valid);
        for (int i = 0; i < n_valid; ++i) {
            read_names.push_back(read_list[valid_idx[i]].read_name);
            for (int j = 0; j < n_valid; ++j) {
                sub_dist(i, j) = full_dist(valid_idx[i], valid_idx[j]);
            }
        }
    } else {
        for (const auto& r : read_list) read_names.push_back(r.read_name);
    }
    const Eigen::MatrixXd& work_dist = is_subset ? sub_dist : full_dist;

    // Build hierarchical clustering tree
    HierarchicalClustering clusterer(linkage_method_);
    Tree tree = clusterer.build_tree(work_dist, read_names);

    if (tree.empty()) return;

    if (output_tree_files_) {
        TreeWriter tree_writer;
        tree_writer.write_newick(tree, clustering_dir + "/tree.nwk");

        if (output_linkage_matrix_) {
            tree_writer.write_linkage_matrix(tree, clustering_dir + "/linkage_matrix.csv");
        }

        // Write leaf order
        std::ofstream order_file(clustering_dir + "/leaf_order.txt");
        if (order_file.is_open()) {
            auto leaves = tree.get_leaves();
            for (const auto& leaf : leaves) {
                order_file << leaf->label << "\n";
            }
            order_file.close();
        }
    }

    if (log_level_ >= LogLevel::LOG_TRACE) {
        LOG_TRACE("  Clustering tree: " + std::to_string(tree.num_leaves()) +
                  " leaves, method=" + HierarchicalClustering::method_to_string(linkage_method_));
    }

    // Significance Analysis
    try {
        auto [best_k, cluster_labels] =
            TreeCutter::find_optimal_clusters(tree, work_dist, 2, std::min(6, n_valid / 2));

        if (best_k < 2 || static_cast<int>(cluster_labels.size()) != n_valid) return;

        // Build FullLabel vector (valid-read space, aligned with cluster_labels / work_dist)
        std::vector<FullLabel> full_labels;
        full_labels.reserve(n_valid);
        for (int k = 0; k < n_valid; ++k) {
            const auto& read = read_list[valid_idx[k]];
            FullLabel label;
            label.allele = read.alt_support;
            label.hp_tag = read.hp_tag;
            label.strand = read.strand;
            label.is_tumor = read.is_tumor;
            full_labels.push_back(label);
        }

        // Check for unbalanced clusters
        std::map<int, std::vector<int>> cluster_read_indices;
        for (size_t i = 0; i < cluster_labels.size(); ++i) {
            cluster_read_indices[cluster_labels[i]].push_back(i);
        }

        int min_cluster_size = std::max(3, n_valid / 20);
        bool unbalanced = false;
        for (const auto& [cid, indices] : cluster_read_indices) {
            if (static_cast<int>(indices.size()) < min_cluster_size) {
                unbalanced = true;
                break;
            }
        }

        // Try increasing k if unbalanced
        if (unbalanced && best_k < 5) {
            int next_k = best_k + 1;
            std::vector<int> next_labels = TreeCutter::cut_by_num_clusters(tree, next_k);

            std::map<int, int> next_counts;
            for (int cl : next_labels) next_counts[cl]++;

            int substantial_clusters = 0;
            for (const auto& [cl, count] : next_counts) {
                if (count >= min_cluster_size) substantial_clusters++;
            }

            if (substantial_clusters >= 2) {
                best_k = next_k;
                cluster_labels = next_labels;
            }
        }

        // Configure and run significance analysis
        SignificanceConfig sig_config;
        sig_config.enable_global_test = true;
        sig_config.enable_local_test = true;
        sig_config.enable_permanova = true;
        sig_config.structure_config.n_permutations = 999;  // gap#5: unify with LabelTest 999 (p-floor 0.001)
        sig_config.seed = 42;  // Fixed top-level seed propagates to all components via init_analyzers()
        sig_config.enable_dispersion = true;  // PERMDISP: flag dispersion-driven PERMANOVA hits (analytic F p)
        sig_config.enable_bootstrap = false;
        sig_config.run_id = vcf_filename_;
        sig_config.vcf_id = vcf_filename_;

        std::string anchor_key = chr_name + "_" + std::to_string(snv.pos);
        SignificanceAnalyzer analyzer(sig_config);

        // Methylation matrix in valid-read space (rows aligned with cluster_labels / work_dist)
        const Eigen::MatrixXd& full_meth = meth_mat.raw_matrix;
        Eigen::MatrixXd sub_meth;
        if (is_subset) {
            sub_meth.resize(n_valid, full_meth.cols());
            for (int i = 0; i < n_valid; ++i) {
                sub_meth.row(i) = full_meth.row(valid_idx[i]);
            }
        }
        const Eigen::MatrixXd& work_meth = is_subset ? sub_meth : full_meth;

        SignificanceResult sig_result =
            analyzer.analyze(cluster_labels, full_labels, work_dist, work_meth, region_id, anchor_key);

        // Store results
        result.significance_computed = true;
        double p_alt = sig_result.global_alt.valid ? sig_result.global_alt.fisher_ffh.p_value : 1.0;
        double p_hp = sig_result.global_hp.valid ? sig_result.global_hp.fisher_ffh.p_value : 1.0;
        double p_hp_family = sig_result.global_hp_family.valid ? sig_result.global_hp_family.fisher_ffh.p_value : 1.0;
        result.global_p_value = std::min({p_alt, p_hp, p_hp_family});

        double v_alt = sig_result.global_alt.cramers_v_reliable ? sig_result.global_alt.cramers_v : 0.0;
        double v_hp = sig_result.global_hp.cramers_v_reliable ? sig_result.global_hp.cramers_v : 0.0;
        double v_hp_family =
            sig_result.global_hp_family.cramers_v_reliable ? sig_result.global_hp_family.cramers_v : 0.0;
        result.cramers_v = std::max({v_alt, v_hp, v_hp_family});

        // Multi-layer HP results
        result.global_hp_family_p = p_hp_family;
        result.global_hp_family_v = v_hp_family;
        result.global_hp_fine_p = sig_result.global_hp_fine.valid ? sig_result.global_hp_fine.fisher_ffh.p_value : 1.0;
        result.global_hp_fine_v =
            sig_result.global_hp_fine.cramers_v_reliable ? sig_result.global_hp_fine.cramers_v : 0.0;
        result.global_hp_fine_n_groups = sig_result.global_hp_fine.n_valid_groups;
        result.local_best_p_value = sig_result.local_result.best_p_value;
        result.local_best_cluster = sig_result.local_result.best_cluster_id;
        result.heuristic_score = sig_result.heuristic_score;
        result.passed_gating = sig_result.passed_gate;
        result.cluster_permanova_f = sig_result.permanova.pseudo_f;
        result.cluster_permanova_p = sig_result.permanova.p_value;
        result.cluster_permanova_valid = sig_result.permanova.valid;
        result.cluster_dispersion_p = sig_result.dispersion.anova_p;
        result.cluster_dispersion_warning = sig_result.dispersion_warning;

        result.label_test_computed = true;

        // Multi-Stage HP Verification results (NEW)
        const auto& hp_ms = sig_result.hp_multistage;

        // Stage 1: HP Family Merged Test
        result.hp_merged_delta = hp_ms.merged.delta;
        result.hp_merged_p = hp_ms.merged.p_value;
        result.hp_merged_sig = hp_ms.merged.significant;
        result.hp_merged_n_hp1 = hp_ms.n_hp1_family;
        result.hp_merged_n_hp2 = hp_ms.n_hp2_family;

        // Stage 2: HP Fine-Grained Test
        result.hp_fine_f = hp_ms.fine_f;
        result.hp_fine_p = hp_ms.fine_p;
        result.hp_fine_sig = hp_ms.fine_sig;
        result.hp_fine_n_groups = hp_ms.fine_n_groups;
        // Fine group counts (HP1=0, HP1-1=1, HP2=2, HP2-1=3)
        if (hp_ms.fine_group_counts.size() >= 4) {
            result.hp_fine_n_hp1 = hp_ms.fine_group_counts[0];
            result.hp_fine_n_hp1s = hp_ms.fine_group_counts[1];
            result.hp_fine_n_hp2 = hp_ms.fine_group_counts[2];
            result.hp_fine_n_hp2s = hp_ms.fine_group_counts[3];
        }
        // Fine pairwise mean distances
        result.hp_fine_d_hp1_hp1s = hp_ms.fine_pairwise.d_hp1_hp1s;
        result.hp_fine_d_hp1_hp2 = hp_ms.fine_pairwise.d_hp1_hp2;
        result.hp_fine_d_hp1_hp2s = hp_ms.fine_pairwise.d_hp1_hp2s;
        result.hp_fine_d_hp1s_hp2 = hp_ms.fine_pairwise.d_hp1s_hp2;
        result.hp_fine_d_hp1s_hp2s = hp_ms.fine_pairwise.d_hp1s_hp2s;
        result.hp_fine_d_hp2_hp2s = hp_ms.fine_pairwise.d_hp2_hp2s;

        // Stage 3: Allele Test
        result.allele_delta = sig_result.label_allele.delta;
        result.allele_p = sig_result.label_allele.p_value;
        result.allele_sig = sig_result.label_allele.significant;
        result.label_hp_permanova_f = sig_result.label_hp_permanova.pseudo_f;
        result.label_hp_permanova_p = sig_result.label_hp_permanova.p_value;
        result.label_hp_permanova_valid = sig_result.label_hp_permanova.valid;
        result.label_hp_dispersion_p = sig_result.label_hp_dispersion.anova_p;
        result.label_hp_dispersion_warning = sig_result.label_hp_dispersion_warning;
        result.label_allele_permanova_f = sig_result.label_allele_permanova.pseudo_f;
        result.label_allele_permanova_p = sig_result.label_allele_permanova.p_value;
        result.label_allele_permanova_valid = sig_result.label_allele_permanova.valid;
        result.label_allele_dispersion_p = sig_result.label_allele_dispersion.anova_p;
        result.label_allele_dispersion_warning = sig_result.label_allele_dispersion_warning;

        // Stage 4: Unassigned Affinity Test
        result.unassigned_affinity = hp_ms.unassigned.affinity_score;
        result.unassigned_affinity_p = hp_ms.unassigned.affinity_p;
        result.unassigned_dir = hp_ms.unassigned.affinity_direction;
        result.unassigned_n_hp3 = hp_ms.unassigned.n_hp3;
        result.unassigned_n_hp0 = hp_ms.unassigned.n_hp0;

        // Sample ASM (Tumor vs Normal) - Label-First
        result.sample_asm_delta = sig_result.label_sample.delta;
        result.sample_asm_p = sig_result.label_sample.p_value;
        result.sample_asm_sig = sig_result.label_sample.significant;
        result.sample_asm_n_tumor = sig_result.label_sample.n_group_a;
        result.sample_asm_n_normal = sig_result.label_sample.n_group_b;

        // Backward compatibility: use Stage 1 result as primary label_delta
        if (hp_ms.merged.valid) {
            result.label_delta = (sig_result.dominant_label_dimension == "hp") ? hp_ms.merged.delta
                                                                               : sig_result.label_allele.delta;
            result.label_p_value = (sig_result.dominant_label_dimension == "hp") ? hp_ms.merged.p_value
                                                                                  : sig_result.label_allele.p_value;
        } else if (sig_result.label_allele.valid) {
            result.label_delta = sig_result.label_allele.delta;
            result.label_p_value = sig_result.label_allele.p_value;
        }
        result.label_significant = sig_result.label_significant;
        result.dominant_label = sig_result.dominant_label_dimension;
        result.verification_class = sig_result.verification_class;
        result.n_clusters = best_k;

        // Multi-Layer Validation: Compute quality metrics (NEW)
        result.hp_ratio = compute_hp_ratio(result.hp_merged_n_hp1, result.hp_merged_n_hp2);
        result.potential_loh = is_potential_loh(result.hp_ratio);
        result.coverage_multiple = compute_coverage_multiple(result.num_reads);
        result.diploid_coverage_used = 75.0;  // Pass 1 default; overwritten by Pass 2
        result.coverage_category = determine_coverage_category(result.coverage_multiple);
        result.loh_subtype = determine_loh_subtype(result.potential_loh, result.verification_class);
        auto qs_weights = normal_bam_path_.empty() ? get_tumor_only_weights() : get_paired_weights();
        result.quality_score = compute_quality_score(
            result.num_reads, result.num_cpgs, result.coverage_multiple,
            result.potential_loh, result.hp_merged_sig, result.allele_sig, result.cramers_v,
            qs_weights);
        result.quality_tier = determine_quality_tier(result.quality_score);

        // LOH BED annotation (Phase C)
        if (loh_annotator_.loaded()) {
            result.loh_bed_overlap = loh_annotator_.overlaps(chr_name, snv.pos);
            LohSource loh_src = loh_annotator_.classify(chr_name, snv.pos, result.potential_loh);
            result.loh_source = loh_source_to_string(loh_src);
            if (result.loh_bed_overlap) {
                result.loh_bed_annotation = loh_annotator_.get_annotation(chr_name, snv.pos);
            }
        }

        // Write significance results to JSON
        std::ofstream sig_file(clustering_dir + "/significance.json");
        if (sig_file.is_open()) {
            sig_file << "{\n";
            sig_file << "  \"anchor_key\": \"" << anchor_key << "\",\n";
            sig_file << "  \"num_reads\": " << read_list.size() << ",\n";
            sig_file << "  \"optimal_k\": " << best_k << ",\n";
            sig_file << "  \"passed_gating\": " << (sig_result.passed_gate ? "true" : "false") << ",\n";
            sig_file << "  \"pairwise_distance\": {\n";
            sig_file << "    \"mean\": " << std::fixed << std::setprecision(4) << result.pairwise_mean_distance
                     << ",\n";
            sig_file << "    \"median\": " << std::fixed << std::setprecision(4) << result.pairwise_median_distance
                     << "\n";
            sig_file << "  },\n";
            sig_file << "  \"global_alt\": {\n";
            sig_file << "    \"p_value\": " << std::scientific << std::setprecision(6)
                     << sig_result.global_alt.fisher_ffh.p_value << ",\n";
            sig_file << "    \"cramers_v\": " << std::fixed << std::setprecision(4) << sig_result.global_alt.cramers_v
                     << "\n";
            sig_file << "  },\n";
            sig_file << "  \"global_hp\": {\n";
            sig_file << "    \"p_value\": " << std::scientific << std::setprecision(6)
                     << sig_result.global_hp.fisher_ffh.p_value << ",\n";
            sig_file << "    \"cramers_v\": " << std::fixed << std::setprecision(4) << sig_result.global_hp.cramers_v
                     << "\n";
            sig_file << "  },\n";
            sig_file << "  \"global_hp_family\": {\n";
            sig_file << "    \"p_value\": " << std::scientific << std::setprecision(6)
                     << sig_result.global_hp_family.fisher_ffh.p_value << ",\n";
            sig_file << "    \"cramers_v\": " << std::fixed << std::setprecision(4)
                     << sig_result.global_hp_family.cramers_v << "\n";
            sig_file << "  },\n";
            sig_file << "  \"global_hp_fine\": {\n";
            sig_file << "    \"valid\": " << (sig_result.global_hp_fine.valid ? "true" : "false") << ",\n";
            sig_file << "    \"p_value\": " << std::scientific << std::setprecision(6)
                     << sig_result.global_hp_fine.fisher_ffh.p_value << ",\n";
            sig_file << "    \"cramers_v\": " << std::fixed << std::setprecision(4)
                     << sig_result.global_hp_fine.cramers_v << "\n";
            sig_file << "  },\n";
            sig_file << "  \"cluster_structure\": {\n";
            sig_file << "    \"permanova_valid\": " << (sig_result.permanova.valid ? "true" : "false") << ",\n";
            sig_file << "    \"permanova_f\": " << std::fixed << std::setprecision(4) << sig_result.permanova.pseudo_f
                     << ",\n";
            sig_file << "    \"permanova_p\": " << std::scientific << std::setprecision(6)
                     << sig_result.permanova.p_value << ",\n";
            sig_file << "    \"dispersion_p\": " << std::scientific << std::setprecision(6)
                     << sig_result.dispersion.anova_p << ",\n";
            sig_file << "    \"dispersion_warning\": "
                     << (sig_result.dispersion_warning ? "true" : "false") << "\n";
            sig_file << "  },\n";
            sig_file << "  \"label_structure\": {\n";
            sig_file << "    \"hp_permanova_valid\": "
                     << (sig_result.label_hp_permanova.valid ? "true" : "false") << ",\n";
            sig_file << "    \"hp_permanova_f\": " << std::fixed << std::setprecision(4)
                     << sig_result.label_hp_permanova.pseudo_f << ",\n";
            sig_file << "    \"hp_permanova_p\": " << std::scientific << std::setprecision(6)
                     << sig_result.label_hp_permanova.p_value << ",\n";
            sig_file << "    \"hp_dispersion_p\": " << std::scientific << std::setprecision(6)
                     << sig_result.label_hp_dispersion.anova_p << ",\n";
            sig_file << "    \"hp_dispersion_warning\": "
                     << (sig_result.label_hp_dispersion_warning ? "true" : "false") << ",\n";
            sig_file << "    \"allele_permanova_valid\": "
                     << (sig_result.label_allele_permanova.valid ? "true" : "false") << ",\n";
            sig_file << "    \"allele_permanova_f\": " << std::fixed << std::setprecision(4)
                     << sig_result.label_allele_permanova.pseudo_f << ",\n";
            sig_file << "    \"allele_permanova_p\": " << std::scientific << std::setprecision(6)
                     << sig_result.label_allele_permanova.p_value << ",\n";
            sig_file << "    \"allele_dispersion_p\": " << std::scientific << std::setprecision(6)
                     << sig_result.label_allele_dispersion.anova_p << ",\n";
            sig_file << "    \"allele_dispersion_warning\": "
                     << (sig_result.label_allele_dispersion_warning ? "true" : "false") << "\n";
            sig_file << "  },\n";
            sig_file << "  \"heuristic_score\": " << std::fixed << std::setprecision(4) << sig_result.heuristic_score
                     << "\n";
            sig_file << "}\n";
            sig_file.close();
        }

    } catch (const std::exception& e) {
        if (log_level_ >= LogLevel::LOG_DEBUG) {
            LOG_DEBUG("Significance analysis skipped for region " + std::to_string(region_id) + ": " + std::string(e.what()));
        }
    }
}

void RegionProcessor::write_strand_specific_trees(const DistanceMatrix& forward_dist,
                                                   const DistanceMatrix& reverse_dist,
                                                   const std::vector<ReadInfo>& read_list,
                                                   const std::string& clustering_dir) {
    HierarchicalClustering clusterer(linkage_method_);
    TreeWriter tree_writer;

    // Forward strand tree
    if (forward_dist.size() >= 2) {
        std::vector<std::string> fwd_names;
        for (const auto& r : read_list) {
            if (r.strand == Strand::FORWARD) {
                fwd_names.push_back(r.read_name);
            }
        }
        if (fwd_names.size() >= 2) {
            Tree fwd_tree = clusterer.build_tree(forward_dist, fwd_names);
            if (!fwd_tree.empty()) {
                tree_writer.write_newick(fwd_tree, clustering_dir + "/tree_forward.nwk");
            }
        }
    }

    // Reverse strand tree
    if (reverse_dist.size() >= 2) {
        std::vector<std::string> rev_names;
        for (const auto& r : read_list) {
            if (r.strand == Strand::REVERSE) {
                rev_names.push_back(r.read_name);
            }
        }
        if (rev_names.size() >= 2) {
            Tree rev_tree = clusterer.build_tree(reverse_dist, rev_names);
            if (!rev_tree.empty()) {
                tree_writer.write_newick(rev_tree, clustering_dir + "/tree_reverse.nwk");
            }
        }
    }
}

}  // namespace InterSubMod
