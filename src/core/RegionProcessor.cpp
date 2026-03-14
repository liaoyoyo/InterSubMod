#include "core/RegionProcessor.hpp"

#include <omp.h>
#include <sys/stat.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
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
#include "core/MethylationMatrix.hpp"
#include "core/SignificanceAnalyzer.hpp"
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

/**
 * @brief Compute coverage multiple relative to expected 75x
 */
double compute_coverage_multiple(int num_reads, double expected_coverage = 75.0) {
    return static_cast<double>(num_reads) / expected_coverage;
}

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
 * @brief Compute comprehensive quality score
 *
 * Scoring components:
 * - Base score: 100
 * - Read count penalty: -20 if < 30, -10 if < 50
 * - CpG count penalty: -15 if < 20, -10 if < 30
 * - Coverage penalty: -30 if CNV_Loss, -15 if Low, -20 if High_Copy
 * - LOH penalty: -25 if potential LOH
 * - Verification bonus: +10 if HP and Allele tests agree, +15 if V >= 0.3
 */
float compute_quality_score(int num_reads, int num_cpgs, double coverage_multiple,
                            bool potential_loh, bool hp_sig, bool allele_sig, double cramers_v) {
    float score = 100.0f;

    // Read count penalty
    if (num_reads < 30) {
        score -= 20.0f;
    } else if (num_reads < 50) {
        score -= 10.0f;
    }

    // CpG count penalty
    if (num_cpgs < 20) {
        score -= 15.0f;
    } else if (num_cpgs < 30) {
        score -= 10.0f;
    }

    // Coverage anomaly penalty
    if (coverage_multiple < 0.5) {
        score -= 30.0f;  // Likely CNV loss
    } else if (coverage_multiple < 0.8) {
        score -= 15.0f;  // Low coverage
    } else if (coverage_multiple > 2.0) {
        score -= 20.0f;  // Very high coverage (potential artifact)
    }

    // LOH penalty
    if (potential_loh) {
        score -= 25.0f;
    }

    // Verification consistency bonus
    if (hp_sig && allele_sig) {
        score += 10.0f;  // Both tests significant
    }

    // Effect size bonus
    if (cramers_v >= 0.3) {
        score += 15.0f;  // Strong effect size
    } else if (cramers_v >= 0.2) {
        score += 5.0f;   // Moderate effect size
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
      clustering_min_reads_(config.clustering_min_reads) {
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

    // Set filter configuration
    filter_config_.min_mapq = config.min_mapq;
    filter_config_.min_read_length = config.min_read_length;
    filter_config_.min_base_quality = config.min_base_quality;
    filter_config_.require_mm_ml = true;

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
            results[i] = process_single_region(snv, i, tumor_reader, ref_reader);

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

    // Write significance summary and print stats
    write_significance_summary(results);
    print_summary(results);

    return results;
}

RegionResult RegionProcessor::process_single_region(const SomaticSnv& snv, int region_id, BamReader& bam_reader,
                                                    FastaReader& fasta_reader) {
    RegionResult result;
    result.region_id = region_id;
    result.snv_id = snv.snv_id;

    auto t_start = std::chrono::high_resolution_clock::now();

    try {
        MatrixBuilder matrix_builder;

        // Define initial region window (centered on SNV)
        int32_t region_start = static_cast<int32_t>(snv.pos) - static_cast<int32_t>(window_size_);
        int32_t region_end = static_cast<int32_t>(snv.pos) + static_cast<int32_t>(window_size_);

        // Get chromosome name and length
        std::string chr_name = chrom_index_.get_name(snv.chr_id);
        int32_t chr_length = fasta_reader.get_chr_length(chr_name);

        // Clamp coordinates
        if (region_start < 1) region_start = 1;
        if (chr_length > 0 && region_end > chr_length) region_end = chr_length;

        // Fetch reads from BAM
        auto reads = bam_reader.fetch_reads(chr_name, region_start, region_end);

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

        // Process reads using helper method
        std::vector<FilteredReadInfo> filtered_reads;
        process_reads(reads, snv, ref_seq, region_start, matrix_builder, filtered_reads, result);

        // Finalize matrix construction
        matrix_builder.finalize();
        result.num_reads = matrix_builder.num_reads();
        result.num_cpgs = matrix_builder.num_cpgs();
        result.num_filtered_reads = filtered_reads.size();

        // Write output to disk
        RegionWriter writer(output_dir_, debug_output_dir_, true, vcf_filename_);
        writer.write_region(snv, chr_name, region_id, region_start, region_end, matrix_builder.get_reads(),
                            matrix_builder.get_cpg_positions(), matrix_builder.get_matrix(), 0.0, 0.0);

        // Write filtered reads in debug mode
        if (output_filtered_reads_ && !filtered_reads.empty()) {
            std::string region_dir = writer.get_region_dir(chr_name, snv.pos, region_start, region_end);
            writer.write_filtered_reads(region_dir, chr_name, filtered_reads);
        }

        // Distance matrix, clustering, and significance analysis
        if (compute_distance_matrix_ && result.num_reads >= 2 && result.num_cpgs >= 1) {
            MethylationMatrix meth_mat = build_methylation_matrix(matrix_builder, region_id);
            const auto& read_list = matrix_builder.get_reads();
            std::string region_dir = writer.get_region_dir(chr_name, snv.pos, region_start, region_end);

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
        }

        // Cleanup BAM records
        for (auto* r : reads) {
            bam_destroy1(r);
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
    csv_file << "RegionID,Chr,Pos,Ref,Alt,NumReads,NumCpGs,GlobalP,CramersV,HeuristicScore,PassedGating,"
                "PairwiseMeanDist,PairwiseMedianDist,"
                "ClusterPermanovaF,ClusterPermanovaP,ClusterPermanovaValid,ClusterDispersionP,ClusterDispersionWarn,"
                // Stage 1: HP Merged
                "HPMergedDelta,HPMergedP,HPMergedSig,HP1FamilyN,HP2FamilyN,"
                // Stage 2: HP Fine-Grained
                "HPFineF,HPFineP,HPFineSig,HPFineNGroups,"
                // Stage 3: Allele
                "AlleleDelta,AlleleP,AlleleSig,"
                "LabelHPPermanovaF,LabelHPPermanovaP,LabelHPPermanovaValid,LabelHPDispersionP,LabelHPDispersionWarn,"
                "LabelAllelePermanovaF,LabelAllelePermanovaP,LabelAllelePermanovaValid,LabelAlleleDispersionP,LabelAlleleDispersionWarn,"
                // Stage 4: Unassigned Affinity
                "UnassignedAffinity,UnassignedAffinityP,UnassignedDir,NHP3,NHP0,"
                // Summary
                "DominantLabel,Stability,VerificationClass,"
                // Multi-Layer Validation Quality Assessment (NEW)
                "HP_Ratio,Potential_LOH,Coverage_Multiple,Coverage_Category,LOH_Subtype,"
                "Quality_Score,Quality_Tier,"
                // Original columns
                "LocalBestCluster,LocalBestP,Significant,SuggestFilter\n";

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
        // Definition: Passed gating AND p-value <= 0.05 AND V >= 0.1 AND depth >= 20 (Round 7)
        // Optimized based on depth analysis: Depth 20-30 has high TP/FP ratio (2.31)
        bool is_significant = r.passed_gating && 
                              (r.global_p_value <= 0.05) && 
                              (r.cramers_v >= 0.1) &&
                              (r.num_reads >= 20);
        
        // Suggest filtering sites with LabelDelta > 0.3 (F1 optimization)
        bool suggest_filter = (r.label_delta > 0.3);

        if (is_significant) {
            total_significant++;
            chr_stats[chr_name].significant++;
        }
        chr_stats[chr_name].total++;

        csv_file << r.region_id << "," << chr_name << "," << snv.pos << "," << snv.ref_base << "," << snv.alt_base
                 << "," << r.num_reads << "," << r.num_cpgs << "," << std::scientific << std::setprecision(6)
                 << r.global_p_value << "," << std::fixed << std::setprecision(4) << r.cramers_v << "," << std::fixed
                 << std::setprecision(4) << r.heuristic_score << "," << (r.passed_gating ? "true" : "false") << ","
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
                 << r.verification_class << ","
                 // Multi-Layer Validation Quality Assessment (NEW)
                 << std::fixed << std::setprecision(4) << r.hp_ratio << ","
                 << (r.potential_loh ? "true" : "false") << ","
                 << std::fixed << std::setprecision(3) << r.coverage_multiple << ","
                 << r.coverage_category << ","
                 << r.loh_subtype << ","
                 << std::fixed << std::setprecision(1) << r.quality_score << ","
                 << r.quality_tier << ","
                 // Original local test columns
                 << r.local_best_cluster << "," << std::scientific
                 << std::setprecision(6) << r.local_best_p_value << ","
                 << (is_significant ? "true" : "false") << ","
                 // SuggestFilter column for F1 optimization
                 << (suggest_filter ? "true" : "false")
                 << "\n";
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

// ============================================================================
// Refactored Helper Methods
// ============================================================================

void RegionProcessor::process_reads(const std::vector<bam1_t*>& reads, const SomaticSnv& snv,
                                    const std::string& ref_seq, int32_t region_start,
                                    MatrixBuilder& matrix_builder, std::vector<FilteredReadInfo>& filtered_reads,
                                    RegionResult& result) {
    ReadParser read_parser(filter_config_);
    MethylationParser methyl_parser;

    int read_count = 0;
    std::unordered_set<std::string> processed_read_names;

    for (auto* b : reads) {
        // Filter and Parse Read Info
        auto [keep, filter_reason] = read_parser.should_keep_with_reason(b);

        if (!keep && !no_filter_output_) {
            if (output_filtered_reads_) {
                FilteredReadInfo filtered = read_parser.create_filtered_info(b, true, filter_reason);
                filtered_reads.push_back(filtered);
            }
            continue;
        }

        // Even in no-filter mode, always skip structurally invalid reads (secondary,
        // supplementary, unmapped, duplicate) because they have no valid sequence
        // for methylation parsing and cause undefined behavior in MethylationParser.
        if (no_filter_output_) {
            uint16_t flag = b->core.flag;
            if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FDUP)) {
                continue;
            }
        }

        ReadInfo info = read_parser.parse(b, read_count, true, snv, ref_seq, region_start);

        // Skip reads that don't cover the SNV
        if (info.alt_support == AltSupport::UNKNOWN && !no_filter_output_) {
            if (output_filtered_reads_) {
                FilterReason reason = FilterReason::SNV_NOT_COVERED;
                FilteredReadInfo filtered = read_parser.create_filtered_info(b, true, reason);
                filtered_reads.push_back(filtered);
            }
            continue;
        }

        // Duplicate read name check
        if (processed_read_names.find(info.read_name) != processed_read_names.end()) {
            continue;
        }
        processed_read_names.insert(info.read_name);

        // Parse Methylation
        auto methyl_calls = methyl_parser.parse_read(b, ref_seq, region_start);

        // Add to Matrix Builder
        matrix_builder.add_read(info, methyl_calls);
        read_count++;

        // Count strand
        if (info.strand == Strand::FORWARD) {
            result.num_forward_reads++;
        } else if (info.strand == Strand::REVERSE) {
            result.num_reverse_reads++;
        }
    }
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
                // Binary threshold
                if (val >= 0.8) {
                    meth_mat.binary_matrix(i, j) = 1;
                } else if (val <= 0.2) {
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

void RegionProcessor::perform_clustering_and_significance(const DistanceMatrix& all_dist,
                                                           const std::vector<ReadInfo>& read_list,
                                                           const MethylationMatrix& meth_mat,
                                                           const std::string& clustering_dir,
                                                           const std::string& chr_name, const SomaticSnv& snv,
                                                           int region_id, RegionResult& result) {
    std::filesystem::create_directories(clustering_dir);

    // Prepare read names
    std::vector<std::string> read_names;
    for (const auto& r : read_list) {
        read_names.push_back(r.read_name);
    }

    // Build hierarchical clustering tree
    HierarchicalClustering clusterer(linkage_method_);
    Tree tree = clusterer.build_tree(all_dist, read_names);

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
        auto [best_k, cluster_labels] = TreeCutter::find_optimal_clusters(
            tree, all_dist.dist_matrix, 2, std::min(6, static_cast<int>(read_list.size()) / 2));

        if (best_k < 2 || cluster_labels.size() != read_list.size()) return;

        // Build FullLabel vector
        std::vector<FullLabel> full_labels;
        full_labels.reserve(read_list.size());
        for (const auto& read : read_list) {
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

        int min_cluster_size = std::max(3, static_cast<int>(read_list.size()) / 20);
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
        sig_config.enable_permanova = false;
        sig_config.enable_dispersion = false;
        sig_config.enable_bootstrap = false;
        sig_config.run_id = vcf_filename_;
        sig_config.vcf_id = vcf_filename_;

        std::string anchor_key = chr_name + "_" + std::to_string(snv.pos);
        SignificanceAnalyzer analyzer(sig_config);
        SignificanceResult sig_result =
            analyzer.analyze(cluster_labels, full_labels, all_dist.dist_matrix, meth_mat.raw_matrix, region_id, anchor_key);

        // Store results
        result.significance_computed = true;
        double p_alt = sig_result.global_alt.valid ? sig_result.global_alt.fisher_ffh.p_value : 1.0;
        double p_hp = sig_result.global_hp.valid ? sig_result.global_hp.fisher_ffh.p_value : 1.0;
        result.global_p_value = std::min(p_alt, p_hp);

        double v_alt = sig_result.global_alt.cramers_v_reliable ? sig_result.global_alt.cramers_v : 0.0;
        double v_hp = sig_result.global_hp.cramers_v_reliable ? sig_result.global_hp.cramers_v : 0.0;
        result.cramers_v = std::max(v_alt, v_hp);
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
        result.coverage_category = determine_coverage_category(result.coverage_multiple);
        result.loh_subtype = determine_loh_subtype(result.potential_loh, result.verification_class);
        result.quality_score = compute_quality_score(
            result.num_reads, result.num_cpgs, result.coverage_multiple,
            result.potential_loh, result.hp_merged_sig, result.allele_sig, result.cramers_v);
        result.quality_tier = determine_quality_tier(result.quality_score);

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
