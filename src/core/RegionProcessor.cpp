#include "core/RegionProcessor.hpp"

#include <omp.h>
#include <sys/stat.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include "utils/Logger.hpp"

#include "core/HierarchicalClustering.hpp"
#include "core/MethylationMatrix.hpp"
#include "io/TreeWriter.hpp"

namespace InterSubMod {

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

    omp_set_num_threads(num_threads_);

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
    
    LOG_INFO("Starting processing of " + std::to_string(num_to_process) + " regions...");

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

            // Using ScopedLogger within process_single_region, but we can also log high level start here
            // Note: Excessive locking might preserve order but slow down things. 
            // The Logger handles locking.
            
            // Process the region using the thread-local readers
            results[i] = process_single_region(snv, i, tumor_reader, ref_reader);

            // Log completion status
            if (results[i].success) {
                std::stringstream ss;
                ss << "Region " << i << " (" << chr_name << ":" << snv.pos << ") completed: " 
                   << results[i].num_reads << " reads, " << results[i].elapsed_ms << " ms";
                LOG_INFO(ss.str());
            } else {
                std::stringstream ss;
                ss << "Region " << i << " (" << chr_name << ":" << snv.pos << ") failed: " << results[i].error_message;
                LOG_ERROR(ss.str());
            }
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_elapsed = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    LOG_INFO("All regions processed in " + std::to_string(total_elapsed) + " ms (" + 
             std::to_string(total_elapsed / num_to_process) + " ms/region)");

    return results;
}

RegionResult RegionProcessor::process_single_region(const SomaticSnv& snv, int region_id, BamReader& bam_reader,
                                                    FastaReader& fasta_reader) {
    RegionResult result;
    result.region_id = region_id;
    result.snv_id = snv.snv_id;
    
    // Get chromosome name for logging
    // std::string chr_name = chrom_index_.get_name(snv.chr_id);
    // Utils::ScopedLogger logger("Process Region " + std::to_string(region_id), LogLevel::LOG_DEBUG);
    
    auto t_start = std::chrono::high_resolution_clock::now();

    try {
        // Initialize parsers with filter configuration
        ReadParser read_parser(filter_config_);
        MethylationParser methyl_parser;
        MatrixBuilder matrix_builder;

        // Define region window
        // Ensure we don't go below 1 or beyond chromosome end (checked later)
        int32_t region_start = static_cast<int32_t>(snv.pos) - static_cast<int32_t>(window_size_);
        int32_t region_end = static_cast<int32_t>(snv.pos) + static_cast<int32_t>(window_size_);

        // Get chromosome name and length
        std::string chr_name = chrom_index_.get_name(snv.chr_id);
        int32_t chr_length = fasta_reader.get_chr_length(chr_name);

        // Clamp coordinates
        if (region_start < 1) region_start = 1;
        if (chr_length > 0 && region_end > chr_length) region_end = chr_length;

        // Fetch reads from BAM
        // This uses the thread-local reader, which is much faster than re-opening
        auto reads = bam_reader.fetch_reads(chr_name, region_start, region_end);

        // Fetch reference sequence for this window
        // Needed for CIGAR parsing and CpG verification
        std::string ref_seq = fasta_reader.fetch_sequence(chr_name, region_start, region_end);

        if (ref_seq.empty()) {
            throw std::runtime_error("Failed to fetch reference sequence");
        }

        // Collect filtered reads for debug output
        std::vector<FilteredReadInfo> filtered_reads;

        // Process each read
        int read_count = 0;
        std::unordered_set<std::string> processed_read_names;

        for (auto* b : reads) {
            // 1. Filter and Parse Read Info
            auto [keep, filter_reason] = read_parser.should_keep_with_reason(b);

            if (!keep && !no_filter_output_) {
                // Record filtered read for debug output
                if (output_filtered_reads_) {
                    FilteredReadInfo filtered = read_parser.create_filtered_info(b, true, filter_reason);
                    filtered_reads.push_back(filtered);
                }
                continue;
            }

            ReadInfo info = read_parser.parse(b, read_count, true, snv, ref_seq, region_start);

            // Skip reads that don't cover the SNV or have low quality at the SNV site
            // (unless no_filter_output_ is enabled)
            if (info.alt_support == AltSupport::UNKNOWN && !no_filter_output_) {
                if (output_filtered_reads_) {
                    // Determine the reason for UNKNOWN support
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

            // 2. Parse Methylation (MM/ML tags)
            auto methyl_calls = methyl_parser.parse_read(b, ref_seq, region_start);

            // 3. Add to Matrix Builder
            matrix_builder.add_read(info, methyl_calls);
            read_count++;

            // Count strand
            if (info.strand == Strand::FORWARD) {
                result.num_forward_reads++;
            } else if (info.strand == Strand::REVERSE) {
                result.num_reverse_reads++;
            }
        }

        // Finalize matrix construction (sort CpGs, fill NaNs)
        matrix_builder.finalize();

        result.num_reads = matrix_builder.num_reads();
        result.num_cpgs = matrix_builder.num_cpgs();
        result.num_filtered_reads = filtered_reads.size();

        // Write output to disk
        RegionWriter writer(output_dir_, debug_output_dir_, true, vcf_filename_);
        writer.write_region(snv, chr_name, region_id, region_start, region_end, matrix_builder.get_reads(),
                            matrix_builder.get_cpg_positions(), matrix_builder.get_matrix(),
                            0.0,  // elapsed_ms will be set below
                            0.0   // peak_memory_mb not tracked yet
        );

        // Write filtered reads in debug mode
        if (output_filtered_reads_ && !filtered_reads.empty()) {
            std::string region_dir = writer.get_region_dir(chr_name, snv.pos, region_start, region_end);
            writer.write_filtered_reads(region_dir, chr_name, filtered_reads);
        }

        // Compute and write distance matrix if enabled
        if (compute_distance_matrix_ && result.num_reads >= 2 && result.num_cpgs >= 1) {
            // Build MethylationMatrix for distance calculation
            MethylationMatrix meth_mat;
            meth_mat.region_id = region_id;

            // Get reads and matrix data
            const auto& read_list = matrix_builder.get_reads();
            const auto& raw_matrix = matrix_builder.get_matrix();
            const auto& cpg_positions = matrix_builder.get_cpg_positions();

            // Set up MethylationMatrix
            meth_mat.read_ids.resize(read_list.size());
            for (size_t i = 0; i < read_list.size(); ++i) {
                meth_mat.read_ids[i] = read_list[i].read_id;
            }

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

            // Create distance calculator
            DistanceCalculator dist_calc(distance_config_);

            // Loop over all requested metrics
            for (auto metric : distance_metrics_) {
                // Update config for this metric
                distance_config_.metric = metric;
                dist_calc = DistanceCalculator(distance_config_);

                // Compute all-reads distance matrix
                DistanceMatrix all_dist = dist_calc.compute(meth_mat, read_list);

                // Update result statistics (only for the first metric to avoid overwriting/confusion in summary)
                if (metric == distance_metrics_[0]) {
                    result.num_valid_pairs = all_dist.num_valid_pairs;
                    result.num_invalid_pairs = all_dist.num_invalid_pairs;
                    result.avg_common_coverage = all_dist.avg_common_coverage;
                }

                // Compute strand-specific matrices if enabled
                DistanceMatrix forward_dist, reverse_dist;
                if (output_strand_distance_matrices_) {
                    std::tie(forward_dist, reverse_dist) = dist_calc.compute_strand_specific(meth_mat, read_list);
                }

                // Write distance matrices
                if (output_distance_matrix_) {
                    std::string region_dir = writer.get_region_dir(chr_name, snv.pos, region_start, region_end);
                    writer.write_distance_matrices(region_dir, all_dist, forward_dist, reverse_dist, metric,
                                                   output_strand_distance_matrices_);
                }

                if (log_level_ >= LogLevel::LOG_DEBUG) {
                    std::stringstream ss;
                    ss << "  Distance matrix (" << DistanceCalculator::metric_to_string(metric)
                       << "): " << all_dist.size() << "x" << all_dist.size()
                       << ", valid pairs: " << all_dist.num_valid_pairs << ", avg coverage: " << std::fixed
                       << std::setprecision(1) << all_dist.avg_common_coverage;
                    LOG_DEBUG(ss.str());
                }

                // === Hierarchical Clustering and Tree Output ===
                // Only run for the first metric (typically NHD) to avoid redundant trees
                if (compute_clustering_ && metric == distance_metrics_[0] &&
                    result.num_reads >= clustering_min_reads_) {
                    std::string region_dir = writer.get_region_dir(chr_name, snv.pos, region_start, region_end);
                    std::string clustering_dir = region_dir + "/clustering";
                    std::filesystem::create_directories(clustering_dir);

                    // Prepare read names (using read_name field for identification)
                    std::vector<std::string> read_names;
                    for (const auto& r : read_list) {
                        read_names.push_back(r.read_name);
                    }

                    // Build hierarchical clustering tree
                    HierarchicalClustering clusterer(linkage_method_);
                    Tree tree = clusterer.build_tree(all_dist, read_names);

                    if (!tree.empty() && output_tree_files_) {
                        TreeWriter tree_writer;

                        // Write Newick tree file
                        std::string tree_path = clustering_dir + "/tree.nwk";
                        tree_writer.write_newick(tree, tree_path);

                        // Write linkage matrix (scipy-compatible format)
                        if (output_linkage_matrix_) {
                            std::string linkage_path = clustering_dir + "/linkage_matrix.csv";
                            tree_writer.write_linkage_matrix(tree, linkage_path);
                        }

                        // Write leaf order for Python visualization
                        std::string order_path = clustering_dir + "/leaf_order.txt";
                        std::ofstream order_file(order_path);
                        if (order_file.is_open()) {
                            auto leaves = tree.get_leaves();
                            for (const auto& leaf : leaves) {
                                order_file << leaf->label << "\n";
                            }
                            order_file.close();
                        }

                        if (log_level_ >= LogLevel::LOG_DEBUG) {
                            LOG_DEBUG("  Clustering tree: " + std::to_string(tree.num_leaves()) + " leaves, method=" +
                                      HierarchicalClustering::method_to_string(linkage_method_));
                        }
                    }

                    // Strand-specific trees (if enabled)
                    if (output_strand_distance_matrices_ && output_tree_files_) {
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
                                    TreeWriter tree_writer;
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
                                    TreeWriter tree_writer;
                                    tree_writer.write_newick(rev_tree, clustering_dir + "/tree_reverse.nwk");
                                }
                            }
                        }
                    }
                }
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
       << "Average reads per region: " << (success_count > 0 ? (total_reads / static_cast<double>(success_count)) : 0) << "\n"
       << "Average CpGs per region: " << (success_count > 0 ? (total_cpgs / static_cast<double>(success_count)) : 0) << "\n";

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

}  // namespace InterSubMod
