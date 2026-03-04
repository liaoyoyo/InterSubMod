#pragma once

#include <string>
#include <vector>

#include "core/DataStructs.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/SomaticSnv.hpp"
#include "core/Types.hpp"

namespace InterSubMod {

/**
 * @brief Handles output of all region data to files
 *
 * Output directory structure:
 * ```
 * output/
 *   vcf_filename/                    # VCF filename
 *     chr1/                          # Chromosome name
 *       chr1_12345/                  # Chromosome_SNV_position
 *         chr1_11345_13345/          # Chromosome_start_end
 *           metadata.txt             # Region and SNV information
 *           reads/                   # Read data directory
 *             reads.tsv              # Read list with label info (incl. strand)
 *           methylation/             # Methylation data directory
 *             cpg_sites.tsv          # CpG site list
 *             methylation.csv        # Read x CpG methylation matrix (CSV)
 *             methylation_forward.csv# Forward strand methylation matrix
 *             methylation_reverse.csv# Reverse strand methylation matrix
 *   debug/                           # Debug mode output
 *     filtered_reads.tsv             # Filtered reads and reasons
 * ```
 */
class RegionWriter {
public:
    /**
     * @brief Construct RegionWriter
     * @param output_dir Output root directory (e.g., "output/")
     * @param debug_output_dir Debug output directory (e.g., "output/debug")
     * @param output_strand_matrices Whether to output strand-separated matrices
     * @param vcf_filename VCF filename (used for creating subdirectory)
     */
    explicit RegionWriter(const std::string& output_dir, const std::string& debug_output_dir = "",
                          bool output_strand_matrices = true, const std::string& vcf_filename = "");

    /**
     * @brief Write complete data for a region
     *
     * @param snv SNV information (defines region center)
     * @param region_id Region ID (used for subdirectory naming)
     * @param region_start Region start coordinate (1-based)
     * @param region_end Region end coordinate (1-based)
     * @param reads Read list (including label information)
     * @param cpg_positions CpG site list (1-based, sorted)
     * @param matrix Methylation matrix (rows=reads, cols=CpGs, -1.0=no coverage)
     * @param elapsed_ms Processing time for this region (milliseconds)
     * @param peak_memory_mb Peak memory for this region (MB)
     */
    void write_region(const SomaticSnv& snv, const std::string& chr_name, int region_id, int32_t region_start,
                      int32_t region_end, const std::vector<ReadInfo>& reads, const std::vector<int32_t>& cpg_positions,
                      const std::vector<std::vector<double>>& matrix, double elapsed_ms = 0.0,
                      double peak_memory_mb = 0.0);

    /**
     * @brief Write filtered reads (for Debug mode)
     *
     * @param region_dir Region directory
     * @param chr_name Chromosome name
     * @param filtered_reads List of filtered reads
     */
    void write_filtered_reads(const std::string& region_dir, const std::string& chr_name,
                              const std::vector<FilteredReadInfo>& filtered_reads);

    /**
     * @brief Get region directory path (for external use)
     */
    std::string get_region_dir(const std::string& chr_name, int32_t snv_pos, int32_t region_start, int32_t region_end);

    /**
     * @brief Write distance matrices (including strand-specific matrices)
     *
     * @param region_dir Region directory
     * @param all_matrix Distance matrix for all reads
     * @param forward_matrix Distance matrix for forward strand reads
     * @param reverse_matrix Distance matrix for reverse strand reads
     * @param output_strand_matrices Whether to output strand-specific matrices
     */
    void write_distance_matrices(const std::string& region_dir, const DistanceMatrix& all_matrix,
                                 const DistanceMatrix& forward_matrix, const DistanceMatrix& reverse_matrix,
                                 DistanceMetricType metric, bool output_strand_matrices = true);

    /**
     * @brief Write a single distance matrix
     *
     * @param filepath Output file path
     * @param matrix Distance matrix
     */
    void write_single_distance_matrix(const std::string& filepath, const DistanceMatrix& matrix);

private:
    std::string output_dir_;
    std::string debug_output_dir_;
    bool output_strand_matrices_;
    std::string vcf_filename_;

    /**
     * @brief Create region subdirectory
     * @return Full path of region subdirectory
     */
    std::string create_region_dir(const std::string& chr_name, int32_t snv_pos, int32_t region_start,
                                  int32_t region_end);

    /**
     * @brief Write metadata.txt
     */
    void write_metadata(const std::string& region_dir, const SomaticSnv& snv, const std::string& chr_name,
                        int region_id, int32_t region_start, int32_t region_end, int num_reads, int num_cpgs,
                        int num_forward, int num_reverse, double elapsed_ms, double peak_memory_mb);

    /**
     * @brief Write reads.tsv (with strand information)
     *
     * Format:
     * read_id  read_name  chr  start  end  mapq  hp  alt_support  is_tumor  strand
     */
    void write_reads(const std::string& region_dir, const std::vector<ReadInfo>& reads, const std::string& chr_name);

    /**
     * @brief Write cpg_sites.tsv
     *
     * Format:
     * cpg_id  chr  position
     */
    void write_cpg_sites(const std::string& region_dir, const std::string& chr_name,
                         const std::vector<int32_t>& cpg_positions);

    /**
     * @brief Write methylation.csv
     *
     * Format: CSV, first row is CpG coordinates, subsequent rows are reads
     * -1.0 is represented as "NA" (no coverage)
     */
    void write_matrix_csv(const std::string& region_dir, const std::vector<std::vector<double>>& matrix,
                          const std::vector<int32_t>& cpg_positions);

    /**
     * @brief Write strand-separated methylation matrices
     *
     * @param region_dir Region directory
     * @param reads Read list (used for strand determination)
     * @param matrix Complete methylation matrix
     * @param cpg_positions CpG site list
     */
    void write_strand_matrices(const std::string& region_dir, const std::vector<ReadInfo>& reads,
                               const std::vector<std::vector<double>>& matrix,
                               const std::vector<int32_t>& cpg_positions);

    /**
     * @brief Helper function: Convert Strand to string
     */
    static std::string strand_to_string(Strand s);
};

}  // namespace InterSubMod
