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
#include "core/SomaticSnv.hpp"
#include "core/TreeStructure.hpp"
#include "io/RegionWriter.hpp"
#include "io/TreeWriter.hpp"
#include "utils/FastaReader.hpp"

namespace InterSubMod {

/**
 * @brief 處理單個 SNV region 的結果統計
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
          avg_common_coverage(0.0) {
    }
};

/**
 * @brief 平行化處理多個 SNV regions 的核心類別
 *
 * 此類別負責：
 * 1. 載入 SNV table
 * 2. 為每個 SNV 定義 region（如 ±2000bp）
 * 3. 使用 OpenMP 平行處理多個 regions
 * 4. 管理 thread-local 資源（BamReader, FastaReader）
 * 5. 收集並報告每個 region 的處理結果
 * 6. 在 debug 模式下記錄被過濾的 reads
 *
 * Thread-safety:
 * - 每個 thread 維護自己的 BAM/FASTA readers
 * - MatrixBuilder 與 RegionWriter 在 critical section 中使用
 * - 結果收集使用 mutex 保護
 */
class RegionProcessor {
public:
    /**
     * @brief 建構 RegionProcessor（簡化版，用於向後相容）
     */
    RegionProcessor(const std::string& tumor_bam_path, const std::string& normal_bam_path,
                    const std::string& ref_fasta_path, const std::string& output_dir, int num_threads = 4,
                    int32_t window_size = 2000);

    /**
     * @brief 建構 RegionProcessor（完整版，使用 Config）
     *
     * @param config Configuration object containing all parameters
     */
    explicit RegionProcessor(const Config& config);

    /**
     * @brief 載入 SNV table（TSV 格式）
     *
     * 格式：chr  pos  ref  alt  qual
     * 範例：chr17  7578000  C  T  100.0
     *
     * @param snv_table_path SNV table 檔案路徑
     * @return 成功載入的 SNV 數量
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
     * @brief 處理所有 SNVs（平行化）
     *
     * @param max_snvs 最多處理幾個 SNVs（0 = 全部）
     * @return 處理結果的 vector
     */
    std::vector<RegionResult> process_all_regions(int max_snvs = 0);

    /**
     * @brief 處理單個 region（由 OpenMP worker thread 呼叫）
     *
     * @param snv SNV 資訊
     * @param region_id Region ID
     * @param bam_reader Thread-local BAM reader
     * @param fasta_reader Thread-local FASTA reader
     * @return RegionResult
     */
    RegionResult process_single_region(const SomaticSnv& snv, int region_id, BamReader& bam_reader,
                                       FastaReader& fasta_reader);

    /**
     * @brief 取得載入的 SNVs 列表
     */
    const std::vector<SomaticSnv>& get_snvs() const {
        return snvs_;
    }

    /**
     * @brief 輸出處理摘要報告
     */
    void print_summary(const std::vector<RegionResult>& results) const;

private:
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

    // Thread-local 資源會在 process_single_region 中建立
};

}  // namespace InterSubMod
