#pragma once

#include <vector>
#include <string>
#include <memory>
#include "core/SomaticSnv.hpp"
#include "core/BamReader.hpp"
#include "utils/FastaReader.hpp"
#include "core/ReadParser.hpp"
#include "core/MethylationParser.hpp"
#include "core/MatrixBuilder.hpp"
#include "io/RegionWriter.hpp"

namespace InterSubMod {

/**
 * @brief 處理單個 SNV region 的結果統計
 */
struct RegionResult {
    int region_id;
    int snv_id;
    int num_reads;
    int num_cpgs;
    double elapsed_ms;
    double peak_memory_mb;
    bool success;
    std::string error_message;
    
    RegionResult() : region_id(-1), snv_id(-1), num_reads(0), num_cpgs(0),
                     elapsed_ms(0.0), peak_memory_mb(0.0), success(false) {}
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
 * 
 * Thread-safety:
 * - 每個 thread 維護自己的 BAM/FASTA readers
 * - MatrixBuilder 與 RegionWriter 在 critical section 中使用
 * - 結果收集使用 mutex 保護
 */
class RegionProcessor {
public:
    /**
     * @brief 建構 RegionProcessor
     * 
     * @param tumor_bam_path Tumor BAM 檔案路徑
     * @param normal_bam_path Normal BAM 檔案路徑（可選）
     * @param ref_fasta_path 參考基因組 FASTA 路徑
     * @param output_dir 輸出根目錄
     * @param num_threads OpenMP 執行緒數量
     * @param window_size Region 窗口大小（以 SNV 為中心，±window_size）
     */
    RegionProcessor(
        const std::string& tumor_bam_path,
        const std::string& normal_bam_path,
        const std::string& ref_fasta_path,
        const std::string& output_dir,
        int num_threads = 4,
        int32_t window_size = 2000
    );
    
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
     * @return RegionResult
     */
    RegionResult process_single_region(const SomaticSnv& snv, int region_id);
    
    /**
     * @brief 取得載入的 SNVs 列表
     */
    const std::vector<SomaticSnv>& get_snvs() const { return snvs_; }
    
    /**
     * @brief 輸出處理摘要報告
     */
    void print_summary(const std::vector<RegionResult>& results) const;
    
private:
    std::string tumor_bam_path_;
    std::string normal_bam_path_;
    std::string ref_fasta_path_;
    std::string output_dir_;
    int num_threads_;
    int32_t window_size_;
    
    std::vector<SomaticSnv> snvs_;
    ChromIndex chrom_index_;  // Manage chromosome name to ID mapping
    
    // Thread-local 資源會在 process_single_region 中建立
};

} // namespace InterSubMod

