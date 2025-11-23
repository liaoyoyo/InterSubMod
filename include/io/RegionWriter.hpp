#pragma once

#include <string>
#include <vector>
#include "core/DataStructs.hpp"
#include "core/SomaticSnv.hpp"

namespace InterSubMod {

/**
 * @brief 負責將 region 的所有資料輸出到檔案
 * 
 * 輸出目錄結構：
 * ```
 * output/
 *   region_0000/
 *     metadata.txt          # Region 與 SNV 資訊
 *     reads.tsv             # Read 列表與標籤資訊
 *     cpg_sites.tsv         # CpG 位點列表
 *     methylation.csv       # Read × CpG 甲基化矩陣（CSV格式）
 *     methylation.npz       # (可選) NumPy binary format
 *   region_0001/
 *     ...
 * ```
 */
class RegionWriter {
public:
    /**
     * @brief 建構 RegionWriter
     * @param output_dir 輸出根目錄（如 "output/"）
     */
    explicit RegionWriter(const std::string& output_dir);
    
    /**
     * @brief 寫出一個 region 的完整資料
     * 
     * @param snv SNV 資訊（定義 region 的中心）
     * @param region_id Region 的 ID（用於命名子目錄）
     * @param region_start Region 的起始座標（1-based）
     * @param region_end Region 的結束座標（1-based）
     * @param reads Read 列表（包含標籤資訊）
     * @param cpg_positions CpG 位點列表（1-based，已排序）
     * @param matrix 甲基化矩陣（rows=reads, cols=CpGs，-1.0=no coverage）
     * @param elapsed_ms 處理此 region 的時間（毫秒）
     * @param peak_memory_mb 處理此 region 的峰值記憶體（MB）
     */
    void write_region(
        const SomaticSnv& snv,
        int region_id,
        int32_t region_start,
        int32_t region_end,
        const std::vector<ReadInfo>& reads,
        const std::vector<int32_t>& cpg_positions,
        const std::vector<std::vector<double>>& matrix,
        double elapsed_ms = 0.0,
        double peak_memory_mb = 0.0
    );
    
private:
    std::string output_dir_;
    
    /**
     * @brief 建立 region 子目錄
     * @return region 子目錄的完整路徑
     */
    std::string create_region_dir(int region_id);
    
    /**
     * @brief 寫出 metadata.txt
     */
    void write_metadata(
        const std::string& region_dir,
        const SomaticSnv& snv,
        int region_id,
        int32_t region_start,
        int32_t region_end,
        int num_reads,
        int num_cpgs,
        double elapsed_ms,
        double peak_memory_mb
    );
    
    /**
     * @brief 寫出 reads.tsv
     * 
     * 格式：
     * read_id  read_name  chr  start  end  mapq  hp  ps  alt_support  tags_encoded
     */
    void write_reads(
        const std::string& region_dir,
        const std::vector<ReadInfo>& reads
    );
    
    /**
     * @brief 寫出 cpg_sites.tsv
     * 
     * 格式：
     * cpg_id  chr  position
     */
    void write_cpg_sites(
        const std::string& region_dir,
        int chr_id,
        const std::vector<int32_t>& cpg_positions
    );
    
    /**
     * @brief 寫出 methylation.csv
     * 
     * 格式：CSV，第一行為 CpG 座標，後續每行為一個 read
     * -1.0 用 "NA" 表示（無覆蓋）
     */
    void write_matrix_csv(
        const std::string& region_dir,
        const std::vector<std::vector<double>>& matrix,
        const std::vector<int32_t>& cpg_positions
    );
};

} // namespace InterSubMod

