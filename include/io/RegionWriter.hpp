#pragma once

#include <string>
#include <vector>

#include "core/DataStructs.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/SomaticSnv.hpp"
#include "core/Types.hpp"

namespace InterSubMod {

/**
 * @brief 負責將 region 的所有資料輸出到檔案
 *
 * 輸出目錄結構：
 * ```
 * output/
 *   vcf_filename/                    # VCF 檔案名稱
 *     chr1/                          # 染色體名稱
 *       chr1_12345/                  # 染色體_SNV位點
 *         chr1_11345_13345/          # 染色體_起始_結束
 *           metadata.txt             # Region 與 SNV 資訊
 *           reads/                   # Read 資料目錄
 *             reads.tsv              # Read 列表與標籤資訊（含 strand）
 *           methylation/             # 甲基化資料目錄
 *             cpg_sites.tsv          # CpG 位點列表
 *             methylation.csv        # Read × CpG 甲基化矩陣（CSV格式）
 *             methylation_forward.csv# Forward strand 甲基化矩陣
 *             methylation_reverse.csv# Reverse strand 甲基化矩陣
 *   debug/                           # Debug 模式輸出
 *     filtered_reads.tsv             # 被過濾的 reads 與原因
 * ```
 */
class RegionWriter {
public:
    /**
     * @brief 建構 RegionWriter
     * @param output_dir 輸出根目錄（如 "output/"）
     * @param debug_output_dir Debug 輸出目錄（如 "output/debug"）
     * @param output_strand_matrices 是否輸出依 strand 分類的矩陣
     * @param vcf_filename VCF 檔案名稱（用於建立子目錄）
     */
    explicit RegionWriter(const std::string& output_dir, const std::string& debug_output_dir = "",
                          bool output_strand_matrices = true, const std::string& vcf_filename = "");

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
    void write_region(const SomaticSnv& snv, const std::string& chr_name, int region_id, int32_t region_start,
                      int32_t region_end, const std::vector<ReadInfo>& reads, const std::vector<int32_t>& cpg_positions,
                      const std::vector<std::vector<double>>& matrix, double elapsed_ms = 0.0,
                      double peak_memory_mb = 0.0);

    /**
     * @brief 寫出被過濾的 reads（Debug 模式使用）
     *
     * @param region_dir Region 目錄
     * @param chr_name 染色體名稱
     * @param filtered_reads 被過濾的 reads 列表
     */
    void write_filtered_reads(const std::string& region_dir, const std::string& chr_name,
                              const std::vector<FilteredReadInfo>& filtered_reads);

    /**
     * @brief 取得 region 目錄路徑（供外部使用）
     */
    std::string get_region_dir(const std::string& chr_name, int32_t snv_pos, int32_t region_start, int32_t region_end);

    /**
     * @brief 寫出距離矩陣（含 strand-specific 矩陣）
     *
     * @param region_dir Region 目錄
     * @param all_matrix 所有 reads 的距離矩陣
     * @param forward_matrix Forward strand reads 的距離矩陣
     * @param reverse_matrix Reverse strand reads 的距離矩陣
     * @param output_strand_matrices 是否輸出 strand-specific 矩陣
     */
    void write_distance_matrices(const std::string& region_dir, const DistanceMatrix& all_matrix,
                                 const DistanceMatrix& forward_matrix, const DistanceMatrix& reverse_matrix,
                                 DistanceMetricType metric, bool output_strand_matrices = true);

    /**
     * @brief 寫出單一距離矩陣
     *
     * @param filepath 輸出檔案路徑
     * @param matrix 距離矩陣
     */
    void write_single_distance_matrix(const std::string& filepath, const DistanceMatrix& matrix);

private:
    std::string output_dir_;
    std::string debug_output_dir_;
    bool output_strand_matrices_;
    std::string vcf_filename_;

    /**
     * @brief 建立 region 子目錄
     * @return region 子目錄的完整路徑
     */
    std::string create_region_dir(const std::string& chr_name, int32_t snv_pos, int32_t region_start,
                                  int32_t region_end);

    /**
     * @brief 寫出 metadata.txt
     */
    void write_metadata(const std::string& region_dir, const SomaticSnv& snv, const std::string& chr_name,
                        int region_id, int32_t region_start, int32_t region_end, int num_reads, int num_cpgs,
                        int num_forward, int num_reverse, double elapsed_ms, double peak_memory_mb);

    /**
     * @brief 寫出 reads.tsv（含 strand 資訊）
     *
     * 格式：
     * read_id  read_name  chr  start  end  mapq  hp  alt_support  is_tumor  strand
     */
    void write_reads(const std::string& region_dir, const std::vector<ReadInfo>& reads, const std::string& chr_name);

    /**
     * @brief 寫出 cpg_sites.tsv
     *
     * 格式：
     * cpg_id  chr  position
     */
    void write_cpg_sites(const std::string& region_dir, const std::string& chr_name,
                         const std::vector<int32_t>& cpg_positions);

    /**
     * @brief 寫出 methylation.csv
     *
     * 格式：CSV，第一行為 CpG 座標，後續每行為一個 read
     * -1.0 用 "NA" 表示（無覆蓋）
     */
    void write_matrix_csv(const std::string& region_dir, const std::vector<std::vector<double>>& matrix,
                          const std::vector<int32_t>& cpg_positions);

    /**
     * @brief 寫出依 strand 分類的甲基化矩陣
     *
     * @param region_dir Region 目錄
     * @param reads Read 列表（用於判定 strand）
     * @param matrix 完整甲基化矩陣
     * @param cpg_positions CpG 位點列表
     */
    void write_strand_matrices(const std::string& region_dir, const std::vector<ReadInfo>& reads,
                               const std::vector<std::vector<double>>& matrix,
                               const std::vector<int32_t>& cpg_positions);

    /**
     * @brief 輔助函式：將 Strand 轉換為字串
     */
    static std::string strand_to_string(Strand s);
};

}  // namespace InterSubMod
