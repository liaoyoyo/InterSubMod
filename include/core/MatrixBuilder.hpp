#pragma once

#include <cstdint>
#include <map>
#include <vector>

#include "core/DataStructs.hpp"
#include "core/MethylationParser.hpp"

namespace InterSubMod {

/**
 * @brief 建構並管理 Read × CpG 甲基化矩陣
 *
 * 此類別負責：
 * 1. 收集所有 reads 的甲基化資訊
 * 2. 建立唯一的 CpG 位點列表（排序後）
 * 3. 建構稀疏矩陣（使用 NaN 表示無數據）
 * 4. 提供矩陣存取與查詢介面
 *
 * 矩陣維度：
 * - Rows: reads (按讀取順序)
 * - Cols: unique CpG sites (按基因組座標排序)
 *
 * 數值含義：
 * - [0.0, 1.0]: 甲基化機率
 * - -1.0: 此 read 未覆蓋此 CpG 位點（使用 -1 代表 NaN）
 *
 * Note: 使用 vector<vector<double>> 而非 Eigen 以避免依賴
 */
class MatrixBuilder {
public:
    MatrixBuilder() = default;

    /**
     * @brief 添加一個 read 的甲基化資訊
     *
     * @param read_info Read 的基本資訊（用於記錄 read metadata）
     * @param methyl_calls 該 read 的所有甲基化 calls
     * @return read_id (矩陣中的 row index)
     */
    int add_read(const ReadInfo& read_info, const std::vector<MethylCall>& methyl_calls);

    /**
     * @brief 完成資料收集，建構最終矩陣
     *
     * 此方法會：
     * 1. 排序所有唯一的 CpG 位點
     * 2. 分配矩陣空間（rows × cols）
     * 3. 填充甲基化數值
     * 4. 設定未覆蓋的位置為 NaN
     */
    void finalize();

    /**
     * @brief 取得最終的甲基化矩陣（只讀）
     * @return 矩陣 reference (rows=reads, cols=CpGs)
     */
    const std::vector<std::vector<double>>& get_matrix() const {
        return matrix_;
    }

    /**
     * @brief 取得所有 Read 資訊（按 row order）
     */
    const std::vector<ReadInfo>& get_reads() const {
        return reads_;
    }

    /**
     * @brief 取得所有 CpG 位點座標（按 column order，已排序）
     */
    const std::vector<int32_t>& get_cpg_positions() const {
        return cpg_positions_;
    }

    /**
     * @brief 取得矩陣維度資訊
     */
    int num_reads() const {
        return reads_.size();
    }
    int num_cpgs() const {
        return cpg_positions_.size();
    }

    /**
     * @brief 清空所有資料（用於處理下一個 region）
     */
    void clear();

private:
    std::vector<ReadInfo> reads_;  ///< Read metadata (row indices)

    // 暫存：read_id -> vector of (pos, prob)
    // Changed from map<int, map> to map<int, vector> for performance
    // Changed from map<int, map> to vector<vector> for performance (O(1) access)
    std::vector<std::vector<std::pair<int32_t, float>>> read_methyl_data_;

    // 最終資料
    std::vector<int32_t> cpg_positions_;       ///< Sorted unique CpG positions (column indices)
    std::vector<std::vector<double>> matrix_;  ///< Final matrix (rows × cols), -1.0 = no coverage

    bool finalized_ = false;
};

}  // namespace InterSubMod
