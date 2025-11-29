# 實作與開發計畫 (詳細版)

## 開發總覽 (Development Overview)

本計畫將「SNV Read 與甲基化資訊擷取」功能拆分為 5 個漸進式階段，每個階段都有明確的目標、實作步驟與驗證方法。預計總開發時間 6 個工作天，最終產出經過完整測試的模組。

**平行化策略**: 從第一階段就考慮執行緒安全設計，到第四階段實現真正的平行處理。

**測試驅動**: 每個階段都包含單元測試，第五階段進行整合測試（32 SNVs, 4 threads）。

---

## 第一階段：基礎架構與 BAM 讀取 (Phase 1) - 第 1-2 天

### 第一階段：目標 (Goal)

建立讀取 BAM 檔案並針對特定區域擷取基本 read 資訊的架構。驗證 HTSlib 的使用正確性與 read 過濾邏輯。

### 第一階段：步驟 (Steps)

#### 1.1 專案設置 (Project Setup)

**任務**:

* 驗證 `CMakeLists.txt` 是否連結 `htslib` 與 `Eigen3`
* 建立新的原始碼檔案

**檔案結構**:

```
include/core/
  ├── BamReader.hpp          # BAM 讀取器介面
  └── ReadParser.hpp         # Read 資訊解析器

src/core/
  ├── BamReader.cpp          # 實作
  └── ReadParser.cpp

tests/
  └── test_bam_reader.cpp    # 單元測試
```

**CMakeLists.txt 檢查項目**:

```cmake
find_package(HTSlib REQUIRED)
find_package(Eigen3 REQUIRED)

target_link_libraries(your_target PRIVATE 
    HTSlib::HTSlib 
    Eigen3::Eigen
)
```

#### 1.2 BamReader 類別實作

**標頭檔** (`include/core/BamReader.hpp`):

```cpp
namespace InterSubMod {

class BamReader {
public:
    BamReader(const std::string& bam_path, int n_threads = 1);
    ~BamReader();
    
    // 禁止拷貝，允許移動
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;
    BamReader(BamReader&&) = default;
    BamReader& operator=(BamReader&&) = default;
    
    // 查詢指定區域的 reads
    std::vector<bam1_t*> fetch_reads(
        const std::string& chr, 
        uint32_t start,  // 0-based
        uint32_t end     // 0-based, exclusive
    );
    
    // 取得 header (用於染色體名稱轉換)
    const bam_hdr_t* get_header() const { return hdr_; }
    
    // 檢查是否成功開啟
    bool is_open() const { return fp_ != nullptr; }

private:
    samFile* fp_;
    bam_hdr_t* hdr_;
    hts_idx_t* idx_;
};

} // namespace InterSubMod
```

**實作重點** (`src/core/BamReader.cpp`):

```cpp
BamReader::BamReader(const std::string& bam_path, int n_threads) {
    fp_ = sam_open(bam_path.c_str(), "r");
    if (!fp_) {
        throw std::runtime_error("Failed to open BAM: " + bam_path);
    }
    
    // 設定執行緒數（用於解壓縮）
    if (n_threads > 1) {
        hts_set_threads(fp_, n_threads);
    }
    
    hdr_ = sam_hdr_read(fp_);
    if (!hdr_) {
        sam_close(fp_);
        throw std::runtime_error("Failed to read BAM header");
    }
    
    idx_ = sam_index_load(fp_, bam_path.c_str());
    if (!idx_) {
        bam_hdr_destroy(hdr_);
        sam_close(fp_);
        throw std::runtime_error("Failed to load BAM index (.bai)");
    }
}

BamReader::~BamReader() {
    if (idx_) hts_idx_destroy(idx_);
    if (hdr_) bam_hdr_destroy(hdr_);
    if (fp_) sam_close(fp_);
}

std::vector<bam1_t*> BamReader::fetch_reads(
    const std::string& chr, 
    uint32_t start, 
    uint32_t end
) {
    std::vector<bam1_t*> reads;
    
    // 建立查詢字串 "chr:start-end"
    std::string region = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    
    hts_itr_t* iter = sam_itr_querys(idx_, hdr_, region.c_str());
    if (!iter) {
        // 區域無效或染色體不存在
        return reads;
    }
    
    bam1_t* b = bam_init1();
    while (sam_itr_next(fp_, iter, b) >= 0) {
        // 分配新的 bam1_t 並複製
        bam1_t* copy = bam_dup1(b);
        reads.push_back(copy);
    }
    
    bam_destroy1(b);
    hts_itr_destroy(iter);
    
    return reads;
}
```

**注意事項**:

* 使用 RAII 管理資源
* 回傳的 `bam1_t*` 需由呼叫者負責釋放
* 執行緒安全：每個執行緒需有自己的 `BamReader` 實例

#### 1.3 ReadParser 類別實作

**功能**: 將 `bam1_t*` 轉換為 `ReadInfo` struct。

**標頭檔** (`include/core/ReadParser.hpp`):

```cpp
namespace InterSubMod {

struct ReadFilterConfig {
    int min_mapq = 20;
    int min_read_length = 1000;
    int min_base_quality = 20;
    bool require_mm_ml = true;
};

class ReadParser {
public:
    ReadParser(const ReadFilterConfig& config = {});
    
    // 過濾函式：檢查 read 是否應該被保留
    bool should_keep(const bam1_t* b) const;
    
    // 轉換為 ReadInfo
    ReadInfo parse(
        const bam1_t* b,
        int read_id,
        bool is_tumor,
        const SomaticSnv& anchor_snv,
        const std::string& ref_seq,
        uint32_t ref_start_pos
    ) const;
    
private:
    ReadFilterConfig config_;
    
    AltSupport determine_alt_support(
        const bam1_t* b,
        const SomaticSnv& snv,
        const std::string& ref_seq,
        uint32_t ref_start_pos
    ) const;
};

} // namespace InterSubMod
```

**實作重點** (`src/core/ReadParser.cpp`):

```cpp
bool ReadParser::should_keep(const bam1_t* b) const {
    // FLAG 檢查
    if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FUNMAP)) {
        return false;
    }
    
    // MAPQ 檢查
    if (b->core.qual < config_.min_mapq) {
        return false;
    }
    
    // Read 長度檢查
    int read_len = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
    if (read_len < config_.min_read_length) {
        return false;
    }
    
    // MM/ML tags 檢查
    if (config_.require_mm_ml) {
        if (!bam_aux_get(b, "MM") || !bam_aux_get(b, "ML")) {
            return false;
        }
    }
    
    return true;
}

ReadInfo ReadParser::parse(
    const bam1_t* b,
    int read_id,
    bool is_tumor,
    const SomaticSnv& anchor_snv,
    const std::string& ref_seq,
    uint32_t ref_start_pos
) const {
    ReadInfo info;
    info.read_id = read_id;
    info.read_name = bam_get_qname(b);
    info.chr_id = anchor_snv.chr_id;
    info.align_start = b->core.pos;
    info.align_end = bam_endpos(b);
    info.mapq = b->core.qual;
    info.is_tumor = is_tumor;
    
    // HP tag
    uint8_t* hp_aux = bam_aux_get(b, "HP");
    info.hp_tag = hp_aux ? bam_aux2i(hp_aux) : 0;
    
    // AltSupport
    info.alt_support = determine_alt_support(b, anchor_snv, ref_seq, ref_start_pos);
    
    return info;
}

AltSupport ReadParser::determine_alt_support(
    const bam1_t* b,
    const SomaticSnv& snv,
    const std::string& ref_seq,
    uint32_t ref_start_pos
) const {
    // 詳細實作見設計文件 section 3.2.2
    // 此處為簡化版本
    
    uint32_t snv_pos_0based = snv.pos - 1;
    uint32_t read_start = b->core.pos;
    uint32_t read_end = bam_endpos(b);
    
    // 1. 檢查覆蓋
    if (snv_pos_0based < read_start || snv_pos_0based >= read_end) {
        return AltSupport::UNKNOWN;
    }
    
    // 2. 透過 CIGAR 找到 SNV 在 read 中的位置
    // (省略詳細實作，見設計文件)
    // ...
    
    return AltSupport::UNKNOWN;  // placeholder
}
```

#### 1.4 單元測試

**測試檔案** (`tests/test_bam_reader.cpp`):

```cpp
#include <gtest/gtest.h>
#include "core/BamReader.hpp"
#include "core/ReadParser.hpp"

using namespace InterSubMod;

TEST(BamReaderTest, OpenAndFetch) {
    // 使用測試資料
    BamReader reader("/path/to/test.bam");
    ASSERT_TRUE(reader.is_open());
    
    // 擷取已知區域（例如 chr17:7577000-7578000）
    auto reads = reader.fetch_reads("chr17", 7577000, 7578000);
    
    // 驗證：與 samtools view -c 的結果比對
    // samtools view -c test.bam chr17:7577000-7578000
    int expected_count = 50;  // 假設值，實際需測量
    EXPECT_EQ(reads.size(), expected_count);
    
    // 清理
    for (auto* b : reads) {
        bam_destroy1(b);
    }
}

TEST(ReadParserTest, Filtering) {
    ReadParser parser;
    
    // 建立模擬 bam1_t (需要手動設定 FLAG, MAPQ 等)
    // ...
    
    // 測試過濾邏輯
    // EXPECT_FALSE(parser.should_keep(secondary_read));
    // EXPECT_TRUE(parser.should_keep(good_read));
}
```

**驗證方法**:

```bash
# 1. 使用 samtools 確認 read 數量
samtools view -c tumor.bam chr17:7577000-7578000

# 2. 執行測試
./build/tests/test_bam_reader

# 3. 檢查輸出
# Expected: 測試通過，讀取的 read 數量與 samtools 一致
```

### 第一階段：交付成果 (Deliverables)

* [x] `BamReader` class 實作完成
* [x] `ReadParser` class 實作完成（AltSupport 判斷可為 placeholder）
* [x] 單元測試通過
* [x] 文件：撰寫 `docs/dev/phase1_report.md`，記錄遇到的問題與解決方法

---

## 第二階段：甲基化與 CIGAR 解析 (Phase 2) - 第 2-3 天

### 第二階段：目標 (Goal)

正確地將 `MM/ML` tags 映射至 hg38 基因組座標，完成 CIGAR 解析邏輯，並進行詳細驗證。

### 第二階段：步驟 (Steps)

#### 2.1 參考基因組載入器 (FastaReader)

**標頭檔** (`include/utils/FastaReader.hpp`):

```cpp
namespace InterSubMod {

class FastaReader {
public:
    explicit FastaReader(const std::string& fasta_path);
    ~FastaReader();
    
    // 取得指定區域的序列
    std::string fetch_sequence(
        const std::string& chr,
        uint32_t start,  // 0-based
        uint32_t end     // 0-based, exclusive
    );
    
    // 檢查是否成功載入
    bool is_loaded() const { return fai_ != nullptr; }

private:
    faidx_t* fai_;
};

} // namespace InterSubMod
```

**實作** (`src/utils/FastaReader.cpp`):

```cpp
FastaReader::FastaReader(const std::string& fasta_path) {
    fai_ = fai_load(fasta_path.c_str());
    if (!fai_) {
        throw std::runtime_error("Failed to load FASTA index: " + fasta_path + ".fai");
    }
}

FastaReader::~FastaReader() {
    if (fai_) fai_destroy(fai_);
}

std::string FastaReader::fetch_sequence(
    const std::string& chr,
    uint32_t start,
    uint32_t end
) {
    int len;
    char* seq = faidx_fetch_seq(fai_, chr.c_str(), start, end - 1, &len);
    if (!seq) {
        return "";  // 區域無效
    }
    
    std::string result(seq, len);
    free(seq);
    
    // 轉為大寫
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    
    return result;
}
```

#### 2.2 甲基化解析器 (MethylationParser)

**標頭檔** (`include/core/MethylationParser.hpp`):

```cpp
namespace InterSubMod {

struct MethylCall {
    uint32_t ref_pos;      // 1-based, hg38 座標
    float probability;    // 0.0 - 1.0
};

class MethylationParser {
public:
    // 解析單一 read 的甲基化資訊
    std::vector<MethylCall> parse_read(
        const bam1_t* b,
        const std::string& ref_seq,  // 該 read 涵蓋範圍的參考序列
        uint32_t ref_start_pos        // ref_seq 的起始座標 (0-based)
    );
    
private:
    // 解析 MM tag，取得 delta-encoded skip counts
    std::vector<int> parse_mm_tag(const char* mm_str, const std::string& mod_code = "C+m?");
    
    // 建立 read sequence index -> reference position 映射
    std::vector<int32_t> build_seq_to_ref_map(const bam1_t* b);
    
    // 驗證 CpG context
    bool is_cpg_site(const std::string& ref_seq, int offset);
};

} // namespace InterSubMod
```

**實作重點** (`src/core/MethylationParser.cpp`):

```cpp
std::vector<MethylCall> MethylationParser::parse_read(
    const bam1_t* b,
    const std::string& ref_seq,
    uint32_t ref_start_pos
) {
    std::vector<MethylCall> calls;
    
    // 1. 取得 MM 與 ML tags
    uint8_t* mm_aux = bam_aux_get(b, "MM");
    uint8_t* ml_aux = bam_aux_get(b, "ML");
    if (!mm_aux || !ml_aux) return calls;
    
    const char* mm_str = bam_aux2Z(mm_aux);
    
    // 2. 解析 ML array
    // ML 格式: 'B' 'C' [len:4 bytes] [data:len bytes]
    if (ml_aux[0] != 'B' || ml_aux[1] != 'C') {
        return calls;  // 格式錯誤
    }
    uint32_t ml_len = le_to_u32(ml_aux + 2);
    const uint8_t* ml_data = ml_aux + 6;
    
    // 3. 解析 MM deltas
    std::vector<int> deltas = parse_mm_tag(mm_str, "C+m?");
    if (deltas.size() != ml_len) {
        // 長度不匹配，可能是錯誤的 MM tag 或包含多種修飾
        return calls;
    }
    
    // 4. 建立 seq -> ref 映射
    std::vector<int32_t> seq_to_ref = build_seq_to_ref_map(b);
    
    // 5. 遍歷 read sequence，找到有甲基化資訊的 'C'
    uint8_t* seq = bam_get_seq(b);
    int seq_len = b->core.l_qseq;
    
    int c_count = 0;       // read 中的 'C' 計數
    int delta_idx = 0;     // deltas 索引
    int next_c_target = (deltas.size() > 0) ? deltas[0] : -1;
    
    for (int seq_idx = 0; seq_idx < seq_len; seq_idx++) {
        char base = seq_nt16_str[bam_seqi(seq, seq_idx)];
        
        if (base == 'C') {
            if (c_count == next_c_target) {
                // 這個 'C' 有甲基化資訊
                int32_t ref_pos = seq_to_ref[seq_idx];
                
                if (ref_pos >= 0) {  // 非 insertion
                    // 驗證參考序列
                    int ref_offset = ref_pos - ref_start_pos;
                    if (ref_offset >= 0 && ref_offset < ref_seq.size()) {
                        if (ref_seq[ref_offset] == 'C' && is_cpg_site(ref_seq, ref_offset)) {
                            float prob = ml_data[delta_idx] / 255.0f;
                            calls.push_back({ref_pos + 1, prob});  // 轉為 1-based
                        }
                    }
                }
                
                // 移動到下一個 delta
                delta_idx++;
                if (delta_idx < deltas.size()) {
                    next_c_target += deltas[delta_idx] + 1;
                } else {
                    next_c_target = -1;
                }
            }
            c_count++;
        }
    }
    
    return calls;
}

std::vector<int> MethylationParser::parse_mm_tag(const char* mm_str, const std::string& mod_code) {
    std::vector<int> deltas;
    
    // MM 格式範例: "C+m?,3,5,0,2;C+h?,..."
    // 我們需要找到 "C+m?" 段落
    
    std::string mm(mm_str);
    size_t pos = mm.find(mod_code);
    if (pos == std::string::npos) return deltas;
    
    // 跳過 "C+m?,"
    pos += mod_code.length();
    if (pos < mm.length() && mm[pos] == ',') pos++;
    
    // 解析數字直到遇到 ';'
    std::stringstream ss(mm.substr(pos));
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty() || token[0] == ';') break;
        try {
            deltas.push_back(std::stoi(token));
        } catch (...) {
            break;
        }
    }
    
    return deltas;
}

std::vector<int32_t> MethylationParser::build_seq_to_ref_map(const bam1_t* b) {
    std::vector<int32_t> seq_to_ref(b->core.l_qseq, -1);
    
    int32_t ref_pos = b->core.pos;  // 0-based
    int seq_pos = 0;
    
    uint32_t* cigar = bam_get_cigar(b);
    for (int i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        
        switch (op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                for (int j = 0; j < len; j++) {
                    seq_to_ref[seq_pos++] = ref_pos++;
                }
                break;
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                seq_pos += len;  // seq 前進，ref 不動
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                ref_pos += len;  // ref 前進，seq 不動
                break;
            case BAM_CHARD_CLIP:
                // 不計入 seq
                break;
        }
    }
    
    return seq_to_ref;
}

bool MethylationParser::is_cpg_site(const std::string& ref_seq, int offset) {
    if (offset + 1 >= ref_seq.size()) return false;
    return (ref_seq[offset] == 'C' && ref_seq[offset + 1] == 'G');
}
```

#### 2.3 補完 AltSupport 判斷邏輯

在 `ReadParser::determine_alt_support` 中實作完整的 CIGAR 遍歷邏輯（參考設計文件 3.2.2 節）。

#### 2.4 驗證測試

**測試策略**:

1. **選擇測試 Read**:

    ```bash
    # 從 BAM 中選擇一條已知的 read
    samtools view tumor.bam chr17:7577000-7578000 | head -1
    ```

2. **手動驗證 CIGAR**:
    * 記錄 read 的 POS 與 CIGAR string
    * 計算 SNV 位置應對應的 read offset
    * 檢查程式輸出是否一致

3. **手動驗證 MM/ML**:
    * 使用 `samtools view` 顯示 MM/ML tags
    * 手動計算第一個 CpG 的參考座標
    * 檢查程式輸出的 `MethylCall.ref_pos` 是否正確

4. **CpG 驗證**:

    ```bash
    # 使用 samtools faidx 取得序列
    samtools faidx hg38.fa chr17:7578000-7578002
    # 輸出應為 "CG"
    ```

**測試檔案** (`tests/test_methylation_parser.cpp`):

```cpp
TEST(MethylationParserTest, ParseKnownRead) {
    // 使用已知的 read (read_name = "xyz")
    BamReader reader("tumor.bam");
    auto reads = reader.fetch_reads("chr17", 7577000, 7578000);
    
    // 找到特定 read
    bam1_t* target_read = nullptr;
    for (auto* r : reads) {
        if (std::string(bam_get_qname(r)) == "target_read_name") {
            target_read = r;
            break;
        }
    }
    ASSERT_NE(target_read, nullptr);
    
    // 取得參考序列
    FastaReader fasta("hg38.fa");
    std::string ref_seq = fasta.fetch_sequence("chr17", 7577000, 7579000);
    
    // 解析甲基化
    MethylationParser parser;
    auto calls = parser.parse_read(target_read, ref_seq, 7577000);
    
    // 驗證：至少應有一些 CpG
    EXPECT_GT(calls.size(), 0);
    
    // 驗證第一個 CpG 的座標（手動計算的預期值）
    int expected_first_cpg = 7577050;  // 假設值
    EXPECT_EQ(calls[0].ref_pos, expected_first_cpg);
    
    // 驗證機率範圍
    for (const auto& call : calls) {
        EXPECT_GE(call.probability, 0.0f);
        EXPECT_LE(call.probability, 1.0f);
    }
    
    // 清理
    for (auto* r : reads) bam_destroy1(r);
}
```

**輸出除錯資訊**:

```cpp
// 在 MethylationParser 中加入 debug 輸出
void MethylationParser::debug_print(const bam1_t* b, const std::vector<MethylCall>& calls) {
    std::cout << "Read: " << bam_get_qname(b) << "\n";
    std::cout << "  Position: " << b->core.pos << " - " << bam_endpos(b) << "\n";
    std::cout << "  CIGAR: ";
    // 印出 CIGAR string
    std::cout << "\n  Methylation sites found: " << calls.size() << "\n";
    for (const auto& call : calls) {
        std::cout << "    CpG " << call.ref_pos << ": " << call.probability << "\n";
    }
}
```

### 第二階段：交付成果 (Deliverables)

* [x] `FastaReader` class 實作完成
* [x] `MethylationParser` class 實作完成
* [x] `ReadParser::determine_alt_support` 完整實作
* [x] CIGAR 解析邏輯通過手動驗證（至少 3 個 reads）
* [x] CpG 座標驗證通過（與參考序列比對）
* [x] 文件：`docs/dev/phase2_cigar_verification.md`，記錄驗證過程與結果

---

## 第三階段：矩陣建構與輸出 (Phase 3) - 第 4 天

### 第三階段：目標 (Goal)

將解析後的 reads 與 CpG 資訊聚合為矩陣，並輸出為標準化的檔案格式，便於後續驗證與分析。

### 第三階段：步驟 (Steps)

#### 3.1 矩陣建構器 (MatrixBuilder)

**標頭檔** (`include/core/MatrixBuilder.hpp`):

```cpp
namespace InterSubMod {

class MatrixBuilder {
public:
    // 建構 region 的甲基化矩陣
    MethylationMatrix build(
        const RegionOfInterest& region,
        const std::vector<ReadInfo>& reads,
        const std::vector<std::vector<MethylCall>>& methyl_calls,  // per-read
        bool use_ref_cpgs = true  // 是否使用參考基因組的 CpG 作為欄位
    );
    
private:
    // 從參考序列提取所有 CpG 位點
    std::vector<CpGSite> extract_ref_cpgs(
        const std::string& chr,
        int chr_id,
        uint32_t start,
        uint32_t end,
        const std::string& ref_seq
    );
    
    // 二值化甲基化機率
    int binarize(float probability, float high_thresh = 0.8f, float low_thresh = 0.2f);
};

} // namespace InterSubMod
```

**實作** (`src/core/MatrixBuilder.cpp`):

```cpp
MethylationMatrix MatrixBuilder::build(
    const RegionOfInterest& region,
    const std::vector<ReadInfo>& reads,
    const std::vector<std::vector<MethylCall>>& methyl_calls,
    bool use_ref_cpgs
) {
    MethylationMatrix matrix;
    matrix.region_id = region.region_id;
    
    // 1. 收集 CpG 位點
    std::set<uint32_t> cpg_positions;
    
    if (use_ref_cpgs) {
        // 從參考基因組提取（更完整）
        FastaReader fasta("hg38.fa");  // 應該作為成員變數或參數傳入
        std::string chr_name = "chr17";  // 應從 ChromIndex 取得
        std::string ref_seq = fasta.fetch_sequence(chr_name, region.win_start_pos, region.win_end_pos);
        
        for (size_t i = 0; i + 1 < ref_seq.size(); i++) {
            if (ref_seq[i] == 'C' && ref_seq[i+1] == 'G') {
                cpg_positions.insert(region.win_start_pos + i + 1);  // 1-based
            }
        }
    } else {
        // 僅使用 reads 中發現的 CpG
        for (const auto& calls : methyl_calls) {
            for (const auto& call : calls) {
                cpg_positions.insert(call.ref_pos);
            }
        }
    }
    
    // 2. 建立 CpGSite 列表
    std::vector<CpGSite> sites;
    int cpg_id = 0;
    for (uint32_t pos : cpg_positions) {
        sites.push_back({cpg_id++, region.chr_id, pos, false, false, false});
    }
    
    // 3. 建立位置 -> 欄位索引映射
    std::unordered_map<int32_t, int> pos_to_col;
    for (size_t col = 0; col < sites.size(); col++) {
        pos_to_col[sites[col].pos] = col;
    }
    
    // 4. 初始化矩陣
    int n_reads = reads.size();
    int n_sites = sites.size();
    
    matrix.read_ids.resize(n_reads);
    matrix.cpg_ids.resize(n_sites);
    
    for (int i = 0; i < n_reads; i++) {
        matrix.read_ids[i] = reads[i].read_id;
    }
    for (int i = 0; i < n_sites; i++) {
        matrix.cpg_ids[i] = sites[i].cpg_id;
    }
    
    matrix.raw_matrix = Eigen::MatrixXd::Constant(n_reads, n_sites, NAN);
    matrix.binary_matrix = Eigen::MatrixXi::Constant(n_reads, n_sites, -1);
    
    // 5. 填入數值
    for (int row = 0; row < n_reads; row++) {
        for (const auto& call : methyl_calls[row]) {
            auto it = pos_to_col.find(call.ref_pos);
            if (it != pos_to_col.end()) {
                int col = it->second;
                matrix.raw_matrix(row, col) = call.probability;
                matrix.binary_matrix(row, col) = binarize(call.probability);
            }
        }
    }
    
    return matrix;
}

int MatrixBuilder::binarize(float probability, float high_thresh, float low_thresh) {
    if (probability >= high_thresh) return 1;
    if (probability <= low_thresh) return 0;
    return -1;  // uncertain
}
```

#### 3.2 輸出格式化器 (RegionWriter)

**標頭檔** (`include/io/RegionWriter.hpp`):

```cpp
namespace InterSubMod {

class RegionWriter {
public:
    explicit RegionWriter(const std::string& output_dir);
    
    // 寫入單一 region 的所有資料
    void write_region(
        const RegionOfInterest& region,
        const MethylationMatrix& matrix,
        const std::vector<ReadInfo>& reads,
        const std::vector<CpGSite>& sites
    );
    
private:
    std::string output_dir_;
    
    void write_matrix_csv(const std::string& path, const Eigen::MatrixXd& mat, const std::vector<int>& row_ids);
    void write_matrix_csv(const std::string& path, const Eigen::MatrixXi& mat, const std::vector<int>& row_ids);
    void write_reads_jsonl(const std::string& path, const std::vector<ReadInfo>& reads);
    void write_sites_csv(const std::string& path, const std::vector<CpGSite>& sites);
};

} // namespace InterSubMod
```

**實作** (`src/io/RegionWriter.cpp`):

```cpp
void RegionWriter::write_region(
    const RegionOfInterest& region,
    const MethylationMatrix& matrix,
    const std::vector<ReadInfo>& reads,
    const std::vector<CpGSite>& sites
) {
    // 建立目錄
    std::string region_dir = output_dir_ + "/region_" + std::to_string(region.region_id);
    std::filesystem::create_directories(region_dir);
    
    // 寫入檔案
    write_matrix_csv(region_dir + "/matrix_raw.csv", matrix.raw_matrix, matrix.read_ids);
    write_matrix_csv(region_dir + "/matrix_binary.csv", matrix.binary_matrix, matrix.read_ids);
    write_reads_jsonl(region_dir + "/reads.jsonl", reads);
    write_sites_csv(region_dir + "/sites.csv", sites);
}

void RegionWriter::write_matrix_csv(
    const std::string& path, 
    const Eigen::MatrixXd& mat, 
    const std::vector<int>& row_ids
) {
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Cannot open: " + path);
    
    // 表頭
    ofs << "read_id";
    for (int col = 0; col < mat.cols(); col++) {
        ofs << ",cpg_" << col;
    }
    ofs << "\n";
    
    // 資料
    for (int row = 0; row < mat.rows(); row++) {
        ofs << row_ids[row];
        for (int col = 0; col < mat.cols(); col++) {
            ofs << ",";
            double val = mat(row, col);
            if (std::isnan(val)) {
                ofs << "NaN";
            } else {
                ofs << val;
            }
        }
        ofs << "\n";
    }
}

void RegionWriter::write_reads_jsonl(const std::string& path, const std::vector<ReadInfo>& reads) {
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Cannot open: " + path);
    
    for (const auto& read : reads) {
        ofs << "{";
        ofs << "\"read_id\":" << read.read_id << ",";
        ofs << "\"read_name\":\"" << read.read_name << "\",";
        ofs << "\"chr_id\":" << read.chr_id << ",";
        ofs << "\"align_start\":" << read.align_start << ",";
        ofs << "\"align_end\":" << read.align_end << ",";
        ofs << "\"mapq\":" << read.mapq << ",";
        ofs << "\"hp\":" << read.hp_tag << ",";
        ofs << "\"is_tumor\":" << (read.is_tumor ? "true" : "false") << ",";
        ofs << "\"alt_support\":\"";
        switch (read.alt_support) {
            case AltSupport::ALT: ofs << "ALT"; break;
            case AltSupport::REF: ofs << "REF"; break;
            case AltSupport::UNKNOWN: ofs << "UNKNOWN"; break;
        }
        ofs << "\"}\n";
    }
}

void RegionWriter::write_sites_csv(const std::string& path, const std::vector<CpGSite>& sites) {
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Cannot open: " + path);
    
    ofs << "cpg_id,chr_id,pos\n";
    for (const auto& site : sites) {
        ofs << site.cpg_id << "," << site.chr_id << "," << site.pos << "\n";
    }
}
```

#### 3.3 測試

**單一 Region 測試**:

```cpp
TEST(MatrixBuilderTest, BuildSmallMatrix) {
    // 模擬 3 reads, 5 CpG sites
    RegionOfInterest region{0, 0, 1, 7577000, 7579000};
    
    std::vector<ReadInfo> reads(3);
    // ... 填入測試資料
    
    std::vector<std::vector<MethylCall>> methyl_calls(3);
    // read 0: CpG at 7577050 (prob=0.9), 7577100 (prob=0.1)
    methyl_calls[0] = {{7577050, 0.9f}, {7577100, 0.1f}};
    // read 1: CpG at 7577050 (prob=0.85)
    methyl_calls[1] = {{7577050, 0.85f}};
    // read 2: 無甲基化資訊
    methyl_calls[2] = {};
    
    MatrixBuilder builder;
    auto matrix = builder.build(region, reads, methyl_calls, false);
    
    // 驗證
    EXPECT_EQ(matrix.raw_matrix.rows(), 3);
    EXPECT_EQ(matrix.raw_matrix.cols(), 2);  // 2 unique CpG
    
    // read 0, cpg 0: 0.9
    EXPECT_NEAR(matrix.raw_matrix(0, 0), 0.9, 0.01);
    // read 0, cpg 1: 0.1
    EXPECT_NEAR(matrix.raw_matrix(0, 1), 0.1, 0.01);
    // read 1, cpg 1: NaN (未覆蓋)
    EXPECT_TRUE(std::isnan(matrix.raw_matrix(1, 1)));
    // read 2, all: NaN
    EXPECT_TRUE(std::isnan(matrix.raw_matrix(2, 0)));
}
```

### 第三階段：交付成果 (Deliverables)

* [x] `MatrixBuilder` class 實作完成
* [x] `RegionWriter` class 實作完成
* [x] 輸出檔案格式符合規範
* [x] 單元測試通過（小規模矩陣）
* [x] 手動檢查輸出的 CSV/JSONL 檔案

---

## 第四階段：平行化與驅動程式 (Phase 4) - 第 5 天

### 第四階段：目標 (Goal)

整合前三階段的模組，建立完整的 region 處理流程，並實現 OpenMP 平行化。

### 第四階段：步驟 (Steps)

#### 4.1 區域處理器 (RegionProcessor)

封裝單一 region 的完整處理流程。

**標頭檔** (`include/core/RegionProcessor.hpp`):

```cpp
namespace InterSubMod {

class RegionProcessor {
public:
    RegionProcessor(
        const std::string& tumor_bam,
        const std::string& normal_bam,
        const std::string& ref_fasta,
        const std::string& output_dir
    );
    
    // 處理單一 region
    void process_region(const RegionOfInterest& region, const SomaticSnv& snv);
    
private:
    std::string tumor_bam_;
    std::string normal_bam_;
    std::string ref_fasta_;
    std::string output_dir_;
    
    // Thread-local 資源（在 process_region 中建立）
    // BamReader, FastaReader 等
};

} // namespace InterSubMod
```

**實作** (`src/core/RegionProcessor.cpp`):

```cpp
void RegionProcessor::process_region(const RegionOfInterest& region, const SomaticSnv& snv) {
    // 1. 建立 thread-local 資源
    thread_local BamReader tumor_reader(tumor_bam_);
    thread_local BamReader normal_reader(normal_bam_);
    thread_local FastaReader fasta_reader(ref_fasta_);
    thread_local ReadParser read_parser;
    thread_local MethylationParser methyl_parser;
    
    // 2. 取得參考序列
    std::string chr_name = "chr17";  // 應從 ChromIndex 取得
    std::string ref_seq = fasta_reader.fetch_sequence(
        chr_name, region.win_start_pos, region.win_end_pos
    );
    
    // 3. 擷取 reads (tumor + normal)
    std::vector<ReadInfo> all_reads;
    std::vector<std::vector<MethylCall>> all_methyl_calls;
    int read_id = 0;
    
    // Tumor reads
    auto tumor_bams = tumor_reader.fetch_reads(chr_name, region.win_start_pos, region.win_end_pos);
    for (auto* b : tumor_bams) {
        if (!read_parser.should_keep(b)) {
            bam_destroy1(b);
            continue;
        }
        
        ReadInfo info = read_parser.parse(b, read_id++, true, snv, ref_seq, region.win_start_pos);
        auto methyl_calls = methyl_parser.parse_read(b, ref_seq, region.win_start_pos);
        
        all_reads.push_back(info);
        all_methyl_calls.push_back(methyl_calls);
        
        bam_destroy1(b);
    }
    
    // Normal reads (類似)
    // ...
    
    // 4. 建構矩陣
    MatrixBuilder builder;
    auto matrix = builder.build(region, all_reads, all_methyl_calls);
    
    // 5. 輸出
    RegionWriter writer(output_dir_);
    // 需要重建 CpGSite 列表（從 matrix.cpg_ids 與座標）
    std::vector<CpGSite> sites;  // ... 填入
    writer.write_region(region, matrix, all_reads, sites);
}
```

#### 4.2 主驅動程式 (Main Driver)

**檔案** (`src/main_extract_reads.cpp`):

```cpp
#include <omp.h>
#include "core/SomaticSnv.hpp"
#include "core/RegionProcessor.hpp"

int main(int argc, char** argv) {
    // 1. 載入 SNV table
    ChromIndex chrom_index;
    SomaticSnvTable snv_table;
    snv_table.load_from_vcf("somatic.vcf", chrom_index);
    
    // 2. 產生 Regions（假設每個 SNV 產生一個 region，±2000bp）
    std::vector<RegionOfInterest> regions;
    int region_id = 0;
    for (const auto& snv : snv_table.all()) {
        RegionOfInterest region;
        region.region_id = region_id++;
        region.snv_id = snv.snv_id;
        region.chr_id = snv.chr_id;
        region.win_start_pos = std::max(0, snv.pos - 2000);
        region.win_end_pos = snv.pos + 2000;
        regions.push_back(region);
    }
    
    // 3. 平行處理
    int n_threads = 4;
    omp_set_num_threads(n_threads);
    
    std::cout << "Processing " << regions.size() << " regions with " << n_threads << " threads...\n";
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < regions.size(); i++) {
        RegionProcessor processor(
            "tumor.bam",
            "normal.bam",
            "hg38.fa",
            "output/"
        );
        
        const auto& region = regions[i];
        const auto& snv = snv_table.all()[region.snv_id];
        
        processor.process_region(region, snv);
        
        #pragma omp critical
        {
            std::cout << "Completed region " << region.region_id << "\n";
        }
    }
    
    std::cout << "All regions processed.\n";
    return 0;
}
```

#### 4.3 關鍵設計考量

**執行緒安全**:

* 每個執行緒必須有自己的 `BamReader`, `FastaReader` 等實例
* 使用 `thread_local` 關鍵字
* 輸出檔案不會衝突（因為每個 region 有獨立的目錄）

**記憶體管理**:

* 在迴圈內明確釋放 `bam1_t*`
* 避免在 parallel region 內建立大型全域資料結構

**負載平衡**:

* 使用 `schedule(dynamic)` 讓較快的執行緒處理更多 regions

### 第四階段：交付成果 (Deliverables)

* [x] `RegionProcessor` class 實作完成
* [x] 主驅動程式實作完成
* [x] 編譯無錯誤，可執行
* [x] 小規模測試（2-3 個 regions, 2 threads）確認平行化正常

---

## 第五階段：「32 SNV」驗證執行 (Phase 5) - 第 6 天

### 第五階段：目標 (Goal)

執行完整的測試案例（前 32 個 SNVs, 4 執行緒），收集性能數據，驗證正確性，撰寫最終報告。

### 第五階段：步驟 (Steps)

#### 5.1 配置測試

**修改主程式**:

```cpp
// 限制只處理前 32 個 SNVs
int max_regions = 32;
regions.resize(std::min(regions.size(), (size_t)max_regions));
```

**設定執行緒數**:

```cpp
int n_threads = 4;
omp_set_num_threads(n_threads);
```

#### 5.2 執行與監控

**執行指令**:

```bash
# 清空輸出目錄
rm -rf output/region_*

# 執行（記錄時間）
time ./build/extract_reads \
    --vcf data/vcf/somatic.vcf \
    --tumor data/bam/tumor.bam \
    --normal data/bam/normal.bam \
    --ref data/ref/hg38.fa \
    --output output/ \
    --threads 4 \
    --max-regions 32

# 監控記憶體（另一個 terminal）
watch -n 1 'ps aux | grep extract_reads | grep -v grep'
```

**資源監控**:

```cpp
// 在主程式中加入
#include "utils/ResourceMonitor.hpp"

int main() {
    ResourceMonitor monitor;
    monitor.start();
    
    // ... 處理 regions ...
    
    monitor.stop();
    monitor.print_report("output/resource_usage.txt");
}
```

#### 5.3 正確性驗證

**檢查清單**:

1. **檔案完整性**:

```bash
# 檢查是否所有 32 個 region 都有輸出
ls output/region_*/matrix_raw.csv | wc -l
# 預期: 32
```

2. **Read 數量驗證** (隨機抽樣 3 個 regions):

```bash
# Region 0: chr17:7577000-7580000
samtools view -c tumor.bam chr17:7577000-7580000
# 比對 output/region_0/reads.jsonl 的行數
wc -l output/region_0/reads.jsonl
```

3. **矩陣維度檢查**:

```python
import pandas as pd
import numpy as np

# 讀取 matrix_raw.csv
df = pd.read_csv("output/region_0/matrix_raw.csv")
print(f"Shape: {df.shape}")  # (n_reads, n_cpgs+1)

# 檢查 NaN 分布
nan_count = df.isna().sum().sum()
total_cells = df.shape[0] * (df.shape[1] - 1)
print(f"NaN ratio: {nan_count / total_cells:.2%}")
# 應該有合理的 NaN 比例（通常 30%-70%）
```

4. **CIGAR 驗證** (詳細檢查 2-3 個 reads):

```bash
# 檢視 debug_cigar.txt
cat output/region_0/debug_cigar.txt
# 手動驗證 SNV 位置與 AltSupport 判斷
```

5. **CpG 座標驗證** (抽樣 5 個 CpG):

```bash
# 從 sites.csv 選擇座標
# 使用 samtools faidx 驗證
samtools faidx hg38.fa chr17:7578000-7578002
# 應輸出 "CG"
```

#### 5.4 性能分析

**預期結果**:

* 單個 region 處理時間: 1-5 秒（取決於 read 數量）
* 總處理時間: 10-40 秒（32 regions, 4 threads）
* Peak RSS: < 2 GB
* 加速比: 2.5-3.5x（相對於單執行緒）

**記錄表格**:

```
| Region ID | Reads | CpGs | Processing Time (ms) | Memory (MB) |
|-----------|-------|------|----------------------|-------------|
| 0         | 45    | 180  | 1200                 | 15          |
| 1         | 52    | 210  | 1400                 | 18          |
| ...       | ...   | ...  | ...                  | ...         |
```

#### 5.5 撰寫驗證報告

**文件** (`docs/dev/phase5_validation_report.md`):

內容應包括:

1. **執行環境**: OS, CPU, 記憶體
2. **測試配置**: 32 SNVs, 4 threads
3. **性能數據**: 總時間, Peak RSS, 加速比
4. **正確性驗證結果**: 各項檢查的結果
5. **發現的問題**: 列出所有問題與解決方法
6. **輸出範例**: 貼上 1-2 個 region 的輸出檔案片段

### 第五階段：交付成果 (Deliverables)

* [x] 32 個 regions 的完整輸出
* [x] 資源使用報告 (`resource_usage.txt`)
* [x] 驗證報告 (`phase5_validation_report.md`)
* [x] 所有驗證清單項目通過

---

## 驗證策略總清單 (Complete Verification Checklist)

執行完第五階段後，使用此清單進行最終確認：

### 功能正確性

* [ ] **BAM Reading**: 是否擷取了視窗內的所有 reads？
  * 方法: 與 `samtools view -c` 比對至少 3 個 regions
  * 容許誤差: ±5% (因過濾條件)

* [ ] **CIGAR**: Read 中的 SNV 位置識別是否正確？
  * 方法: 手動驗證 3 個 reads 的 AltSupport 判斷
  * 預期: 100% 正確

* [ ] **Methylation**: 映射的 CpG 位點在參考序列中是否真的有 'CG'？
  * 方法: 隨機抽樣 5 個 CpG，用 `samtools faidx` 驗證
  * 預期: 100% 為 CG dinucleotide

* [ ] **Matrix**: 空白欄位是否為 NaN？數值是否為 0-1？
  * 方法: 用 Python 檢查 `matrix_raw.csv`
  * 預期: 所有非 NaN 數值在 [0.0, 1.0] 範圍內

### 性能

* [ ] **Parallelism**: 使用 4 執行緒是否比 1 執行緒快？
  * 方法: 分別執行並記錄時間
  * 預期: 加速比 > 2.0x

* [ ] **Memory**: 是否有洩漏？
  * 方法: 使用 Valgrind 或 ASan 在小資料集上執行
  * 預期: 無記憶體洩漏報告

### 輸出格式

* [ ] CSV 檔案包含正確的表頭
* [ ] JSONL 每行都是合法的 JSON
* [ ] 所有 region 目錄結構一致
* [ ] 檔案可被後續模組正確讀取

### 文件

* [ ] `phase1_report.md` 完成
* [ ] `phase2_cigar_verification.md` 完成
* [ ] `phase5_validation_report.md` 完成

---

## 後續工作 (Next Steps)

完成本計畫後，輸出的資料將作為後續模組的輸入：

1. **距離計算模組**: 讀取 `matrix_binary.csv`，計算 read 間的 NHD 距離
2. **聚類模組**: 對距離矩陣進行階層聚類
3. **關聯分析模組**: 結合 `reads.jsonl` 的標籤進行統計檢定

**介面文件**: 需撰寫 `docs/dev/interface_spec.md`，詳細說明輸出格式與後續模組的讀取方式。
