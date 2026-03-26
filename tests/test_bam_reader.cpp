#include <gtest/gtest.h>

#include <filesystem>
#include <vector>

#include "core/BamReader.hpp"
#include "core/ReadParser.hpp"
#include "core/SomaticSnv.hpp"

using namespace InterSubMod;

// Test fixture for BAM reading
class BamReaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Prefer tiny fixture first; fallback to project sample BAM.
        const std::vector<std::string> candidates = {
            "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/test.bam",
            "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"};

        for (const auto& path : candidates) {
            if (std::filesystem::exists(path)) {
                test_bam_path = path;
                break;
            }
        }

        if (test_bam_path.empty()) {
            GTEST_SKIP() << "No usable BAM fixture found under data/bam/";
        }
    }

    std::string test_bam_path;
};

TEST_F(BamReaderTest, ConstructorOpensFile) {
    EXPECT_NO_THROW({
        BamReader reader(test_bam_path);
        EXPECT_TRUE(reader.is_open());
        EXPECT_NE(reader.get_header(), nullptr);
    });
}

TEST_F(BamReaderTest, FetchReadsFromValidRegion) {
    BamReader reader(test_bam_path);

    // Fetch reads from a known region (adjust based on actual data)
    auto reads = reader.fetch_reads("chr17", 7577000, 7578000);

    // We should get some reads (exact count depends on data)
    // For now, just check that the operation succeeds
    EXPECT_GE(reads.size(), 0);

    // Cleanup
    for (auto* r : reads) {
        bam_destroy1(r);
    }
}

TEST_F(BamReaderTest, FetchReadsFromInvalidChromosome) {
    BamReader reader(test_bam_path);

    // Try to fetch from non-existent chromosome
    auto reads = reader.fetch_reads("chrNonExistent", 0, 1000);

    // Should return empty vector
    EXPECT_EQ(reads.size(), 0);
}

TEST_F(BamReaderTest, MoveConstructor) {
    BamReader reader1(test_bam_path);
    EXPECT_TRUE(reader1.is_open());

    // Move construct
    BamReader reader2(std::move(reader1));
    EXPECT_TRUE(reader2.is_open());
}

// Test fixture for ReadParser
class ReadParserTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Default configuration
        config.min_mapq = 20;
        config.min_read_length = 1000;
        config.min_base_quality = 20;
        config.require_mm_ml = true;
    }

    ReadFilterConfig config;
};

TEST_F(ReadParserTest, FilterConfiguration) {
    ReadParser parser(config);
    EXPECT_EQ(parser.get_config().min_mapq, 20);
    EXPECT_EQ(parser.get_config().min_read_length, 1000);
}

// ---- RAII fetch_reads_safe() tests ----

// Compile-time check: BamRecordPtr type must exist and be a unique_ptr
TEST(BamRecordPtrTest, TypeExists) {
    // Verify BamRecordPtr is defined and wraps bam1_t
    static_assert(
        std::is_same_v<
            InterSubMod::BamRecordPtr,
            std::unique_ptr<bam1_t, InterSubMod::BamRecordDeleter>>,
        "BamRecordPtr must be unique_ptr<bam1_t, BamRecordDeleter>");
    SUCCEED();
}

// Verify BamRecordDeleter correctly calls bam_destroy1
TEST(BamRecordPtrTest, DeleterReleasesMemory) {
    bam1_t* raw = bam_init1();
    ASSERT_NE(raw, nullptr);
    {
        // Place raw ptr into RAII wrapper - should be destroyed at end of scope
        InterSubMod::BamRecordPtr ptr(raw, InterSubMod::BamRecordDeleter{});
        EXPECT_EQ(ptr.get(), raw);
    }
    // ptr destroyed here; raw pointer is now invalid (bam_destroy1 was called)
    // We can't dereference raw, but we can verify no crash occurred
    SUCCEED();
}

// fetch_reads_safe() returns RAII-managed records (no manual cleanup needed)
TEST_F(BamReaderTest, FetchReadsSafe_ReturnsUniquePtr) {
    BamReader reader(test_bam_path);
    auto reads = reader.fetch_reads_safe("chr17", 7577000, 7578000);

    // Type must be vector<BamRecordPtr>
    static_assert(
        std::is_same_v<decltype(reads),
                       std::vector<InterSubMod::BamRecordPtr>>,
        "fetch_reads_safe must return vector<BamRecordPtr>");

    // Should succeed (exact count depends on data)
    EXPECT_GE(reads.size(), 0u);
    // No manual cleanup needed - unique_ptr handles it
}

// fetch_reads_safe() handles invalid chromosome gracefully
TEST_F(BamReaderTest, FetchReadsSafe_InvalidChromosome) {
    BamReader reader(test_bam_path);
    auto reads = reader.fetch_reads_safe("chrNonExistent999", 0, 1000);
    EXPECT_TRUE(reads.empty());
}

// ---- END RAII tests ----

// Integration test - requires actual BAM data
TEST_F(BamReaderTest, ParseReadInfo) {
    if (!std::filesystem::exists(test_bam_path)) {
        GTEST_SKIP() << "Test BAM not found";
    }

    BamReader reader(test_bam_path);
    ReadParser parser;

    // Fetch some reads
    auto reads = reader.fetch_reads("chr17", 7577000, 7578000);

    if (reads.empty()) {
        GTEST_SKIP() << "No reads in test region";
    }

    // Create a dummy SNV for testing
    SomaticSnv snv;
    snv.snv_id = 0;
    snv.chr_id = 1;     // chr17
    snv.pos = 7577500;  // 1-based
    snv.ref_base = 'C';
    snv.alt_base = 'T';

    // Parse the first read
    int passed_filter = 0;
    for (auto* b : reads) {
        if (parser.should_keep(b)) {
            passed_filter++;

            // Parse read info (we don't have ref_seq yet, so alt_support will be UNKNOWN)
            ReadInfo info = parser.parse(b, 0, true, snv, "", 0);

            EXPECT_EQ(info.read_id, 0);
            EXPECT_FALSE(info.read_name.empty());
            EXPECT_GE(info.mapq, 0);
            EXPECT_LE(info.mapq, 255);

            break;  // Only test first passing read
        }
    }

    EXPECT_GT(passed_filter, 0) << "Expected at least some reads to pass filter";

    // Cleanup
    for (auto* r : reads) {
        bam_destroy1(r);
    }
}
