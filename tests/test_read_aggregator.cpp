/**
 * @file test_read_aggregator.cpp
 * @brief Unit tests for ReadAggregator and RegionBounds.
 *
 * Tests cover:
 * 1. RegionBounds::compute_region_bounds() — normal, clamped start/end, unknown length
 * 2. RegionBounds::is_valid() and width()
 * 3. ReadAggregator::aggregate() — empty input, low-quality filter, secondary skip
 */

#include <gtest/gtest.h>

#include <htslib/sam.h>

#include "core/ReadAggregator.hpp"
#include "core/RegionBounds.hpp"
#include "core/SomaticSnv.hpp"

using namespace InterSubMod;

// ============================================================================
// RegionBounds tests
// ============================================================================

TEST(RegionBoundsTest, ComputeRegionBounds_Normal) {
    auto b = compute_region_bounds("chr17", 7577000, 5000, 80000000);
    EXPECT_EQ(b.chr, "chr17");
    EXPECT_EQ(b.start, 7572000);
    EXPECT_EQ(b.end, 7582000);
    EXPECT_TRUE(b.is_valid());
    EXPECT_EQ(b.width(), 10000);
}

TEST(RegionBoundsTest, ComputeRegionBounds_ClampStart) {
    // SNV near chromosome start — start must not go below 1
    auto b = compute_region_bounds("chr1", 100, 5000, 100000000);
    EXPECT_EQ(b.start, 1);
    EXPECT_EQ(b.end, 5100);
    EXPECT_TRUE(b.is_valid());
}

TEST(RegionBoundsTest, ComputeRegionBounds_ClampEnd) {
    // SNV near chromosome end — end must not exceed chr_length
    // snv_pos + window_size = 99999000 + 5000 = 100004000 > chr_length → clamped
    auto b = compute_region_bounds("chr1", 99999000, 5000, 100000000);
    EXPECT_EQ(b.start, 99994000);
    EXPECT_EQ(b.end, 100000000);
    EXPECT_TRUE(b.is_valid());
}

TEST(RegionBoundsTest, ComputeRegionBounds_UnknownLength) {
    // chr_length == 0 means unknown: end is not clamped
    auto b = compute_region_bounds("chrX", 50000, 2500, 0);
    EXPECT_EQ(b.start, 47500);
    EXPECT_EQ(b.end, 52500);
    EXPECT_TRUE(b.is_valid());
}

TEST(RegionBoundsTest, IsValid_ZeroStartInvalid) {
    RegionBounds b;
    b.chr   = "chr1";
    b.start = 0;
    b.end   = 1000;
    EXPECT_FALSE(b.is_valid());
}

TEST(RegionBoundsTest, IsValid_EndLeStartInvalid) {
    RegionBounds b;
    b.chr   = "chr1";
    b.start = 1000;
    b.end   = 500;
    EXPECT_FALSE(b.is_valid());
}

// ============================================================================
// ReadAggregator tests
// ============================================================================

TEST(ReadAggregatorTest, EmptyInput_ReturnsEmptyResult) {
    ReadAggregator agg;
    SomaticSnv     snv{};
    auto result = agg.aggregate({}, snv, "", 0);
    EXPECT_EQ(result.matrix_builder.num_reads(), 0);
    EXPECT_EQ(result.num_forward_reads, 0);
    EXPECT_EQ(result.num_reverse_reads, 0);
    EXPECT_TRUE(result.filtered_reads.empty());
}

TEST(ReadAggregatorTest, LowMapqRead_DefaultConfig_SilentlyFiltered) {
    // A bare bam_init1() has qual=0, l_qseq=0 — fails MAPQ and read-length filters.
    // With collect_filtered_reads=false (default), no crash occurs.
    ReadAggregator agg;
    SomaticSnv     snv{};

    bam1_t* b = bam_init1();
    auto result = agg.aggregate({b}, snv, "", 0);

    EXPECT_EQ(result.matrix_builder.num_reads(), 0);
    EXPECT_TRUE(result.filtered_reads.empty());

    bam_destroy1(b);
}

TEST(ReadAggregatorTest, SecondaryRead_NoFilterMode_Skipped) {
    // In no_filter_output mode, secondary/supplementary reads must still be skipped
    // to avoid undefined behaviour in methylation parsing.
    ReadAggregateConfig cfg;
    cfg.no_filter_output = true;
    ReadAggregator agg(cfg);
    SomaticSnv     snv{};

    bam1_t* b   = bam_init1();
    b->core.flag = BAM_FSECONDARY;

    auto result = agg.aggregate({b}, snv, "", 0);

    EXPECT_EQ(result.matrix_builder.num_reads(), 0);

    bam_destroy1(b);
}

TEST(ReadAggregatorTest, SupplementaryRead_NoFilterMode_Skipped) {
    ReadAggregateConfig cfg;
    cfg.no_filter_output = true;
    ReadAggregator agg(cfg);
    SomaticSnv     snv{};

    bam1_t* b   = bam_init1();
    b->core.flag = BAM_FSUPPLEMENTARY;

    auto result = agg.aggregate({b}, snv, "", 0);
    EXPECT_EQ(result.matrix_builder.num_reads(), 0);

    bam_destroy1(b);
}

TEST(ReadAggregatorTest, MultipleFilteredReads_AllRejected) {
    // Three reads with qual=0; all filtered, matrix stays empty.
    ReadAggregator agg;
    SomaticSnv     snv{};

    std::vector<bam1_t*> reads;
    for (int i = 0; i < 3; ++i) {
        reads.push_back(bam_init1());
    }

    auto result = agg.aggregate(reads, snv, "", 0);
    EXPECT_EQ(result.matrix_builder.num_reads(), 0);
    EXPECT_EQ(result.num_forward_reads, 0);
    EXPECT_EQ(result.num_reverse_reads, 0);

    for (auto* b : reads) bam_destroy1(b);
}
