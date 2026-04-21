/**
 * @file test_read_parser.cpp
 * @brief Unit tests for ReadParser HP tag handling, including the
 *        germline_hp_only fallback that demotes somatic phase-block
 *        tags (HP:i:11/21/33) to unphased while preserving raw values.
 */

#include <gtest/gtest.h>

#include <htslib/sam.h>

#include <cstring>
#include <string>
#include <vector>

#include "core/ReadParser.hpp"
#include "core/SomaticSnv.hpp"

using namespace InterSubMod;

namespace {

// Build a minimal BAM record carrying a single HP integer aux tag.
// The record is otherwise empty (no sequence/qual/cigar); parse()'s
// alt-support branch short-circuits when anchor_snv.pos == 0 so the
// skeletal record is sufficient for HP-tag assertions.
bam1_t* make_bam_with_hp_int(const char* qname, int32_t hp_value) {
    bam1_t* b = bam_init1();
    if (bam_set1(b, std::strlen(qname), qname,
                 0,           // flag
                 0,           // tid
                 100,         // pos (0-based)
                 60,          // mapq
                 0, nullptr,  // cigar
                 -1, -1, 0,   // mtid, mpos, isize
                 0, nullptr, nullptr,  // seq, qual
                 0) < 0) {
        bam_destroy1(b);
        return nullptr;
    }
    bam_aux_append(b, "HP", 'i', sizeof(hp_value),
                   reinterpret_cast<uint8_t*>(&hp_value));
    return b;
}

bam1_t* make_bam_with_hp_str(const char* qname, const char* hp_value) {
    bam1_t* b = bam_init1();
    if (bam_set1(b, std::strlen(qname), qname, 0, 0, 100, 60,
                 0, nullptr, -1, -1, 0, 0, nullptr, nullptr, 0) < 0) {
        bam_destroy1(b);
        return nullptr;
    }
    // 'Z' aux: length includes the trailing NUL.
    bam_aux_append(b, "HP", 'Z',
                   static_cast<int>(std::strlen(hp_value) + 1),
                   reinterpret_cast<const uint8_t*>(hp_value));
    return b;
}

bam1_t* make_bam_without_hp(const char* qname) {
    bam1_t* b = bam_init1();
    if (bam_set1(b, std::strlen(qname), qname, 0, 0, 100, 60,
                 0, nullptr, -1, -1, 0, 0, nullptr, nullptr, 0) < 0) {
        bam_destroy1(b);
        return nullptr;
    }
    return b;
}

// Parse a record with a given germline_hp_only flag and return the populated ReadInfo.
ReadInfo parse_single(bam1_t* b, bool germline_hp_only) {
    ReadFilterConfig cfg;
    cfg.germline_hp_only = germline_hp_only;
    ReadParser parser(cfg);
    SomaticSnv snv{};  // pos == 0 short-circuits alt-support evaluation
    return parser.parse(b, /*read_id=*/0, /*is_tumor=*/true, snv,
                        /*ref_seq=*/"", /*ref_start_pos=*/0);
}

}  // namespace

// ============================================================================
// HP tag parsing — germline values (1/2) and unphased (0 / missing)
// ============================================================================

TEST(ReadParserHPTagTest, GermlineHP1_FlagOff_PassesThrough) {
    bam1_t* b = make_bam_with_hp_int("r1", 1);
    ReadInfo info = parse_single(b, /*germline_hp_only=*/false);
    EXPECT_EQ(info.hp_tag_raw, "1");
    EXPECT_EQ(info.hp_tag, "1");
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, GermlineHP1_FlagOn_PassesThrough) {
    bam1_t* b = make_bam_with_hp_int("r1", 1);
    ReadInfo info = parse_single(b, /*germline_hp_only=*/true);
    EXPECT_EQ(info.hp_tag_raw, "1");
    EXPECT_EQ(info.hp_tag, "1");
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, GermlineHP2_FlagOn_PassesThrough) {
    bam1_t* b = make_bam_with_hp_int("r1", 2);
    ReadInfo info = parse_single(b, /*germline_hp_only=*/true);
    EXPECT_EQ(info.hp_tag_raw, "2");
    EXPECT_EQ(info.hp_tag, "2");
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, NoHPTag_DefaultsToZero_BothFlags) {
    bam1_t* b = make_bam_without_hp("r1");
    ReadInfo info_off = parse_single(b, false);
    ReadInfo info_on = parse_single(b, true);
    EXPECT_EQ(info_off.hp_tag_raw, "0");
    EXPECT_EQ(info_off.hp_tag, "0");
    EXPECT_EQ(info_on.hp_tag_raw, "0");
    EXPECT_EQ(info_on.hp_tag, "0");
    bam_destroy1(b);
}

// ============================================================================
// HP tag parsing — somatic phase-block tags (11/21/33)
// The core fix: when germline_hp_only is set, these are demoted to "0"
// downstream, but the raw value must always be preserved for audit.
// ============================================================================

TEST(ReadParserHPTagTest, Somatic11_FlagOff_KeepsOneDashOne) {
    bam1_t* b = make_bam_with_hp_int("r1", 11);
    ReadInfo info = parse_single(b, false);
    EXPECT_EQ(info.hp_tag_raw, "1-1");
    EXPECT_EQ(info.hp_tag, "1-1");
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, Somatic11_FlagOn_DemotedToUnphased) {
    bam1_t* b = make_bam_with_hp_int("r1", 11);
    ReadInfo info = parse_single(b, true);
    EXPECT_EQ(info.hp_tag_raw, "1-1");  // audit value retained
    EXPECT_EQ(info.hp_tag, "0");        // demoted for downstream features
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, Somatic21_FlagOn_DemotedToUnphased) {
    bam1_t* b = make_bam_with_hp_int("r1", 21);
    ReadInfo info = parse_single(b, true);
    EXPECT_EQ(info.hp_tag_raw, "2-1");
    EXPECT_EQ(info.hp_tag, "0");
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, Somatic33_FlagOn_DemotedToUnphased) {
    bam1_t* b = make_bam_with_hp_int("r1", 33);
    ReadInfo info = parse_single(b, true);
    EXPECT_EQ(info.hp_tag_raw, "3");
    EXPECT_EQ(info.hp_tag, "0");
    bam_destroy1(b);
}

// ============================================================================
// String-format HP tags (longphase-s output) — same demotion rules apply.
// ============================================================================

TEST(ReadParserHPTagTest, StringHP1Dash1_FlagOn_Demoted) {
    bam1_t* b = make_bam_with_hp_str("r1", "1-1");
    ReadInfo info = parse_single(b, true);
    EXPECT_EQ(info.hp_tag_raw, "1-1");
    EXPECT_EQ(info.hp_tag, "0");
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, StringHP3_FlagOn_Demoted) {
    bam1_t* b = make_bam_with_hp_str("r1", "3");
    ReadInfo info = parse_single(b, true);
    EXPECT_EQ(info.hp_tag_raw, "3");
    EXPECT_EQ(info.hp_tag, "0");
    bam_destroy1(b);
}

TEST(ReadParserHPTagTest, StringHP1_FlagOn_NotDemoted) {
    bam1_t* b = make_bam_with_hp_str("r1", "1");
    ReadInfo info = parse_single(b, true);
    EXPECT_EQ(info.hp_tag_raw, "1");
    EXPECT_EQ(info.hp_tag, "1");
    bam_destroy1(b);
}

// ============================================================================
// Aggregated-distribution smoke test mirroring the plan's Fixture A.
// 20 reads with HP values [1x5, 2x5, 11x3, 21x3, 33x2, 0x2].
// flag=off: six distinct hp_tag categories.
// flag=on:  three categories (1, 2, 0) — somatic tags folded into 0.
// ============================================================================

TEST(ReadParserHPTagTest, Fixture20Reads_DistributionMatchesExpectation) {
    std::vector<int32_t> hp_seq = {1, 1, 1, 1, 1,
                                    2, 2, 2, 2, 2,
                                    11, 11, 11,
                                    21, 21, 21,
                                    33, 33};
    std::vector<bam1_t*> unphased = {make_bam_without_hp("u1"),
                                      make_bam_without_hp("u2")};

    auto tally = [](const std::vector<ReadInfo>& infos, const std::string& key) {
        int n = 0;
        for (const auto& i : infos) {
            if (i.hp_tag == key) ++n;
        }
        return n;
    };
    auto tally_raw = [](const std::vector<ReadInfo>& infos, const std::string& key) {
        int n = 0;
        for (const auto& i : infos) {
            if (i.hp_tag_raw == key) ++n;
        }
        return n;
    };

    std::vector<bam1_t*> all;
    for (int32_t v : hp_seq) {
        all.push_back(make_bam_with_hp_int("r", v));
    }
    for (auto* b : unphased) all.push_back(b);

    // flag=off: pass-through behaviour
    {
        std::vector<ReadInfo> infos;
        infos.reserve(all.size());
        for (auto* b : all) infos.push_back(parse_single(b, false));
        EXPECT_EQ(tally(infos, "1"), 5);
        EXPECT_EQ(tally(infos, "2"), 5);
        EXPECT_EQ(tally(infos, "1-1"), 3);
        EXPECT_EQ(tally(infos, "2-1"), 3);
        EXPECT_EQ(tally(infos, "3"), 2);
        EXPECT_EQ(tally(infos, "0"), 2);
    }

    // flag=on: somatic tags fold into "0"; raw values still recoverable
    {
        std::vector<ReadInfo> infos;
        infos.reserve(all.size());
        for (auto* b : all) infos.push_back(parse_single(b, true));
        EXPECT_EQ(tally(infos, "1"), 5);
        EXPECT_EQ(tally(infos, "2"), 5);
        EXPECT_EQ(tally(infos, "0"), 2 + 3 + 3 + 2);  // 2 unphased + 11/21/33
        EXPECT_EQ(tally(infos, "1-1"), 0);
        EXPECT_EQ(tally(infos, "2-1"), 0);
        EXPECT_EQ(tally(infos, "3"), 0);
        // Raw values remain untouched — critical for NHP_Somatic audit counters
        EXPECT_EQ(tally_raw(infos, "1-1"), 3);
        EXPECT_EQ(tally_raw(infos, "2-1"), 3);
        EXPECT_EQ(tally_raw(infos, "3"), 2);
    }

    for (auto* b : all) bam_destroy1(b);
}
