/**
 * @file test_methylation_parser.cpp
 * @brief Characterization tests for MethylationParser MM/ML -> read x CpG calls.
 *
 * This is the data-ingestion path that feeds EVERY downstream methylation
 * number (Delta-beta, distance matrix, clustering). It was previously
 * UNTESTED (audit 2026-06-10, finding TESTS-1). These tests pin the
 * spec-correct behavior derived independently from the SAM MM/ML
 * specification, NOT from the implementation:
 *
 *   - MM "C+m?,d0,d1,..." delta decoding (skip-counts between modified C's)
 *   - ML[i]/255 -> probability
 *   - a modified C is kept ONLY in CpG context (ref C followed by ref G)
 *   - 1-based genomic coordinate reporting (MethylCall::ref_pos doc contract)
 *   - reverse-strand reads collapse onto the SAME forward-CpG coordinate
 *     as the forward strand (two-strand CpG collapse) — the most error-prone
 *     branch (backward iteration over BAM SEQ).
 *
 * If any of these FAIL, that is a real bug in a load-bearing parser, not a
 * stale test — investigate before "fixing" the assertion.
 */

#include <gtest/gtest.h>

#include <htslib/sam.h>

#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "core/MethylationParser.hpp"

using namespace InterSubMod;

namespace {

// Build an in-memory aligned BAM record: full-match (LxM) read at `pos`
// (0-based), with the given BAM SEQ, MM (Z) and ML (B,C) tags.
// `flag` carries strand (BAM_FREVERSE for reverse). BAM SEQ is always stored
// in reference-forward orientation, so for a perfect match seq == ref window.
bam1_t* make_meth_bam(const char* qname, uint16_t flag, int32_t pos, const std::string& seq,
                      const char* mm, const std::vector<uint8_t>& ml) {
    bam1_t* b = bam_init1();
    uint32_t cigar = bam_cigar_gen(static_cast<uint32_t>(seq.size()), BAM_CMATCH);
    if (bam_set1(b, std::strlen(qname), qname, flag, /*tid=*/0, pos, /*mapq=*/60,
                 /*n_cigar=*/1, &cigar, /*mtid=*/-1, /*mpos=*/-1, /*isize=*/0,
                 seq.size(), seq.c_str(), /*qual=*/nullptr, /*l_aux=*/0) < 0) {
        bam_destroy1(b);
        return nullptr;
    }
    bam_aux_append(b, "MM", 'Z', static_cast<int>(std::strlen(mm) + 1),
                   reinterpret_cast<const uint8_t*>(mm));
    // ML is a B-array of uint8 ('C'); bam_aux_update_array appends if absent.
    bam_aux_update_array(b, "ML", 'C', static_cast<uint32_t>(ml.size()),
                         const_cast<uint8_t*>(ml.data()));
    return b;
}

}  // namespace

// ============================================================================
// Forward strand: delta decode + CpG filter + probability from ML
// ============================================================================

TEST(MethylationParserTest, ForwardStrand_DecodesDeltas_FiltersNonCpg_ProbFromML) {
    // ref window @ 0-based 1000:  A C G A A C G T T C
    //                  offset:    0 1 2 3 4 5 6 7 8 9
    // CpGs (C then G): offset 1 and offset 5. offset 9 is a C but at the window
    // edge (no offset 10) so it is NOT a valid CpG -> must be dropped.
    const std::string ref = "ACGAACGTTC";
    // Read == ref (perfect forward match). C's are at seq idx 1,5,9.
    // MM "C+m?,0,0,0" marks all three C's modified; ML gives their probs.
    std::vector<uint8_t> ml = {230, 130, 255};
    bam1_t* b = make_meth_bam("fwd", /*flag=*/0, /*pos=*/1000, ref, "C+m?,0,0,0;", ml);
    ASSERT_NE(b, nullptr);

    MethylationParser parser;
    auto calls = parser.parse_read(b, ref, /*ref_start_pos=*/1000);

    // Two kept (offsets 1,5 are CpG); the edge C at offset 9 is dropped.
    ASSERT_EQ(calls.size(), 2u);
    // 1-based coordinate of the forward C = ref_start + offset + 1.
    EXPECT_EQ(calls[0].ref_pos, 1002);  // offset 1 -> 1000+1+1
    EXPECT_NEAR(calls[0].probability, 230.0f / 255.0f, 1e-5f);
    EXPECT_EQ(calls[1].ref_pos, 1006);  // offset 5 -> 1000+5+1
    EXPECT_NEAR(calls[1].probability, 130.0f / 255.0f, 1e-5f);

    bam_destroy1(b);
}

// ============================================================================
// Reverse strand collapses onto the SAME CpG coordinate as the forward strand.
// This pins the two-strand CpG collapse + the backward-iteration branch.
// ============================================================================

TEST(MethylationParserTest, ReverseStrand_CollapsesToForwardCpgCoordinate) {
    // ref window @ 0-based 1000:  A A C G T T A C G T
    //                  offset:    0 1 2 3 4 5 6 7 8 9
    // CpGs (C then G): offset 2 and offset 7.
    const std::string ref = "AACGTTACGT";

    // --- Forward read over both CpGs ---
    bam1_t* fwd = make_meth_bam("fwd", /*flag=*/0, 1000, ref, "C+m?,0,0;", {200, 100});
    ASSERT_NE(fwd, nullptr);
    MethylationParser parser;
    auto fcalls = parser.parse_read(fwd, ref, 1000);
    ASSERT_EQ(fcalls.size(), 2u);
    // forward C at offset 2 -> 1003 ; offset 7 -> 1008
    EXPECT_EQ(fcalls[0].ref_pos, 1003);
    EXPECT_EQ(fcalls[1].ref_pos, 1008);

    // --- Reverse read over the same CpGs (BAM SEQ identical, FREVERSE set) ---
    bam1_t* rev = make_meth_bam("rev", /*flag=*/BAM_FREVERSE, 1000, ref, "C+m?,0,0;", {200, 100});
    ASSERT_NE(rev, nullptr);
    auto rcalls = parser.parse_read(rev, ref, 1000);
    ASSERT_EQ(rcalls.size(), 2u);

    // The KEY property: reverse calls land on the SAME forward-CpG coordinates
    // {1003, 1008} as the forward strand (two-strand collapse), so both strands
    // of a CpG aggregate to one column in the read x CpG matrix.
    std::map<int32_t, float> rev_by_pos;
    for (const auto& c : rcalls) rev_by_pos[c.ref_pos] = c.probability;
    ASSERT_EQ(rev_by_pos.count(1003), 1u);
    ASSERT_EQ(rev_by_pos.count(1008), 1u);
    // Probability mapping on the reverse strand (MM order is original-read 5'->3',
    // iterated backward over BAM SEQ): pos 1008 <- ML[0]=200, pos 1003 <- ML[1]=100.
    EXPECT_NEAR(rev_by_pos[1008], 200.0f / 255.0f, 1e-5f);
    EXPECT_NEAR(rev_by_pos[1003], 100.0f / 255.0f, 1e-5f);

    bam_destroy1(fwd);
    bam_destroy1(rev);
}

// ============================================================================
// No MM tag -> no calls (documented: may be empty).
// ============================================================================

TEST(MethylationParserTest, NoMmTag_ReturnsEmpty) {
    const std::string ref = "ACGACG";
    bam1_t* b = bam_init1();
    uint32_t cigar = bam_cigar_gen(static_cast<uint32_t>(ref.size()), BAM_CMATCH);
    ASSERT_GE(bam_set1(b, 3, "nom", 0, 0, 1000, 60, 1, &cigar, -1, -1, 0, ref.size(), ref.c_str(),
                       nullptr, 0),
              0);
    // No MM/ML appended.
    MethylationParser parser;
    auto calls = parser.parse_read(b, ref, 1000);
    EXPECT_TRUE(calls.empty());
    bam_destroy1(b);
}
