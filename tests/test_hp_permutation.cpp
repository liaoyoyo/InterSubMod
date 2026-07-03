#include <gtest/gtest.h>

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "core/DataStructs.hpp"
#include "utils/HpPermutation.hpp"

using namespace InterSubMod;

namespace {

std::vector<ReadInfo> make_reads(const std::vector<std::string>& tags) {
    std::vector<ReadInfo> reads;
    int id = 0;
    for (const auto& t : tags) {
        ReadInfo r{};
        r.read_id = id++;
        r.hp_tag = t;
        r.hp_tag_raw = t;  // mirror raw = demoted for the simple cases
        reads.push_back(r);
    }
    return reads;
}

std::map<std::string, int> tag_counts(const std::vector<ReadInfo>& reads) {
    std::map<std::string, int> m;
    for (const auto& r : reads) {
        m[r.hp_tag]++;
    }
    return m;
}

}  // namespace

// seed <= 0 must be a no-op: the default pipeline path (--permute-hp-seed 0) is byte-identical.
TEST(HpPermutationTest, SeedZeroIsNoOp) {
    auto reads = make_reads({"1", "2", "1-1", "2-1", "3"});
    auto before = reads;
    permute_hp_tags_per_region(reads, 0, "chr1", 12345);
    for (std::size_t i = 0; i < reads.size(); ++i) {
        EXPECT_EQ(reads[i].hp_tag, before[i].hp_tag);
        EXPECT_EQ(reads[i].hp_tag_raw, before[i].hp_tag_raw);
    }
}

TEST(HpPermutationTest, NegativeSeedIsNoOp) {
    auto reads = make_reads({"1", "2", "1-1"});
    auto before = reads;
    permute_hp_tags_per_region(reads, -5, "chr1", 100);
    for (std::size_t i = 0; i < reads.size(); ++i) {
        EXPECT_EQ(reads[i].hp_tag, before[i].hp_tag);
    }
}

TEST(HpPermutationTest, FewerThanTwoReadsIsNoOp) {
    auto reads = make_reads({"1-1"});
    permute_hp_tags_per_region(reads, 7, "chr1", 100);
    EXPECT_EQ(reads[0].hp_tag, "1-1");
}

// Permutation preserves the per-region HP marginal (same multiset of tags).
TEST(HpPermutationTest, PreservesMarginal) {
    auto reads = make_reads({"1", "1", "2", "1-1", "1-1", "1-1", "2-1", "3", "0", "0"});
    auto before = tag_counts(reads);
    permute_hp_tags_per_region(reads, 42, "chr7", 55555);
    auto after = tag_counts(reads);
    EXPECT_EQ(before, after);
}

// Same (seed, chrom, pos) yields the same permutation (reproducible run).
TEST(HpPermutationTest, Reproducible) {
    auto a = make_reads({"1", "2", "1-1", "2-1", "3", "0", "1", "2"});
    auto b = a;
    permute_hp_tags_per_region(a, 123, "chr3", 999);
    permute_hp_tags_per_region(b, 123, "chr3", 999);
    for (std::size_t i = 0; i < a.size(); ++i) {
        EXPECT_EQ(a[i].hp_tag, b[i].hp_tag);
        EXPECT_EQ(a[i].hp_tag_raw, b[i].hp_tag_raw);
    }
}

// hp_tag and hp_tag_raw are shuffled as a unit: every read's post-permutation pair must be one of
// the original (hp_tag, hp_tag_raw) pairs.
TEST(HpPermutationTest, TagPairStaysUnited) {
    std::vector<ReadInfo> reads;
    const std::vector<std::pair<std::string, std::string>> pairs = {
        {"1", "1-1"}, {"2", "2-1"}, {"0", "3"}, {"1", "1"}, {"2", "2"}};
    for (const auto& p : pairs) {
        ReadInfo r{};
        r.hp_tag = p.first;
        r.hp_tag_raw = p.second;
        reads.push_back(r);
    }
    permute_hp_tags_per_region(reads, 55, "chr9", 42);
    std::set<std::pair<std::string, std::string>> orig(pairs.begin(), pairs.end());
    for (const auto& r : reads) {
        EXPECT_TRUE(orig.count({r.hp_tag, r.hp_tag_raw}) > 0);
    }
}
