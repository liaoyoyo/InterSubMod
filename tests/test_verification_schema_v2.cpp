#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "core/RegionProcessor.hpp"

using namespace InterSubMod;

namespace {

struct VerificationCase {
    VerificationV2Input input;
    std::string expected_class;
    std::string expected_path;
    bool expected_significant;
};

VerificationV2Input base_input(const std::string& legacy = "Weak") {
    VerificationV2Input input;
    input.legacy_class = legacy;
    input.pairwise_mean_distance = 0.25;
    input.per_cpg_epipoly_hp1 = 0.30;
    return input;
}

}  // namespace

TEST(VerificationSchemaV2Test, AuthoritativeElevenClassTruthTable) {
    std::vector<VerificationCase> cases;

    auto input = base_input("Strong");
    cases.push_back({input, "Strong_Bidirectional", "BIDIRECTIONAL", true});
    input = base_input("Subclone");
    cases.push_back({input, "ClusterFirstOnly", "CLUSTER_FIRST_ONLY", true});
    input = base_input();
    input.potential_loh = true;
    input.hp_auc_struct = true;
    cases.push_back({input, "LOH-Structure", "LOH_STRUCTURE", true});
    input = base_input();
    input.within_hp = true;
    cases.push_back({input, "MultiGroupNoLabel", "WITHIN_HP_MULTIGROUP", true});
    input = base_input();
    input.dbeta_sig = true;
    cases.push_back({input, "LabelShift", "LABEL_SHIFT", true});
    input = base_input();
    input.clean_location_permanova = true;
    cases.push_back({input, "PermanovaLocation", "PERMANOVA_LOCATION", true});
    input = base_input();
    input.hp_auc_struct = true;
    cases.push_back({input, "StructureNoLabel", "HP_AUC_STRUCTURE_NO_LABEL", true});
    input = base_input();
    input.dispersion_structure = true;
    input.dispersion_warning = true;
    cases.push_back({input, "DispersionStructure", "DISPERSION_STRUCTURE", false});
    input = base_input("Noise");
    input.pairwise_mean_distance = 0.10;
    cases.push_back({input, "Noise_Uniform", "NOISE_UNIFORM", false});
    input = base_input("Noise");
    input.pairwise_mean_distance = 0.40;
    cases.push_back({input, "Noise_Chaotic", "NOISE_CHAOTIC", false});
    input = base_input("Noise");
    cases.push_back({input, "Noise_Uncorrelated", "NOISE_UNCORRELATED", false});

    ASSERT_EQ(cases.size(), 11U);
    for (const auto& fixture : cases) {
        const auto decision = classify_verification_v2(fixture.input);
        EXPECT_EQ(decision.schema_version, 2);
        EXPECT_EQ(decision.verification_class, fixture.expected_class);
        EXPECT_EQ(decision.evidence_path, fixture.expected_path);
        EXPECT_EQ(decision.significant, fixture.expected_significant);
        EXPECT_EQ(decision.verification_class_legacy, fixture.input.legacy_class);
        EXPECT_EQ(decision.label_first_support,
                  fixture.input.legacy_class == "Strong" || fixture.input.legacy_class == "Weak");
        EXPECT_EQ(decision.cluster_first_support,
                  fixture.input.legacy_class == "Strong" || fixture.input.legacy_class == "Subclone");
        EXPECT_EQ(decision.within_hp_support, fixture.input.within_hp);
        EXPECT_EQ(decision.dispersion_warning, fixture.input.dispersion_warning);
        EXPECT_EQ(decision.evidence_derivation, "LIVE");
        const std::string expected_v1 =
            fixture.expected_class == "Strong_Bidirectional" || fixture.expected_class == "ClusterFirstOnly"
                ? "Strong"
                : fixture.expected_class;
        EXPECT_EQ(decision.verification_class_v1_deprecated, expected_v1);
    }
}

TEST(VerificationSchemaV2Test, OrderedPrecedenceIsStableUnderCollisions) {
    auto input = base_input("Subclone");
    input.potential_loh = true;
    input.hp_auc_struct = true;
    input.within_hp = true;
    input.dbeta_sig = true;
    input.clean_location_permanova = true;
    input.dispersion_structure = true;
    EXPECT_EQ(classify_verification_v2(input).verification_class, "ClusterFirstOnly");

    input = base_input("Weak");
    input.potential_loh = true;
    input.hp_auc_struct = true;
    input.within_hp = true;
    input.dbeta_sig = true;
    EXPECT_EQ(classify_verification_v2(input).verification_class, "LOH-Structure");

    input = base_input("Weak");
    input.within_hp = true;
    input.dbeta_sig = true;
    input.clean_location_permanova = true;
    EXPECT_EQ(classify_verification_v2(input).verification_class, "MultiGroupNoLabel");

    input = base_input("Weak");
    input.dbeta_sig = true;
    input.clean_location_permanova = true;
    input.hp_auc_struct = true;
    input.dispersion_structure = true;
    EXPECT_EQ(classify_verification_v2(input).verification_class, "LabelShift");

    input = base_input("Weak");
    input.clean_location_permanova = true;
    input.hp_auc_struct = true;
    EXPECT_EQ(classify_verification_v2(input).verification_class, "PermanovaLocation");
}

TEST(VerificationSchemaV2Test, LegacyLohSubtypeAndDeprecatedAliasSourceAreExact) {
    EXPECT_EQ(determine_loh_subtype_legacy_vc(false, "Strong"), "None");
    EXPECT_EQ(determine_loh_subtype_legacy_vc(true, "Noise"), "LOH_Noise");
    EXPECT_EQ(determine_loh_subtype_legacy_vc(true, "Weak"), "LOH_Weak");
    EXPECT_EQ(determine_loh_subtype_legacy_vc(true, "Strong"), "LOH_Strong");
    EXPECT_EQ(determine_loh_subtype_legacy_vc(true, "Subclone"), "LOH_Subclone");
    EXPECT_THROW(determine_loh_subtype_legacy_vc(true, "Unknown"), std::invalid_argument);
    EXPECT_THROW(determine_loh_subtype_legacy_vc(false, "Unknown"), std::invalid_argument);
}

TEST(VerificationSchemaV2Test, UnknownLegacyClassFailsClosedBeforeClassification) {
    auto input = base_input("Unknown");
    input.clean_location_permanova = true;
    EXPECT_THROW(classify_verification_v2(input), std::invalid_argument);
}

TEST(VerificationSchemaV2Test, ResolverPreservesLegacyWhenCurrentAlreadyUsesSchemaV2) {
    EXPECT_EQ(resolve_verification_legacy_input("Strong", "Noise"), "Strong");
    EXPECT_EQ(resolve_verification_legacy_input("Noise_Uncorrelated", "Noise"), "Noise");
    EXPECT_EQ(resolve_verification_legacy_input("MultiGroupNoLabel", "Weak"), "Weak");
    EXPECT_THROW(resolve_verification_legacy_input("Noise_Uncorrelated", "FutureLegacyClass"),
                 std::invalid_argument);
}
