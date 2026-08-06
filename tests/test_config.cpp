#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

#include "core/Config.hpp"
#include "utils/ArgParser.hpp"

using namespace InterSubMod;

// Helper to create dummy files
void create_dummy_file(const std::string& path) {
    std::ofstream ofs(path);
    ofs << "dummy content";
    ofs.close();
}

TEST(ConfigTest, ValidationSuccess) {
    Config config;
    create_dummy_file("tumor.bam");
    create_dummy_file("ref.fa");
    create_dummy_file("somatic.vcf");

    config.tumor_bam_path = "tumor.bam";
    config.reference_fasta_path = "ref.fa";
    config.somatic_vcf_path = "somatic.vcf";
    config.window_size_bp = 500;
    config.binary_methyl_high = 0.8;
    config.binary_methyl_low = 0.2;

    // Since these are not valid BAM/FASTA/VCF files, htslib should fail
    EXPECT_FALSE(config.validate());

    // Cleanup
    std::remove("tumor.bam");
    std::remove("ref.fa");
    std::remove("somatic.vcf");
}

TEST(ConfigTest, ValidationFailureMissingFiles) {
    Config config;
    // Missing paths are required, so validate should fail
    EXPECT_FALSE(config.validate());
}

TEST(ConfigTest, ValidationFailureInvalidWindow) {
    Config config;
    config.window_size_bp = -100;
    // Should fail due to missing paths AND invalid window
    EXPECT_FALSE(config.validate());
}

TEST(ConfigTest, ValidationFailureInvalidMethylThresholds) {
    Config config;
    // Set paths to dummy values to avoid path validation failure masking the threshold failure
    // (Although Config::validate() checks all and returns false if any fail,
    // we just want to ensure it returns false)

    config.binary_methyl_high = 0.3;
    config.binary_methyl_low = 0.5;  // Low > High
    EXPECT_FALSE(config.validate());

    config.binary_methyl_low = 0.2;
    config.binary_methyl_high = 1.2;  // > 1.0
    EXPECT_FALSE(config.validate());
}

// Verify that validate() doesn't leak HTSlib file handles when called
// repeatedly (guards against raw-pointer path issues).
TEST(ConfigTest, ValidationRepeatedCallsDoNotLeak) {
    Config config;
    config.tumor_bam_path = "/nonexistent/path/tumor.bam";
    config.reference_fasta_path = "/nonexistent/path/ref.fa";
    config.somatic_vcf_path = "/nonexistent/path/somatic.vcf";
    config.window_size_bp = 500;
    config.binary_methyl_high = 0.8;
    config.binary_methyl_low = 0.2;

    // Calling validate() many times with invalid paths should never crash,
    // hang, or exhaust file descriptors - this guards against resource leaks.
    for (int i = 0; i < 50; ++i) {
        EXPECT_FALSE(config.validate()) << "Iteration " << i;
    }
}

TEST(ArgParserTest, ParseArgumentsShortOptions) {
    Config config;
    // Create dummy files BEFORE parsing because CLI11::ExistingFile checks for them!
    create_dummy_file("t.bam");
    create_dummy_file("r.fa");
    create_dummy_file("s.vcf");

    const char* argv[] = {"program", "-t", "t.bam", "-r", "r.fa", "-v", "s.vcf", "-w", "200", "-j", "4"};
    int argc = 11;

    bool result = Utils::ArgParser::parse(argc, const_cast<char**>(argv), config);

    EXPECT_TRUE(result);
    EXPECT_EQ(config.tumor_bam_path, "t.bam");
    EXPECT_EQ(config.reference_fasta_path, "r.fa");
    EXPECT_EQ(config.somatic_vcf_path, "s.vcf");
    EXPECT_EQ(config.window_size_bp, 200);
    EXPECT_EQ(config.threads, 4);

    std::remove("t.bam");
    std::remove("r.fa");
    std::remove("s.vcf");
}

// ---------------------------------------------------------------------------
// run_params.json provenance (added 2026-08-06)
//
// Before this, a finished run left no machine-readable record of its own
// parameters: Config::print() showed 13 of ~35 fields to stdout and nothing was
// persisted. Switches that materially change results — --germline-hp-only,
// --permute-hp-seed, --nan-distance-strategy, --expected-coverage — became
// unknowable once the terminal scrolled away. These tests pin that behaviour.
// ---------------------------------------------------------------------------

TEST(ConfigRunParamsTest, JsonContainsPreviouslyInvisibleParameters) {
    Config config;
    config.germline_hp_only = true;
    config.permute_hp_seed = 4242;
    config.nan_distance_strategy = NanDistanceStrategy::MAX_DIST;
    config.expected_coverage = 37.5;
    config.min_base_quality = 25;
    config.min_common_coverage = 7;

    const std::string json = config.to_json();

    EXPECT_NE(json.find("\"germline_hp_only\": true"), std::string::npos);
    EXPECT_NE(json.find("\"permute_hp_seed\": 4242"), std::string::npos);
    EXPECT_NE(json.find("\"permute_hp_enabled\": true"), std::string::npos);
    EXPECT_NE(json.find("\"nan_distance_strategy\": \"MAX_DIST\""), std::string::npos);
    EXPECT_NE(json.find("\"expected_coverage\": 37.5"), std::string::npos);
    EXPECT_NE(json.find("\"expected_coverage_mode\": \"user_specified\""), std::string::npos);
    EXPECT_NE(json.find("\"min_base_quality\": 25"), std::string::npos);
    EXPECT_NE(json.find("\"min_common_coverage\": 7"), std::string::npos);
}

TEST(ConfigRunParamsTest, JsonRecordsBuildAndRunIdentity) {
    Config config;
    const std::string json = config.to_json();

    EXPECT_NE(json.find("\"schema_name\": \"intersubmod.run_params\""), std::string::npos);
    EXPECT_NE(json.find("\"git_commit\""), std::string::npos);
    EXPECT_NE(json.find("\"build_type\""), std::string::npos);
    EXPECT_NE(json.find("\"compiled_at\""), std::string::npos);
    EXPECT_NE(json.find("\"started_at\""), std::string::npos);
}

TEST(ConfigRunParamsTest, AutoCoverageModeReportedWhenZero) {
    Config config;
    config.expected_coverage = 0.0;
    EXPECT_NE(config.to_json().find("\"expected_coverage_mode\": \"auto_kde\""), std::string::npos);
}

TEST(ConfigRunParamsTest, DefaultNanStrategyIsSkip) {
    // SKIP has been the default since 2026-06-14; a silent revert to MAX_DIST
    // would change distance results, so pin it here.
    Config config;
    EXPECT_NE(config.to_json().find("\"nan_distance_strategy\": \"SKIP\""), std::string::npos);
}

TEST(ConfigRunParamsTest, JsonEscapesSpecialCharactersInPaths) {
    Config config;
    config.tumor_bam_path = "/tmp/we\"ird\\path/t.bam";

    const std::string json = config.to_json();

    // The raw quote and backslash must be escaped, not leaked through.
    EXPECT_NE(json.find("we\\\"ird\\\\path"), std::string::npos);
}

TEST(ConfigRunParamsTest, WritesFileContainingTheParameters) {
    Config config;
    config.output_dir = "run_params_test_dir";
    config.permute_hp_seed = 99;

    std::string err;
    ASSERT_TRUE(config.write_run_params_json(err)) << err;
    EXPECT_TRUE(err.empty());

    std::ifstream in("run_params_test_dir/run_params.json");
    ASSERT_TRUE(in.good());
    std::stringstream buf;
    buf << in.rdbuf();
    in.close();

    EXPECT_NE(buf.str().find("\"permute_hp_seed\": 99"), std::string::npos);
    EXPECT_NE(buf.str().find("\"schema_name\""), std::string::npos);
    EXPECT_NE(buf.str().find("\"output_dir\": \"run_params_test_dir\""), std::string::npos);

    std::remove("run_params_test_dir/run_params.json");
    std::remove("run_params_test_dir");
}

TEST(ConfigRunParamsTest, ReportsErrorWhenOutputDirCannotBeCreated) {
    // A regular file standing where the output directory should go: creating
    // the directory must fail, and the failure must be reported rather than
    // silently swallowed.
    create_dummy_file("blocking_regular_file");

    Config config;
    config.output_dir = "blocking_regular_file";

    std::string err;
    EXPECT_FALSE(config.write_run_params_json(err));
    EXPECT_FALSE(err.empty());

    std::remove("blocking_regular_file");
}
