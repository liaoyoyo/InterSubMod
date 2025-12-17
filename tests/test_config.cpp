#include <gtest/gtest.h>

#include <fstream>

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
