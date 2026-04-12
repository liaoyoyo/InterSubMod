#include "core/SubcloneAnalyzer.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

#include "core/RegionProcessor.hpp"

namespace InterSubMod {

// Thresholds for subclone stratification
static constexpr double LOH_HP_RATIO_THRESHOLD = 0.1;  // hp_ratio < 0.1 or > 0.9
static constexpr double HIGH_HP_ASM_THRESHOLD = 0.05;   // |hp_asm_delta| > 0.05
static constexpr double HIGH_SAMPLE_ASM_THRESHOLD = 0.1; // |sample_asm_delta| > 0.1
static constexpr int MIN_READS_FOR_PROFILE = 10;
static constexpr int MIN_CPGS_FOR_PROFILE = 3;

SubcloneResult SubcloneAnalyzer::analyze(const std::vector<RegionResult>& results, int min_regions) const {
    SubcloneResult result;

    // Extract valid profiles
    auto profiles = extract_profiles(results);
    result.total_regions_analyzed = static_cast<int>(profiles.size());

    if (static_cast<int>(profiles.size()) < min_regions) {
        result.valid = false;
        return result;
    }

    // Assign subclones
    auto assignments = assign_subclones(profiles);

    // Count actual subclone groups used
    int max_id = *std::max_element(assignments.begin(), assignments.end());
    result.n_subclones = max_id + 1;
    result.region_assignments = assignments;

    // Count assigned regions (exclude -1)
    result.total_regions_assigned = 0;
    for (int a : assignments) {
        if (a >= 0) ++result.total_regions_assigned;
    }

    // Compute per-subclone summaries
    result.summaries = compute_summaries(profiles, assignments, result.n_subclones);

    // Compute overall statistics
    double sum_hp = 0.0, sum_sample = 0.0;
    int n_loh = 0;
    for (const auto& p : profiles) {
        sum_hp += std::abs(p.hp_asm_delta);
        sum_sample += std::abs(p.sample_asm_delta);
        if (p.potential_loh || p.loh_bed_overlap) ++n_loh;
    }
    result.mean_hp_asm_all = sum_hp / profiles.size();
    result.mean_sample_asm_all = sum_sample / profiles.size();
    result.loh_fraction = static_cast<double>(n_loh) / profiles.size();

    result.valid = true;
    return result;
}

std::vector<RegionAsmProfile> SubcloneAnalyzer::extract_profiles(
    const std::vector<RegionResult>& results) const {
    std::vector<RegionAsmProfile> profiles;
    profiles.reserve(results.size());

    for (const auto& r : results) {
        if (!r.success) continue;
        if (r.num_reads < MIN_READS_FOR_PROFILE) continue;
        if (r.num_cpgs < MIN_CPGS_FOR_PROFILE) continue;
        if (!r.significance_computed) continue;

        RegionAsmProfile p;
        p.region_id = r.region_id;
        p.hp_asm_delta = r.hp_merged_delta;
        p.sample_asm_delta = r.sample_asm_delta;
        p.hp_fine_f = r.hp_fine_f;
        p.allele_delta = r.allele_delta;
        p.hp_ratio = r.hp_ratio;
        p.coverage_multiple = r.coverage_multiple;
        p.potential_loh = r.potential_loh;
        p.loh_bed_overlap = r.loh_bed_overlap;
        p.loh_source = r.loh_source;
        p.verification_class = r.verification_class;
        p.quality_score = r.quality_score;
        p.num_reads = r.num_reads;
        p.num_cpgs = r.num_cpgs;

        profiles.push_back(p);
    }

    return profiles;
}

std::vector<int> SubcloneAnalyzer::assign_subclones(
    const std::vector<RegionAsmProfile>& profiles) const {
    std::vector<int> assignments(profiles.size());

    for (size_t i = 0; i < profiles.size(); ++i) {
        const auto& p = profiles[i];

        // Stratification logic:
        bool is_loh = p.potential_loh || p.loh_bed_overlap;
        bool high_hp_asm = std::abs(p.hp_asm_delta) > HIGH_HP_ASM_THRESHOLD;
        bool high_sample_asm = std::abs(p.sample_asm_delta) > HIGH_SAMPLE_ASM_THRESHOLD;

        if (is_loh) {
            // LOH regions: separate group
            assignments[i] = 2;
        } else if (high_sample_asm) {
            // High tumor vs normal difference
            assignments[i] = 3;
        } else if (high_hp_asm) {
            // Epigenetic heterogeneity (HP ASM without LOH)
            assignments[i] = 1;
        } else {
            // Normal diploid, low ASM
            assignments[i] = 0;
        }
    }

    return assignments;
}

std::vector<SubcloneSummary> SubcloneAnalyzer::compute_summaries(
    const std::vector<RegionAsmProfile>& profiles,
    const std::vector<int>& assignments,
    int n_subclones) const {
    std::vector<SubcloneSummary> summaries(n_subclones);

    // Initialize
    for (int i = 0; i < n_subclones; ++i) {
        summaries[i].subclone_id = i;
        summaries[i].n_regions = 0;
        summaries[i].mean_hp_asm = 0.0;
        summaries[i].mean_sample_asm = 0.0;
        summaries[i].mean_hp_ratio = 0.0;
        summaries[i].mean_hp_fine_f = 0.0;
        summaries[i].mean_coverage_multiple = 0.0;
        summaries[i].n_loh_bed = 0;
        summaries[i].n_loh_ratio = 0;
        summaries[i].n_loh_both = 0;
        summaries[i].n_strong = 0;
        summaries[i].n_subclone = 0;
        summaries[i].n_weak = 0;
        summaries[i].n_noise = 0;
    }

    // Accumulate
    for (size_t i = 0; i < profiles.size(); ++i) {
        int sid = assignments[i];
        if (sid < 0 || sid >= n_subclones) continue;

        auto& s = summaries[sid];
        const auto& p = profiles[i];

        s.n_regions++;
        s.mean_hp_asm += std::abs(p.hp_asm_delta);
        s.mean_sample_asm += std::abs(p.sample_asm_delta);
        s.mean_hp_ratio += p.hp_ratio;
        s.mean_hp_fine_f += p.hp_fine_f;
        s.mean_coverage_multiple += p.coverage_multiple;

        // LOH source counting
        if (p.loh_source == "bed_only") ++s.n_loh_bed;
        else if (p.loh_source == "ratio_only") ++s.n_loh_ratio;
        else if (p.loh_source == "both") { ++s.n_loh_bed; ++s.n_loh_ratio; ++s.n_loh_both; }

        // Verification class counting
        if (p.verification_class == "Strong") ++s.n_strong;
        else if (p.verification_class == "Subclone") ++s.n_subclone;
        else if (p.verification_class == "Weak") ++s.n_weak;
        else ++s.n_noise;
    }

    // Compute means and dominant class
    for (auto& s : summaries) {
        if (s.n_regions > 0) {
            s.mean_hp_asm /= s.n_regions;
            s.mean_sample_asm /= s.n_regions;
            s.mean_hp_ratio /= s.n_regions;
            s.mean_hp_fine_f /= s.n_regions;
            s.mean_coverage_multiple /= s.n_regions;
        }

        // Determine dominant verification class
        int max_count = s.n_noise;
        s.dominant_verification_class = "Noise";
        if (s.n_strong > max_count) { max_count = s.n_strong; s.dominant_verification_class = "Strong"; }
        if (s.n_subclone > max_count) { max_count = s.n_subclone; s.dominant_verification_class = "Subclone"; }
        if (s.n_weak > max_count) { s.dominant_verification_class = "Weak"; }
    }

    return summaries;
}

void SubcloneAnalyzer::write_report(const SubcloneResult& result, const std::string& output_path) {
    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        std::cerr << "[SubcloneAnalyzer] Failed to open report file: " << output_path << "\n";
        return;
    }

    ofs << "=== InterSubMod Cross-Region Subclone Analysis ===\n\n";

    if (!result.valid) {
        ofs << "Analysis: INVALID (insufficient regions)\n";
        ofs << "Regions analyzed: " << result.total_regions_analyzed << "\n";
        ofs.close();
        return;
    }

    ofs << "Total regions analyzed: " << result.total_regions_analyzed << "\n";
    ofs << "Total regions assigned: " << result.total_regions_assigned << "\n";
    ofs << "Number of subclone groups: " << result.n_subclones << "\n";
    ofs << "Overall mean |HP ASM|: " << std::fixed << std::setprecision(4) << result.mean_hp_asm_all << "\n";
    ofs << "Overall mean |Sample ASM|: " << std::fixed << std::setprecision(4) << result.mean_sample_asm_all << "\n";
    ofs << "LOH fraction: " << std::fixed << std::setprecision(3) << result.loh_fraction
        << " (" << static_cast<int>(result.loh_fraction * result.total_regions_analyzed)
        << "/" << result.total_regions_analyzed << ")\n\n";

    // Subclone group names
    static const char* group_names[] = {
        "Normal Diploid (low ASM)",
        "Epigenetic Heterogeneity (high HP ASM)",
        "LOH Regions",
        "Tumor-Specific Changes (high Sample ASM)"
    };

    ofs << "--- Per-Subclone Summary ---\n\n";
    for (const auto& s : result.summaries) {
        const char* name = (s.subclone_id < 4) ? group_names[s.subclone_id] : "Unknown";
        ofs << "Subclone " << s.subclone_id << ": " << name << "\n";
        ofs << "  Regions: " << s.n_regions << " ("
            << std::fixed << std::setprecision(1)
            << (100.0 * s.n_regions / result.total_regions_analyzed) << "%)\n";
        ofs << "  Mean |HP ASM|: " << std::setprecision(4) << s.mean_hp_asm << "\n";
        ofs << "  Mean |Sample ASM|: " << std::setprecision(4) << s.mean_sample_asm << "\n";
        ofs << "  Mean HP ratio: " << std::setprecision(3) << s.mean_hp_ratio << "\n";
        ofs << "  Mean HP fine F: " << std::setprecision(2) << s.mean_hp_fine_f << "\n";
        ofs << "  Mean coverage multiple: " << std::setprecision(2) << s.mean_coverage_multiple << "\n";
        ofs << "  LOH: bed=" << s.n_loh_bed << " ratio=" << s.n_loh_ratio << " both=" << s.n_loh_both << "\n";
        ofs << "  Verification: Strong=" << s.n_strong << " Subclone=" << s.n_subclone
            << " Weak=" << s.n_weak << " Noise=" << s.n_noise
            << " (dominant: " << s.dominant_verification_class << ")\n\n";
    }

    ofs.close();
}

void SubcloneAnalyzer::write_assignments_tsv(const std::vector<RegionAsmProfile>& profiles,
                                              const std::vector<int>& assignments,
                                              const std::string& output_path) {
    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        std::cerr << "[SubcloneAnalyzer] Failed to open assignments file: " << output_path << "\n";
        return;
    }

    ofs << "RegionID\tSubcloneID\tHP_ASM_Delta\tSample_ASM_Delta\tHP_Fine_F\t"
        << "HP_Ratio\tCoverage_Multiple\tPotential_LOH\tLOH_Bed_Overlap\tLOH_Source\t"
        << "VerificationClass\tQualityScore\tNumReads\tNumCpGs\n";

    for (size_t i = 0; i < profiles.size(); ++i) {
        const auto& p = profiles[i];
        int sid = (i < assignments.size()) ? assignments[i] : -1;

        ofs << p.region_id << "\t" << sid << "\t"
            << std::fixed << std::setprecision(4) << p.hp_asm_delta << "\t"
            << p.sample_asm_delta << "\t"
            << std::setprecision(2) << p.hp_fine_f << "\t"
            << std::setprecision(3) << p.hp_ratio << "\t"
            << std::setprecision(2) << p.coverage_multiple << "\t"
            << (p.potential_loh ? "true" : "false") << "\t"
            << (p.loh_bed_overlap ? "true" : "false") << "\t"
            << p.loh_source << "\t"
            << p.verification_class << "\t"
            << std::setprecision(1) << p.quality_score << "\t"
            << p.num_reads << "\t" << p.num_cpgs << "\n";
    }

    ofs.close();
}

}  // namespace InterSubMod
