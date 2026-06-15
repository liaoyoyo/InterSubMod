#pragma once

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "core/BamReader.hpp"
#include "core/Config.hpp"
#include "core/LohBedAnnotator.hpp"
#include "core/PerCpgAsm.hpp"
#include "core/SubcloneAnalyzer.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/HierarchicalClustering.hpp"
#include "core/MatrixBuilder.hpp"
#include "core/MethylationMatrix.hpp"
#include "core/MethylationParser.hpp"
#include "core/ReadAggregator.hpp"
#include "core/ReadParser.hpp"
#include "core/RegionBounds.hpp"
#include "core/SignificanceAnalyzer.hpp"
#include "core/SomaticSnv.hpp"
#include "core/TreeStructure.hpp"
#include "io/RegionWriter.hpp"
#include "io/TreeWriter.hpp"
#include "utils/FastaReader.hpp"

namespace InterSubMod {

/**
 * @brief Statistics result for processing a single SNV region
 */
struct RegionResult {
    int region_id;
    int snv_id;
    int num_reads;
    int num_cpgs;
    int num_reads_valid = 0;  ///< Reads used for clustering after NaN-pair filtering (SKIP); == num_reads under MAX_DIST
    double hp_auc_normal = -1.0;  ///< HP-AUC: distance recovers germline-HP on NORMAL reads, P(diff>same). -1=undefined
    double hp_auc_tumor = -1.0;   ///< HP-AUC on TUMOR reads (-1 common: tumor single-HP at somatic site)
    double hp_auc_all = -1.0;     ///< HP-AUC on ALL HP-labeled reads
    int num_forward_reads;   ///< Forward strand reads
    int num_reverse_reads;   ///< Reverse strand reads
    int num_filtered_reads;  ///< Reads filtered out (debug mode)
    int num_tumor_reads;     ///< Tumor reads in region (= num_reads when no normal BAM)
    int num_normal_reads;    ///< Normal reads in region (0 when no normal BAM)
    double elapsed_ms;
    double peak_memory_mb;
    bool success;
    std::string error_message;

    // Distance matrix statistics
    int num_valid_pairs;         ///< Number of valid distance pairs
    int num_invalid_pairs;       ///< Number of invalid pairs (insufficient overlap)
    double avg_common_coverage;  ///< Average common CpG coverage per pair
    double pairwise_mean_distance;    ///< Mean pairwise distance across valid pairs
    double pairwise_median_distance;  ///< Median pairwise distance across valid pairs

    // Significance analysis results (Cluster-First)
    bool significance_computed;  ///< Whether significance analysis was run
    double global_p_value;       ///< Global association p-value (FFH)
    double cramers_v;            ///< Effect size (Cramér's V)
    double global_hp_family_p;   ///< Layer 1: HP family p-value
    double global_hp_family_v;   ///< Layer 1: HP family Cramér's V
    double global_hp_fine_p;     ///< Layer 3: HP fine-grained p-value
    double global_hp_fine_v;     ///< Layer 3: HP fine-grained Cramér's V
    int global_hp_fine_n_groups; ///< Number of valid HP fine-grained groups
    double local_best_p_value;   ///< Best cluster's local p-value
    int local_best_cluster;      ///< Best cluster ID
    double heuristic_score;      ///< Combined heuristic score [0-1]
    bool passed_gating;          ///< Whether passed gating (global_p <= 0.1)
    double cluster_permanova_f;      ///< Cluster-based PERMANOVA pseudo-F
    double cluster_permanova_p;      ///< Cluster-based PERMANOVA p-value
    bool cluster_permanova_valid;    ///< Whether cluster-based PERMANOVA is valid
    double cluster_dispersion_p;     ///< Cluster-based dispersion p-value
    bool cluster_dispersion_warning; ///< Whether cluster-based dispersion warns

    // Label-First verification results
    bool label_test_computed;    ///< Whether label-first test was run
    double label_delta;          ///< Delta = d_between - d_within (deprecated, use Stage 1)
    double label_p_value;        ///< Label permutation test p-value (deprecated, use Stage 1)
    bool label_significant;      ///< Whether label test is significant (p <= 0.05)
    std::string dominant_label;  ///< Which label dimension is most significant ("hp", "allele", "none")
    double label_hp_permanova_f;      ///< Label-based HP PERMANOVA pseudo-F
    double label_hp_permanova_p;      ///< Label-based HP PERMANOVA p-value
    bool label_hp_permanova_valid;    ///< Whether label-based HP PERMANOVA is valid
    double label_hp_dispersion_p;     ///< Label-based HP dispersion p-value
    bool label_hp_dispersion_warning; ///< Whether label-based HP dispersion warns
    double label_allele_permanova_f;      ///< Label-based allele PERMANOVA pseudo-F
    double label_allele_permanova_p;      ///< Label-based allele PERMANOVA p-value
    bool label_allele_permanova_valid;    ///< Whether label-based allele PERMANOVA is valid
    double label_allele_dispersion_p;     ///< Label-based allele dispersion p-value
    bool label_allele_dispersion_warning; ///< Whether label-based allele dispersion warns

    // Multi-Stage HP Verification results (NEW)
    // Stage 1: HP Family Merged Test
    double hp_merged_delta;      ///< Delta for (HP1+HP1-1) vs (HP2+HP2-1)
    double hp_merged_p;          ///< P-value for merged HP test
    bool hp_merged_sig;          ///< Significant at alpha
    int hp_merged_n_hp1;         ///< HP1-family count (HP1 + HP1-1)
    int hp_merged_n_hp2;         ///< HP2-family count (HP2 + HP2-1)

    // Stage 2: HP Fine-Grained Test
    double hp_fine_f;            ///< Pseudo-F statistic
    double hp_fine_p;            ///< P-value for fine-grained test
    bool hp_fine_sig;            ///< Significant at alpha
    int hp_fine_n_groups;        ///< Number of valid groups (max 4)
    int hp_fine_n_hp1;           ///< HP1 (germline) read count in fine test
    int hp_fine_n_hp1s;          ///< HP1-1 (somatic from HP1) read count
    int hp_fine_n_hp2;           ///< HP2 (germline) read count in fine test
    int hp_fine_n_hp2s;          ///< HP2-1 (somatic from HP2) read count
    // Fine-grained pairwise mean distances (6 pairs from 4 groups)
    double hp_fine_d_hp1_hp1s;   ///< d(HP1, HP1-1): same haplotype, germline vs somatic
    double hp_fine_d_hp1_hp2;    ///< d(HP1, HP2): different haplotype, both germline
    double hp_fine_d_hp1_hp2s;   ///< d(HP1, HP2-1): cross haplotype, germline vs somatic
    double hp_fine_d_hp1s_hp2;   ///< d(HP1-1, HP2): cross haplotype, somatic vs germline
    double hp_fine_d_hp1s_hp2s;  ///< d(HP1-1, HP2-1): different haplotype, both somatic
    double hp_fine_d_hp2_hp2s;   ///< d(HP2, HP2-1): same haplotype, germline vs somatic

    // Stage 3: Allele Test
    double allele_delta;         ///< Delta for ALT vs REF
    double allele_p;             ///< P-value for allele test
    bool allele_sig;             ///< Significant at alpha

    // Stage 4: Unassigned Affinity Test
    double unassigned_affinity;  ///< Affinity score for HP3/HP0 reads
    double unassigned_affinity_p;///< P-value for affinity test
    std::string unassigned_dir;  ///< Affinity direction ("HP1", "HP2", "None")
    int unassigned_n_hp3;        ///< Number of HP3 reads
    int unassigned_n_hp0;        ///< Number of HP0/unphased reads

    // Sample ASM (Tumor vs Normal) - Label-First distance test
    double sample_asm_delta;     ///< Delta = d_between - d_within for Tumor vs Normal
    double sample_asm_p;         ///< Permutation p-value for Sample ASM
    bool sample_asm_sig;         ///< Significant at alpha
    int sample_asm_n_tumor;      ///< Tumor reads in Sample ASM test
    int sample_asm_n_normal;     ///< Normal reads in Sample ASM test

    // Normal baseline statistics (Phase B)
    double normal_baseline_mean;      ///< Overall mean normal methylation
    double normal_baseline_coverage;  ///< Mean per-CpG normal read coverage

    // Somatic HP ASM (tumor-only vs normal-only HP comparison)
    double hp_residual_delta;    ///< Somatic HP ASM = tumor_hp_delta - normal_hp_delta
    double hp_residual_p;        ///< P-value for somatic HP ASM (tumor HP test)
    bool hp_residual_sig;        ///< Significant at alpha
    // Diagnostic fields for HP subset analysis
    double tumor_hp_delta;       ///< HP delta on tumor-only subset (NaN if invalid)
    bool tumor_hp_valid;         ///< Whether tumor-only HP test had sufficient groups
    double normal_hp_delta;      ///< HP delta on normal-only subset (NaN if invalid)
    bool normal_hp_valid;        ///< Whether normal-only HP test had sufficient groups
    int tumor_hp1_count;         ///< Tumor reads in HP1
    int tumor_hp2_count;         ///< Tumor reads in HP2
    int normal_hp1_count;        ///< Normal reads in HP1
    int normal_hp2_count;        ///< Normal reads in HP2
    // Audit counters for somatic HP tags (HP:i:11/21/33 in raw BAM).
    // Always populated from hp_tag_raw regardless of the germline_hp_only flag,
    // so users can diagnose self-phasing extent independently of the downstream remapping.
    int n_hp_somatic_11;         ///< Reads with raw HP tag == "1-1" (integer 11)
    int n_hp_somatic_21;         ///< Reads with raw HP tag == "2-1" (integer 21)
    int n_hp_somatic_33;         ///< Reads with raw HP tag == "3"   (integer 33)

    // Signed HP methylation delta (mean per-CpG methylation difference)
    // Unlike distance-based delta, these capture DIRECTION: HP1 meth > HP2 → positive
    double tumor_hp_signed_delta;    ///< mean(tumor_HP1_meth) - mean(tumor_HP2_meth), NaN if invalid
    double normal_hp_signed_delta;   ///< mean(normal_HP1_meth) - mean(normal_HP2_meth), NaN if invalid
    double hp_signed_residual;       ///< tumor_signed - normal_signed (somatic directional ASM change)
    double combined_hp_signed_delta; ///< mean(all_HP1_meth) - mean(all_HP2_meth) (full matrix)
    // [Δβ module] somatic residual Δβ + permutation test (read-level; #3 修法, 取代有缺陷的 hp_residual_sig)
    double somatic_residual_dbeta;    ///< (tumor: mean β HP1 - mean β HP2) - (normal: ...), read-level. NaN=invalid
    double somatic_residual_dbeta_p;  ///< permutation p: shuffle T/N sample label (HP fixed), two-sided |perm|>=|obs|
    bool somatic_residual_dbeta_sig;  ///< p <= 0.05 = tumor HP-ASM 顯著異於 normal germline baseline (somatic ASM)
    // [Δβ module stage 2] germline ASM Δβ (normal HP1f−HP2f) + fine same-hap subclone Δβ (tumor germline vs carrier)
    double germline_asm_dbeta;     ///< normal: mean β(HP1f) − mean β(HP2f), read-level. goal1 germline ASM. NaN=invalid
    double germline_asm_dbeta_p;   ///< perm p: shuffle HP label among normal reads, two-sided |perm|>=|obs|
    bool germline_asm_dbeta_sig;   ///< p <= 0.05 = germline (parental) ASM 顯著 on normal baseline
    double subclone_dbeta_hp1;     ///< tumor: mean β(HP1 germline) − mean β(HP1-1 carrier). goal3 subclone. NaN=invalid
    double subclone_dbeta_hp1_p;   ///< perm p: shuffle germline/carrier label within tumor HP1, two-sided
    bool subclone_dbeta_hp1_sig;   ///< p<=0.05 AND min(germ,carrier)>=min_group (tiny-group guard). same-hap subclone
    int subclone_dbeta_hp1_n_germ;    ///< tumor HP1 germline read count (transparency: small => fragile subclone call)
    int subclone_dbeta_hp1_n_carrier; ///< tumor HP1-1 somatic-carrier read count
    double subclone_dbeta_hp2;     ///< tumor: mean β(HP2 germline) − mean β(HP2-1 carrier). goal3 subclone. NaN=invalid
    double subclone_dbeta_hp2_p;   ///< perm p: shuffle germline/carrier label within tumor HP2, two-sided
    bool subclone_dbeta_hp2_sig;   ///< p<=0.05 AND min(germ,carrier)>=min_group (tiny-group guard). same-hap subclone
    int subclone_dbeta_hp2_n_germ;    ///< tumor HP2 germline read count
    int subclone_dbeta_hp2_n_carrier; ///< tumor HP2-1 somatic-carrier read count

    // LOH BED annotation (Phase C)
    bool loh_bed_overlap;        ///< Whether SNV position overlaps a LOH BED region
    std::string loh_source;      ///< LOH source classification: "none", "bed_only", "ratio_only", "both"
    std::string loh_bed_annotation; ///< Annotation from BED file (4th column)

    // Per-CpG ASM and epiallele metrics (Phase E)
    bool per_cpg_asm_valid;          ///< Whether per-CpG ASM computation succeeded
    int per_cpg_fisher_n_sig;        ///< Number of CpGs with FDR < 0.05
    double per_cpg_fisher_frac_sig;  ///< Fraction of tested CpGs significant
    int per_cpg_fisher_n_tested;     ///< Number of CpGs tested
    double per_cpg_fisher_max_neg_log_fdr;  ///< max(-log10(FDR))
    double per_cpg_nme_hp1;          ///< NME for HP1-family
    double per_cpg_nme_hp2;          ///< NME for HP2-family
    double per_cpg_entropy_imbalance; ///< |NME_HP1 - NME_HP2|
    double per_cpg_epipoly_hp1;      ///< Epipolymorphism for HP1-family
    double per_cpg_epipoly_hp2;      ///< Epipolymorphism for HP2-family
    double per_cpg_epipoly_delta;    ///< |epipoly_hp1 - epipoly_hp2|

    // Cross-region subclone assignment (Phase D)
    int subclone_id;             ///< Subclone group assignment (-1 = unassigned)

    // Cluster stability results (NEW)
    double cluster_stability;    ///< Stability score from subsampling [0-1]
    bool has_outlier_cluster;    ///< True if any cluster has < 3 reads after retry
    int n_clusters;              ///< Number of clusters found

    // Bidirectional verification classification (NEW)
    std::string verification_class;  ///< "Strong", "Subclone", "Weak", or "Noise"

    // Multi-Layer Validation Quality Metrics (NEW - Phase 5)
    double hp_ratio;                 ///< HP1/(HP1+HP2), range [0,1]
    bool potential_loh;              ///< True if HP ratio < 0.1 or > 0.9
    double coverage_multiple;        ///< NumReads / diploid_coverage (auto-estimated per sample)
    double diploid_coverage_used;    ///< Actual diploid coverage baseline used for CovM (audit column)
    std::string coverage_category;   ///< "Normal", "Low", "High", "CNV_Loss", "CNV_Gain", "High_Copy"
    std::string loh_subtype;         ///< "None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"
    float quality_score;             ///< Composite quality score [0-100]
    std::string quality_tier;        ///< "High" (>=70), "Medium" (40-69), "Low" (<40)

    RegionResult()
        : region_id(-1),
          snv_id(-1),
          num_reads(0),
          num_cpgs(0),
          num_forward_reads(0),
          num_reverse_reads(0),
          num_filtered_reads(0),
          num_tumor_reads(0),
          num_normal_reads(0),
          elapsed_ms(0.0),
          peak_memory_mb(0.0),
          success(false),
          num_valid_pairs(0),
          num_invalid_pairs(0),
          avg_common_coverage(0.0),
          pairwise_mean_distance(0.0),
          pairwise_median_distance(0.0),
          significance_computed(false),
          global_p_value(1.0),
          cramers_v(0.0),
          global_hp_family_p(1.0),
          global_hp_family_v(0.0),
          global_hp_fine_p(1.0),
          global_hp_fine_v(0.0),
          global_hp_fine_n_groups(0),
          local_best_p_value(1.0),
          local_best_cluster(-1),
          heuristic_score(0.0),
          passed_gating(false),
          cluster_permanova_f(0.0),
          cluster_permanova_p(1.0),
          cluster_permanova_valid(false),
          cluster_dispersion_p(1.0),
          cluster_dispersion_warning(false),
          label_test_computed(false),
          label_delta(0.0),
          label_p_value(1.0),
          label_significant(false),
          dominant_label("none"),
          label_hp_permanova_f(0.0),
          label_hp_permanova_p(1.0),
          label_hp_permanova_valid(false),
          label_hp_dispersion_p(1.0),
          label_hp_dispersion_warning(false),
          label_allele_permanova_f(0.0),
          label_allele_permanova_p(1.0),
          label_allele_permanova_valid(false),
          label_allele_dispersion_p(1.0),
          label_allele_dispersion_warning(false),
          hp_merged_delta(0.0),
          hp_merged_p(1.0),
          hp_merged_sig(false),
          hp_merged_n_hp1(0),
          hp_merged_n_hp2(0),
          hp_fine_f(0.0),
          hp_fine_p(1.0),
          hp_fine_sig(false),
          hp_fine_n_groups(0),
          hp_fine_n_hp1(0),
          hp_fine_n_hp1s(0),
          hp_fine_n_hp2(0),
          hp_fine_n_hp2s(0),
          hp_fine_d_hp1_hp1s(NAN),
          hp_fine_d_hp1_hp2(NAN),
          hp_fine_d_hp1_hp2s(NAN),
          hp_fine_d_hp1s_hp2(NAN),
          hp_fine_d_hp1s_hp2s(NAN),
          hp_fine_d_hp2_hp2s(NAN),
          allele_delta(0.0),
          allele_p(1.0),
          allele_sig(false),
          unassigned_affinity(0.0),
          unassigned_affinity_p(1.0),
          unassigned_dir("None"),
          unassigned_n_hp3(0),
          unassigned_n_hp0(0),
          sample_asm_delta(0.0),
          sample_asm_p(1.0),
          sample_asm_sig(false),
          sample_asm_n_tumor(0),
          sample_asm_n_normal(0),
          normal_baseline_mean(0.0),
          normal_baseline_coverage(0.0),
          hp_residual_delta(0.0),
          hp_residual_p(1.0),
          hp_residual_sig(false),
          tumor_hp_delta(std::numeric_limits<double>::quiet_NaN()),
          tumor_hp_valid(false),
          normal_hp_delta(std::numeric_limits<double>::quiet_NaN()),
          normal_hp_valid(false),
          tumor_hp1_count(0),
          tumor_hp2_count(0),
          normal_hp1_count(0),
          normal_hp2_count(0),
          n_hp_somatic_11(0),
          n_hp_somatic_21(0),
          n_hp_somatic_33(0),
          tumor_hp_signed_delta(std::numeric_limits<double>::quiet_NaN()),
          normal_hp_signed_delta(std::numeric_limits<double>::quiet_NaN()),
          hp_signed_residual(std::numeric_limits<double>::quiet_NaN()),
          combined_hp_signed_delta(std::numeric_limits<double>::quiet_NaN()),
          somatic_residual_dbeta(std::numeric_limits<double>::quiet_NaN()),
          somatic_residual_dbeta_p(1.0),
          somatic_residual_dbeta_sig(false),
          germline_asm_dbeta(std::numeric_limits<double>::quiet_NaN()),
          germline_asm_dbeta_p(1.0),
          germline_asm_dbeta_sig(false),
          subclone_dbeta_hp1(std::numeric_limits<double>::quiet_NaN()),
          subclone_dbeta_hp1_p(1.0),
          subclone_dbeta_hp1_sig(false),
          subclone_dbeta_hp1_n_germ(0),
          subclone_dbeta_hp1_n_carrier(0),
          subclone_dbeta_hp2(std::numeric_limits<double>::quiet_NaN()),
          subclone_dbeta_hp2_p(1.0),
          subclone_dbeta_hp2_sig(false),
          subclone_dbeta_hp2_n_germ(0),
          subclone_dbeta_hp2_n_carrier(0),
          loh_bed_overlap(false),
          loh_source("none"),
          per_cpg_asm_valid(false),
          per_cpg_fisher_n_sig(0),
          per_cpg_fisher_frac_sig(0.0),
          per_cpg_fisher_n_tested(0),
          per_cpg_fisher_max_neg_log_fdr(0.0),
          per_cpg_nme_hp1(std::numeric_limits<double>::quiet_NaN()),
          per_cpg_nme_hp2(std::numeric_limits<double>::quiet_NaN()),
          per_cpg_entropy_imbalance(std::numeric_limits<double>::quiet_NaN()),
          per_cpg_epipoly_hp1(std::numeric_limits<double>::quiet_NaN()),
          per_cpg_epipoly_hp2(std::numeric_limits<double>::quiet_NaN()),
          per_cpg_epipoly_delta(std::numeric_limits<double>::quiet_NaN()),
          subclone_id(-1),
          cluster_stability(0.0),
          has_outlier_cluster(false),
          n_clusters(0),
          verification_class("Noise"),
          hp_ratio(0.5),
          potential_loh(false),
          coverage_multiple(1.0),
          diploid_coverage_used(75.0),
          coverage_category("Normal"),
          loh_subtype("None"),
          quality_score(50.0f),
          quality_tier("Medium") {
    }
};

/**
 * @brief Core class for parallel processing of multiple SNV regions
 *
 * This class is responsible for:
 * 1. Loading SNV table
 * 2. Defining region for each SNV (e.g., +/-2000bp)
 * 3. Parallel processing of multiple regions using OpenMP
 * 4. Managing thread-local resources (BamReader, FastaReader)
 * 5. Collecting and reporting processing results for each region
 * 6. Recording filtered reads in debug mode
 *
 * Thread-safety:
 * - Each thread maintains its own BAM/FASTA readers
 * - MatrixBuilder and RegionWriter are used within critical sections
 * - Result collection is mutex-protected
 */
class RegionProcessor {
public:
    /**
     * @brief Construct RegionProcessor (simplified version, for backward compatibility)
     */
    RegionProcessor(const std::string& tumor_bam_path, const std::string& normal_bam_path,
                    const std::string& ref_fasta_path, const std::string& output_dir, int num_threads = 4,
                    int32_t window_size = 2000);

    /**
     * @brief Construct RegionProcessor (full version, using Config)
     *
     * @param config Configuration object containing all parameters
     */
    explicit RegionProcessor(const Config& config);

    /**
     * @brief Load SNV table (TSV format)
     *
     * Format: chr  pos  ref  alt  qual
     * Example: chr17  7578000  C  T  100.0
     *
     * @param snv_table_path SNV table file path
     * @return Number of SNVs successfully loaded
     */
    int load_snvs(const std::string& snv_table_path);

    /**
     * @brief Load SNVs from VCF file
     *
     * @param vcf_path Path to VCF file
     * @return Number of SNVs loaded
     */
    int load_snvs_from_vcf(const std::string& vcf_path);

    /**
     * @brief Process all SNVs (parallelized)
     *
     * @param max_snvs Maximum number of SNVs to process (0 = all)
     * @return Vector of processing results
     */
    std::vector<RegionResult> process_all_regions(int max_snvs = 0);

    /**
     * @brief Process a single region (called by OpenMP worker thread)
     *
     * @param snv SNV information
     * @param region_id Region ID
     * @param bam_reader Thread-local BAM reader
     * @param fasta_reader Thread-local FASTA reader
     * @return RegionResult
     */
    RegionResult process_single_region(const SomaticSnv& snv, int region_id, BamReader& bam_reader,
                                       FastaReader& fasta_reader, BamReader* normal_reader = nullptr);

    /**
     * @brief Get loaded SNVs list
     */
    const std::vector<SomaticSnv>& get_snvs() const {
        return snvs_;
    }

    /**
     * @brief Output processing summary report
     */
    void print_summary(const std::vector<RegionResult>& results) const;

private:
    /**
     * @brief Write significance summary CSV and statistics report
     */
    void write_significance_summary(const std::vector<RegionResult>& results) const;

    // ========== Refactored Helper Methods ==========

    /**
     * @brief Build MethylationMatrix from MatrixBuilder for distance calculation
     *
     * @param matrix_builder Source matrix builder
     * @param region_id Region identifier
     * @return Constructed MethylationMatrix
     */
    MethylationMatrix build_methylation_matrix(const MatrixBuilder& matrix_builder, int region_id);

    /**
     * @brief Compute distance matrices and optionally write to disk
     *
     * @param meth_mat MethylationMatrix input
     * @param read_list Read information list
     * @param region_dir Region output directory
     * @param metric Distance metric to use
     * @param result RegionResult to update with statistics
     * @return Computed distance matrix
     */
    DistanceMatrix compute_and_write_distance_matrix(const MethylationMatrix& meth_mat,
                                                      const std::vector<ReadInfo>& read_list,
                                                      const std::string& region_dir, DistanceMetricType metric,
                                                      RegionResult& result);

    /**
     * @brief Perform hierarchical clustering and significance analysis
     *
     * @param all_dist Distance matrix for all reads
     * @param read_list Read information list
     * @param meth_mat MethylationMatrix for analysis
     * @param clustering_dir Directory for clustering output
     * @param chr_name Chromosome name
     * @param snv SNV information
     * @param region_id Region identifier
     * @param result RegionResult to update with analysis results
     */
    void perform_clustering_and_significance(const DistanceMatrix& all_dist, const std::vector<ReadInfo>& read_list,
                                              const MethylationMatrix& meth_mat, const std::string& clustering_dir,
                                              const std::string& chr_name, const SomaticSnv& snv, int region_id,
                                              RegionResult& result);

    /**
     * @brief Retain reads forming a complete (NaN-free) distance sub-matrix.
     *
     * Greedily drops the read with the fewest non-NaN partners until no NaN pair
     * remains. Mirrors StructureTest::filter_reads_for_complete_matrix so clustering
     * and the downstream PERMANOVA use the same valid read subset. For NaN-free input
     * (MAX_DIST) all reads are retained (no-op). Empty / size < 2 => region unclusterable.
     */
    static void extract_complete_submatrix_indices(const Eigen::MatrixXd& dist, std::vector<int>& out_indices);

    /**
     * @brief HP-AUC = P(dist(different-HP read pair) > dist(same-HP pair)) over a read subset.
     *
     * Ground truth = HP family label (germline haplotype from longphase). 1.0 = perfect HP
     * separation, 0.5 = distance unrelated to HP, -1.0 = undefined (no same+diff pairs).
     * NaN distances (SKIP invalid pairs) are skipped. Rank-based, O(P log P) over pair lists.
     * @param hp_fam per-read HP family (0=HP1, 1=HP2, <0=excluded), indexed like the matrix.
     * @param idx subset of read indices to evaluate (e.g. normal-only / tumor-only / all).
     */
    static double compute_hp_auc(const Eigen::MatrixXd& dist, const std::vector<int>& hp_fam,
                                 const std::vector<int>& idx);

    /**
     * @brief somatic residual Δβ + permutation test (Δβ module, #3 hp_residual 修法).
     *
     * residual = (tumor: mean β(HP1) − mean β(HP2)) − (normal: same), read-level (per-read mean β).
     * Tests somatic HP-ASM: does tumor's HP methylation difference exceed normal's germline baseline?
     * Null: tumor & normal share HP-ASM → shuffle sample(T/N) label (HP fixed). Two-sided |perm|>=|obs|.
     * Outputs dbeta / p / sig (NaN/1.0/false if any of the 4 (sample×HP) groups is empty).
     * @param min_group sig requires every one of the 4 (sample×HP) groups to have >= min_group reads
     *        (guards against 1-read groups driving a spurious "significant" residual).
     */
    static void compute_somatic_residual_dbeta_test(const Eigen::MatrixXd& raw, const std::vector<int>& hp_fam,
                                                    const std::vector<bool>& is_tumor, int n_perm, uint64_t seed,
                                                    int min_group, double& out_dbeta, double& out_p, bool& out_sig);

    /**
     * @brief Generic two-group Δβ + label-shuffle permutation test (Δβ module stage 2 primitive).
     *
     * Δβ = mean β(group 0) − mean β(group 1), read-level (per-read mean β over valid CpG).
     * Null: group label exchangeable → shuffle the 0/1 label among labeled reads. Two-sided |perm|>=|obs|.
     * @param group per-read group id: 0 / 1 = the two compared groups, anything else = excluded.
     * Used for: germline ASM (normal HP1f vs HP2f) and fine same-hap subclone (tumor germline vs somatic-carrier).
     * @param min_group sig requires min(group0_n, group1_n) >= min_group (guards 1-read groups; subclone tiny-group
     *        artifact: ~25% of raw subclone "sig" had a group <3 reads). out_n0/out_n1 report the two group sizes.
     * Outputs dbeta / p / sig (NaN/1.0/false if either group is empty) + the two group read counts.
     */
    static void compute_group_dbeta_test(const Eigen::MatrixXd& raw, const std::vector<int>& group, int n_perm,
                                         uint64_t seed, int min_group, double& out_dbeta, double& out_p, bool& out_sig,
                                         int& out_n0, int& out_n1);

    /**
     * @brief Write strand-specific clustering trees
     *
     * @param forward_dist Forward strand distance matrix
     * @param reverse_dist Reverse strand distance matrix
     * @param read_list Read information list
     * @param clustering_dir Directory for clustering output
     */
    void write_strand_specific_trees(const DistanceMatrix& forward_dist, const DistanceMatrix& reverse_dist,
                                      const std::vector<ReadInfo>& read_list, const std::string& clustering_dir);

    std::string tumor_bam_path_;
    std::string normal_bam_path_;
    std::string ref_fasta_path_;
    std::string output_dir_;
    std::string debug_output_dir_;
    std::string vcf_filename_;  ///< VCF filename (without path and extension)
    int num_threads_;
    int32_t window_size_;

    // Configuration
    LogLevel log_level_;
    bool output_filtered_reads_;
    bool no_filter_output_;
    ReadFilterConfig filter_config_;
    bool use_full_read_span_;  ///< If true, dynamically expand window to cover full span of reads
    double expected_coverage_;  ///< User-specified diploid coverage; 0 = auto-estimate

    // Distance matrix configuration
    bool compute_distance_matrix_;
    bool output_distance_matrix_;
    bool output_strand_distance_matrices_;
    DistanceConfig distance_config_;
    std::vector<DistanceMetricType> distance_metrics_;

    // Hierarchical clustering configuration
    bool compute_clustering_;
    bool output_tree_files_;
    bool output_linkage_matrix_;
    LinkageMethod linkage_method_;
    int clustering_min_reads_;
    double binary_methyl_high_;  ///< Threshold for methylated (1) call in binary matrix
    double binary_methyl_low_;   ///< Threshold for unmethylated (0) call in binary matrix

    std::vector<SomaticSnv> snvs_;
    ChromIndex chrom_index_;  // Manage chromosome name to ID mapping

    // LOH BED annotation (Phase C)
    LohBedAnnotator loh_annotator_;  ///< LOH BED file annotator

    // Thread-local resources are created in process_single_region
};

// Free functions exposed for testing
double compute_coverage_multiple(int num_reads, double expected_coverage = 75.0);

struct DiploidEstimate {
    double value;          ///< Estimated diploid coverage (or 75.0 fallback)
    bool used_fallback;    ///< True if any of the three fallback conditions fired
};
DiploidEstimate estimate_diploid_coverage(const std::vector<RegionResult>& results);

}  // namespace InterSubMod
