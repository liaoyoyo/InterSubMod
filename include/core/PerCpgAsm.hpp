#pragma once

/**
 * @file PerCpgAsm.hpp
 * @brief Per-CpG ASM detection and epiallele heterogeneity metrics
 *
 * Implements three literature-validated metric families:
 * 1. Per-CpG Fisher's exact test + BH-FDR (DAMEfinder, pycoMeth)
 * 2. NME - Normalized Methylation Entropy (CPEL, Jenkinson 2020)
 * 3. Epipolymorphism (methclone 2014, Metheor 2023)
 *
 * All metrics operate on ISM's existing MethylationMatrix (reads x CpGs)
 * and HP labels (HP1-family vs HP2-family).
 */

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "core/DataStructs.hpp"
#include "core/FisherExact.hpp"
#include "core/MethylationMatrix.hpp"

namespace InterSubMod {

/**
 * @brief Result of per-CpG ASM and epiallele metrics computation
 *
 * 10 fields chosen from Python PoC validation (2026-04-15):
 * - Fisher (4): n_sig, frac_sig, n_tested, max_neg_log_fdr
 * - NME (3): nme_hp1, nme_hp2, entropy_imbalance
 * - Epipolymorphism (3): epipoly_hp1, epipoly_hp2, epipoly_delta
 */
struct PerCpgAsmResult {
    bool valid = false;  ///< Whether computation succeeded

    // Per-CpG Fisher's exact test (HP1-family vs HP2-family)
    int fisher_n_sig = 0;          ///< Number of CpGs with FDR < 0.05
    double fisher_frac_sig = 0.0;  ///< Fraction of tested CpGs that are significant
    int fisher_n_tested = 0;       ///< Number of CpGs tested (both groups have >=1 read)
    double fisher_max_neg_log_fdr = 0.0;  ///< max(-log10(FDR)) across CpGs

    // NME - Normalized Methylation Entropy per HP group
    double nme_hp1 = std::numeric_limits<double>::quiet_NaN();  ///< NME for HP1-family
    double nme_hp2 = std::numeric_limits<double>::quiet_NaN();  ///< NME for HP2-family
    double entropy_imbalance = std::numeric_limits<double>::quiet_NaN();  ///< |NME_HP1 - NME_HP2|

    // Epipolymorphism (4-CpG sliding window)
    double epipoly_hp1 = std::numeric_limits<double>::quiet_NaN();  ///< Epipolymorphism for HP1-family
    double epipoly_hp2 = std::numeric_limits<double>::quiet_NaN();  ///< Epipolymorphism for HP2-family
    double epipoly_delta = std::numeric_limits<double>::quiet_NaN();  ///< |epipoly_hp1 - epipoly_hp2|
};

/**
 * @brief Compute per-CpG ASM metrics on a MethylationMatrix
 *
 * @param meth_mat MethylationMatrix (reads x CpGs, raw probabilities)
 * @param read_list Read information (for HP tags)
 * @param min_reads_per_group Minimum reads in each HP group (default 5)
 * @return PerCpgAsmResult with all 10 metrics
 */
PerCpgAsmResult compute_per_cpg_asm(const MethylationMatrix& meth_mat,
                                     const std::vector<ReadInfo>& read_list,
                                     int min_reads_per_group = 5);

/**
 * @brief Benjamini-Hochberg FDR correction
 *
 * @param p_values Raw p-values (may contain NaN for untested CpGs)
 * @return FDR-corrected q-values (NaN preserved)
 */
std::vector<double> bh_fdr_correction(const std::vector<double>& p_values);

/**
 * @brief Compute NME (Normalized Methylation Entropy) for a group of reads
 *
 * NME = H / H_max where H = Shannon entropy of binary methylation patterns,
 * H_max = log2(min(n_reads, 2^n_cpgs)).
 *
 * @param binary_rows Binary methylation patterns (n_reads x n_cpgs, 0/1/-1)
 * @param min_cpgs Minimum CpGs required (default 2)
 * @return NME value [0,1], or NaN if insufficient data
 */
double compute_nme(const Eigen::MatrixXi& binary_rows, int min_cpgs = 2);

/**
 * @brief Compute epipolymorphism for a group of reads
 *
 * Uses 4-CpG sliding windows. Epipoly = 1 - sum(p_i^2) where p_i are
 * pattern frequencies in each window. Returns mean across all windows.
 *
 * @param binary_rows Binary methylation patterns (n_reads x n_cpgs, 0/1/-1)
 * @param window_size Sliding window size (default 4)
 * @return Epipolymorphism [0,1], or NaN if < window_size valid CpGs
 */
double compute_epipolymorphism(const Eigen::MatrixXi& binary_rows, int window_size = 4);

}  // namespace InterSubMod
