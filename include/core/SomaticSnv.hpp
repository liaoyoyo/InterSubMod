#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <optional>

namespace InterSubMod {

/**
 * @brief Maps chromosome names (string) to internal integer IDs.
 * Optimized for memory efficiency and fast lookups.
 */
class ChromIndex {
public:
    /**
     * @brief Gets existing ID or creates a new one for the given chromosome name.
     */
    int get_or_create_id(const std::string& chr_name);
    
    /**
     * @brief Finds the ID for a chromosome name.
     * @return ID if found, -1 if not found.
     */
    int find_id(const std::string& chr_name) const;
    
    /**
     * @brief Gets the chromosome name for a given ID.
     */
    std::string get_name(int chr_id) const;

private:
    std::map<std::string, int> name_to_id_;
    std::vector<std::string> id_to_name_;
};

/**
 * @brief Represents a somatic single nucleotide variant.
 */
struct SomaticSnv {
    int snv_id;             ///< Internal SNV ID
    int chr_id;             ///< Chromosome ID
    int32_t pos;            ///< 1-based Position
    char ref_base;          ///< Reference allele
    char alt_base;          ///< Alternate allele
    float qual;             ///< Variant Quality Score
    bool is_pass_filter;    ///< True if FILTER=PASS
    float somatic_conf;     ///< Somatic confidence score
    std::string info_flags; ///< Raw INFO field string
};

/**
 * @brief Container for all loaded somatic SNVs.
 */
class SomaticSnvTable {
public:
    /**
     * @brief Adds a generic SNV to the table.
     * @return The assigned internal SNV ID.
     */
    int add_snv(const SomaticSnv& snv);
    
    /**
     * @brief Returns total number of SNVs.
     */
    size_t size() const;
    
    /**
     * @brief Accessor to all stored SNVs.
     */
    const std::vector<SomaticSnv>& all() const;
    
    /**
     * @brief Exports the table to a TSV file for debugging/logging.
     */
    bool save_to_tsv(const std::string& path, const ChromIndex& chrom_index) const;
    
    /**
     * @brief Loads somatic variants from a VCF file.
     * Parses the VCF using htslib and populates the table.
     */
    bool load_from_vcf(const std::string& vcf_path, ChromIndex& chrom_index);

private:
    std::vector<SomaticSnv> snvs_;
};

} // namespace InterSubMod
