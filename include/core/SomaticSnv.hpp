#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <optional>

namespace InterSubMod {

// Optimized chromosome string <-> ID mapping
class ChromIndex {
public:
    int get_or_create_id(const std::string& chr_name);
    int find_id(const std::string& chr_name) const;
    std::string get_name(int chr_id) const;

private:
    std::map<std::string, int> name_to_id_;
    std::vector<std::string> id_to_name_;
};

struct SomaticSnv {
    int snv_id;
    int chr_id;
    int32_t pos; // 1-based
    char ref_base;
    char alt_base;
    float qual;
    bool is_pass_filter;
    float somatic_conf;
    std::string info_flags;
};

class SomaticSnvTable {
public:
    int add_snv(const SomaticSnv& snv);
    size_t size() const;
    const std::vector<SomaticSnv>& all() const;
    bool save_to_tsv(const std::string& path, const ChromIndex& chrom_index) const;
    bool load_from_vcf(const std::string& vcf_path, ChromIndex& chrom_index);

private:
    std::vector<SomaticSnv> snvs_;
};

} // namespace InterSubMod

