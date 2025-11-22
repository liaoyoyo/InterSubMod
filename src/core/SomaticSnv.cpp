#include "core/SomaticSnv.hpp"
#include "utils/Logger.hpp"
#include <htslib/vcf.h>
#include <iostream>
#include <fstream>
#include <algorithm>

namespace InterSubMod {

// ==================================================
// ChromIndex Implementation
// ==================================================

int ChromIndex::get_or_create_id(const std::string& chr_name) {
    auto it = name_to_id_.find(chr_name);
    if (it != name_to_id_.end()) {
        return it->second;
    }
    
    int new_id = static_cast<int>(id_to_name_.size());
    name_to_id_[chr_name] = new_id;
    id_to_name_.push_back(chr_name);
    return new_id;
}

int ChromIndex::find_id(const std::string& chr_name) const {
    auto it = name_to_id_.find(chr_name);
    if (it != name_to_id_.end()) {
        return it->second;
    }
    return -1;
}

std::string ChromIndex::get_name(int chr_id) const {
    if (chr_id >= 0 && chr_id < static_cast<int>(id_to_name_.size())) {
        return id_to_name_[chr_id];
    }
    return "";
}

// ==================================================
// SomaticSnvTable Implementation
// ==================================================

int SomaticSnvTable::add_snv(const SomaticSnv& snv) {
    int id = static_cast<int>(snvs_.size());
    SomaticSnv snv_copy = snv;
    snv_copy.snv_id = id; // Ensure ID is consistent with index
    snvs_.push_back(snv_copy);
    return id;
}

size_t SomaticSnvTable::size() const {
    return snvs_.size();
}

const std::vector<SomaticSnv>& SomaticSnvTable::all() const {
    return snvs_;
}

bool SomaticSnvTable::save_to_tsv(const std::string& path, const ChromIndex& chrom_index) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        Utils::Logger::error("Failed to open file for writing: " + path);
        return false;
    }

    out << "snv_id\tchr\tpos\tref\talt\tqual\tfilter\tsomatic_conf\n";
    for (const auto& snv : snvs_) {
        out << snv.snv_id << "\t"
            << chrom_index.get_name(snv.chr_id) << "\t"
            << snv.pos << "\t"
            << snv.ref_base << "\t"
            << snv.alt_base << "\t"
            << snv.qual << "\t"
            << (snv.is_pass_filter ? "PASS" : "FAIL") << "\t"
            << snv.somatic_conf << "\n";
    }
    
    out.close();
    return true;
}

bool SomaticSnvTable::load_from_vcf(const std::string& vcf_path, ChromIndex& chrom_index) {
    Utils::Logger::info("Starting to load SNVs from VCF: " + vcf_path);

    // 1. Open file
    htsFile *fp = vcf_open(vcf_path.c_str(), "r");
    if (!fp) {
        Utils::Logger::error("Failed to open VCF file: " + vcf_path);
        return false;
    }

    // 2. Read Header
    bcf_hdr_t *hdr = vcf_hdr_read(fp);
    if (!hdr) {
        Utils::Logger::error("Failed to read VCF header: " + vcf_path);
        vcf_close(fp);
        return false;
    }

    // Pre-populate ChromIndex from header to ensure ID consistency
    int nseq = 0;
    const char **seqnames = bcf_hdr_seqnames(hdr, &nseq);
    if (seqnames) {
        for (int i = 0; i < nseq; i++) {
            chrom_index.get_or_create_id(seqnames[i]);
        }
        free(seqnames);
    }
    Utils::Logger::info("Loaded " + std::to_string(nseq) + " contigs from header.");

    // Verify PASS filter exists in header (warning only)
    if (bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS") < 0) {
        Utils::Logger::warning("Filter 'PASS' not found in VCF header definitions.");
    }

    // 3. Read Records
    bcf1_t *rec = bcf_init();
    int snv_count = 0;
    int skipped_count = 0;

    while (vcf_read(fp, hdr, rec) >= 0) {
        // Unpack everything including FORMAT
        // Note: For performance, if we don't need FORMAT every time, we could unpack selectively.
        // But we need AF from FORMAT.
        bcf_unpack(rec, BCF_UN_ALL);

        // Check FILTER == PASS
        // bcf_has_filter returns 1 if present, 0 if not.
        // Note: If no filters are present (.), it returns 0 for specific query? 
        // Actually bcf_has_filter logic: if filter is "PASS", the id in rec->d.flt is the id of PASS.
        // Standard VCF: PASS is usually ID 0 or defined in header.
        // bcf_has_filter handles string lookup.
        if (bcf_has_filter(hdr, rec, const_cast<char*>("PASS")) != 1) {
            skipped_count++;
            continue;
        }

        // Check if it's a biallelic SNP
        if (!bcf_is_snp(rec) || rec->n_allele != 2) {
            skipped_count++;
            continue;
        }

        // Get somatic confidence (AF from Tumor)
        // In ClairS, it's a single sample VCF but with Tumor stats.
        float tumor_vaf = 0.0f;
        float *af_ptr = NULL;
        int n_af = 0;
        int n_ret = bcf_get_format_float(hdr, rec, "AF", &af_ptr, &n_af);
        
        if (n_ret > 0 && af_ptr) {
            // usually 1 value per sample. Since it's single sample VCF (or we assume the first one)
            tumor_vaf = af_ptr[0];
        }
        if (af_ptr) free(af_ptr);

        // Create SomaticSnv object
        SomaticSnv snv;
        // Since we pre-populated ChromIndex in same order as header, 
        // rec->rid should map to our internal IDs if we added them in order 0..N.
        // However, ChromIndex::get_or_create_id ensures uniqueness. 
        // To be perfectly safe, we map by name again or rely on the fact we added all from header.
        // bcf_hdr_id2name(hdr, rec->rid) gives string.
        const char* chrom_name = bcf_hdr_id2name(hdr, rec->rid);
        snv.chr_id = chrom_index.get_or_create_id(chrom_name);
        
        snv.pos = rec->pos + 1; // 0-based to 1-based
        snv.ref_base = rec->d.allele[0][0];
        snv.alt_base = rec->d.allele[1][0];
        snv.qual = rec->qual;
        snv.is_pass_filter = true;
        snv.somatic_conf = tumor_vaf;
        // snv.info_flags could be parsed from INFO if needed

        this->add_snv(snv);
        snv_count++;
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    vcf_close(fp);

    Utils::Logger::info("Finished loading VCF. Loaded: " + std::to_string(snv_count) + ", Skipped: " + std::to_string(skipped_count));
    return true;
}

} // namespace InterSubMod

