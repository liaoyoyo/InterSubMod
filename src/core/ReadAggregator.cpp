#include "core/ReadAggregator.hpp"

#include <htslib/sam.h>
#include <unordered_set>

#include "core/MethylationParser.hpp"
#include "core/Types.hpp"

namespace InterSubMod {

ReadAggregator::ReadAggregator(const ReadAggregateConfig& config) : config_(config) {}

ReadAggregateResult ReadAggregator::aggregate(const std::vector<bam1_t*>& reads,
                                               const SomaticSnv&            snv,
                                               const std::string&           ref_seq,
                                               int32_t                      region_start) {
    ReadAggregateResult out;
    ReadParser read_parser(config_.filter_config);
    MethylationParser methyl_parser;
    std::unordered_set<std::string> seen_read_names;
    int read_count = 0;

    for (auto* b : reads) {
        auto [keep, filter_reason] = read_parser.should_keep_with_reason(b);

        if (!keep && !config_.no_filter_output) {
            if (config_.collect_filtered_reads) {
                out.filtered_reads.push_back(
                    read_parser.create_filtered_info(b, true, filter_reason));
            }
            continue;
        }

        // Even in no-filter mode, skip structurally invalid records that would
        // cause undefined behaviour in methylation parsing.
        if (config_.no_filter_output) {
            uint16_t flag = b->core.flag;
            if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FDUP)) {
                continue;
            }
        }

        ReadInfo info = read_parser.parse(b, read_count, true, snv, ref_seq, region_start);

        if (info.alt_support == AltSupport::UNKNOWN && !config_.no_filter_output) {
            if (config_.collect_filtered_reads) {
                out.filtered_reads.push_back(
                    read_parser.create_filtered_info(b, true, FilterReason::SNV_NOT_COVERED));
            }
            continue;
        }

        if (seen_read_names.count(info.read_name)) {
            continue;
        }
        seen_read_names.insert(info.read_name);

        auto methyl_calls = methyl_parser.parse_read(b, ref_seq, region_start);
        out.matrix_builder.add_read(info, methyl_calls);
        read_count++;

        if (info.strand == Strand::FORWARD) {
            out.num_forward_reads++;
        } else if (info.strand == Strand::REVERSE) {
            out.num_reverse_reads++;
        }
    }

    return out;
}

}  // namespace InterSubMod
