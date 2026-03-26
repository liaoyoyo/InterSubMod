#pragma once

/**
 * @file ReadAggregator.hpp
 * @brief Encapsulates the read-filtering → methylation-parsing → matrix-building pipeline.
 *
 * Extracted from RegionProcessor::process_reads() to give it a clear, testable interface.
 */

#include <cstdint>
#include <string>
#include <vector>

#include "core/MatrixBuilder.hpp"
#include "core/ReadParser.hpp"
#include "core/SomaticSnv.hpp"

struct bam1_t;

namespace InterSubMod {

/** Configuration for ReadAggregator */
struct ReadAggregateConfig {
    ReadFilterConfig filter_config;
    bool no_filter_output       = false;  ///< Bypass quality filters (debug/no-filter mode)
    bool collect_filtered_reads = false;  ///< Capture rejected reads for downstream output
};

/** Output of ReadAggregator::aggregate() */
struct ReadAggregateResult {
    MatrixBuilder                 matrix_builder;
    std::vector<FilteredReadInfo> filtered_reads;
    int                           num_forward_reads = 0;
    int                           num_reverse_reads = 0;
};

/**
 * @brief Aggregates BAM records into a methylation MatrixBuilder.
 *
 * Single responsibility: given raw bam1_t pointers, run read filtering
 * (ReadParser), methylation parsing (MethylationParser), and populate a
 * MatrixBuilder. Caller retains ownership of all bam1_t pointers.
 *
 * This class is stateless across calls — each call to aggregate() creates
 * fresh ReadParser and MethylationParser instances.
 */
class ReadAggregator {
public:
    explicit ReadAggregator(const ReadAggregateConfig& config = ReadAggregateConfig());

    /**
     * @brief Process BAM records into a methylation matrix.
     *
     * @param reads         Raw bam1_t pointers (ownership retained by caller)
     * @param snv           Anchor somatic SNV for alt-support detection
     * @param ref_seq       Reference sequence covering the region
     * @param region_start  1-based start coordinate of ref_seq window
     * @return Aggregated result (matrix builder, filtered reads, strand counts)
     */
    ReadAggregateResult aggregate(const std::vector<bam1_t*>& reads,
                                   const SomaticSnv&            snv,
                                   const std::string&           ref_seq,
                                   int32_t                      region_start);

private:
    ReadAggregateConfig config_;
};

}  // namespace InterSubMod
