#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct ReadLengthSummary {
    std::vector<uint64_t> counts;
    uint64_t records_seen = 0;
    uint64_t records_used = 0;
    uint64_t records_filtered = 0;
    uint64_t records_without_sequence = 0;
    uint64_t total_bases = 0;
    uint32_t min_length = std::numeric_limits<uint32_t>::max();
    uint32_t max_length = 0;
};

uint32_t weighted_quantile_length(const ReadLengthSummary& summary, double retained_base_fraction) {
    if (summary.total_bases == 0) {
        return 0;
    }

    const long double target =
        static_cast<long double>(summary.total_bases) * static_cast<long double>(retained_base_fraction);
    long double cumulative = 0;
    for (size_t length = summary.counts.size(); length-- > 0;) {
        cumulative += static_cast<long double>(length) * static_cast<long double>(summary.counts[length]);
        if (cumulative >= target) {
            return static_cast<uint32_t>(length);
        }
    }
    return 0;
}

uint32_t median_read_length(const ReadLengthSummary& summary) {
    if (summary.records_used == 0) {
        return 0;
    }

    const uint64_t target = (summary.records_used + 1) / 2;
    uint64_t cumulative = 0;
    for (size_t length = 0; length < summary.counts.size(); ++length) {
        cumulative += summary.counts[length];
        if (cumulative >= target) {
            return static_cast<uint32_t>(length);
        }
    }
    return 0;
}

void print_usage(const char* program) {
    std::cerr << "Usage: " << program
              << " <input.bam> [threads=8] [exclude_flags=0x900] [progress_records=5000000]\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 5) {
        print_usage(argv[0]);
        return 2;
    }

    const std::string bam_path = argv[1];
    const int threads = argc >= 3 ? std::max(0, std::atoi(argv[2])) : 8;
    const uint32_t exclude_flags =
        argc >= 4 ? static_cast<uint32_t>(std::strtoul(argv[3], nullptr, 0)) : 0x900;
    const uint64_t progress_records =
        argc >= 5 ? static_cast<uint64_t>(std::strtoull(argv[4], nullptr, 0)) : 5000000;

    samFile* input = sam_open(bam_path.c_str(), "r");
    if (input == nullptr) {
        std::cerr << "ERROR\tcannot_open_bam\t" << bam_path << '\n';
        return 1;
    }
    if (threads > 0 && hts_set_threads(input, threads) != 0) {
        std::cerr << "ERROR\tcannot_set_threads\t" << threads << '\n';
        sam_close(input);
        return 1;
    }

    sam_hdr_t* header = sam_hdr_read(input);
    if (header == nullptr) {
        std::cerr << "ERROR\tcannot_read_header\t" << bam_path << '\n';
        sam_close(input);
        return 1;
    }

    bam1_t* record = bam_init1();
    if (record == nullptr) {
        std::cerr << "ERROR\tcannot_allocate_record\n";
        sam_hdr_destroy(header);
        sam_close(input);
        return 1;
    }

    ReadLengthSummary summary;
    int read_status = 0;
    while ((read_status = sam_read1(input, header, record)) >= 0) {
        ++summary.records_seen;
        if ((record->core.flag & exclude_flags) != 0) {
            ++summary.records_filtered;
            continue;
        }

        const uint32_t length = record->core.l_qseq;
        if (length == 0) {
            ++summary.records_without_sequence;
            continue;
        }
        if (summary.counts.size() <= length) {
            summary.counts.resize(static_cast<size_t>(length) + 1, 0);
        }
        ++summary.counts[length];
        ++summary.records_used;
        summary.total_bases += length;
        summary.min_length = std::min(summary.min_length, length);
        summary.max_length = std::max(summary.max_length, length);

        if (progress_records > 0 && summary.records_seen % progress_records == 0) {
            std::cerr << "PROGRESS\trecords_seen\t" << summary.records_seen << '\n';
        }
    }

    bam_destroy1(record);
    sam_hdr_destroy(header);
    const int close_status = sam_close(input);
    if (read_status < -1 || close_status != 0) {
        std::cerr << "ERROR\tbam_read_failed\tstatus=" << read_status << "\tclose_status=" << close_status << '\n';
        return 1;
    }
    if (summary.records_used == 0) {
        std::cerr << "ERROR\tno_records_with_sequence_after_filtering\n";
        return 1;
    }

    const uint32_t read_n50 = weighted_quantile_length(summary, 0.50);
    const uint32_t read_n90 = weighted_quantile_length(summary, 0.90);
    const uint32_t median = median_read_length(summary);
    const long double mean =
        static_cast<long double>(summary.total_bases) / static_cast<long double>(summary.records_used);

    std::cout << "metric\tvalue\n"
              << "input_bam\t" << bam_path << '\n'
              << "exclude_flags\t0x" << std::hex << std::uppercase << exclude_flags << std::dec << '\n'
              << "records_seen\t" << summary.records_seen << '\n'
              << "records_used\t" << summary.records_used << '\n'
              << "records_filtered\t" << summary.records_filtered << '\n'
              << "records_without_sequence\t" << summary.records_without_sequence << '\n'
              << "total_bases\t" << summary.total_bases << '\n'
              << "minimum_read_length_bp\t" << summary.min_length << '\n'
              << "maximum_read_length_bp\t" << summary.max_length << '\n'
              << "mean_read_length_bp\t" << std::fixed << std::setprecision(2) << mean << '\n'
              << "median_read_length_bp\t" << median << '\n'
              << "read_N50_bp\t" << read_n50 << '\n'
              << "read_N90_bp\t" << read_n90 << '\n';
    return 0;
}
