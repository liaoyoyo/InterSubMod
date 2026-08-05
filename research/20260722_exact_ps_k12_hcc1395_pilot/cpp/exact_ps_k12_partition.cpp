// Independent C++17 parity kernel for exact-PS bounded contiguous partitioning.

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace fs = std::filesystem;
using BigInt = boost::multiprecision::cpp_int;

namespace {

constexpr const char* kRetained = "retained";
constexpr const char* kCut = "cut";
constexpr const char* kUnavoidable = "unavoidable_span_gt_max_block_size";

const std::vector<std::string> kUnitHeader = {
    "dataset",       "chrom",         "unit_id",     "component_id", "linkage_basis",
    "hp_family",     "phase_set",     "threshold",   "k",            "positions1",
};

const std::vector<std::string> kConstraintHeader = {
    "dataset",          "chrom",           "unit_id",        "constraint_id",
    "hp_family",        "phase_set",       "positions1",     "call_codes",
    "molecule_weight",  "pattern_count",
};

struct CliOptions {
    fs::path units_path;
    fs::path constraints_path;
    fs::path output_dir;
    std::size_t max_block_size = 12;
};

struct Constraint {
    std::string dataset;
    std::string chrom;
    std::string unit_id;
    std::string constraint_id;
    std::string hp_family;
    std::string phase_set;
    std::vector<std::int64_t> positions;
    std::string call_codes;
    BigInt molecule_weight = 0;
    BigInt pattern_count = 0;
    std::size_t lo = 0;
    std::size_t hi = 0;
};

struct Unit {
    std::string dataset;
    std::string chrom;
    std::string unit_id;
    std::string component_id;
    std::string linkage_basis;
    std::string hp_family;
    std::string phase_set;
    std::uint64_t threshold = 0;
    std::vector<std::int64_t> positions;
    std::vector<Constraint> constraints;
};

struct Metric {
    BigInt weight = 0;
    BigInt pattern_count = 0;
};

struct State {
    BigInt retained_weight = 0;
    BigInt retained_pattern_count = 0;
    std::size_t block_count = 0;
    BigInt cut_gap_sum = 0;
    std::vector<std::size_t> cuts;
};

struct Disposition {
    const Constraint* constraint = nullptr;
    std::string status;
    std::size_t span_sites = 0;
    std::size_t crossed_cut_count = 0;
    std::optional<std::size_t> block_index;
};

struct Block {
    std::size_t block_index = 0;
    std::size_t start_index = 0;
    std::size_t end_index_exclusive = 0;
    std::vector<std::string> retained_constraint_ids;
    BigInt retained_weight = 0;
    BigInt retained_pattern_count = 0;
};

struct PartitionResult {
    const Unit* unit = nullptr;
    std::vector<std::size_t> cut_indices;
    std::vector<BigInt> cut_gaps;
    BigInt cut_gap_sum = 0;
    std::vector<Block> blocks;
    std::vector<Disposition> dispositions;
    BigInt total_weight = 0;
    BigInt retained_weight = 0;
    BigInt lost_weight = 0;
    BigInt total_pattern_count = 0;
    BigInt retained_pattern_count = 0;
    BigInt lost_pattern_count = 0;
};

struct AggregateSummary {
    std::size_t units = 0;
    std::size_t sites = 0;
    std::size_t constraints = 0;
    std::size_t blocks = 0;
    std::size_t cuts = 0;
    std::size_t retained_constraints = 0;
    std::size_t cut_constraints = 0;
    std::size_t unavoidable_constraints = 0;
    BigInt total_weight = 0;
    BigInt retained_weight = 0;
    BigInt lost_weight = 0;
    BigInt total_pattern_count = 0;
    BigInt retained_pattern_count = 0;
    BigInt lost_pattern_count = 0;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> split_preserving_empty(const std::string& value, char delimiter) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t found = value.find(delimiter, start);
        if (found == std::string::npos) {
            fields.push_back(value.substr(start));
            break;
        }
        fields.push_back(value.substr(start, found - start));
        start = found + 1;
    }
    return fields;
}

std::string join_strings(const std::vector<std::string>& values, const std::string& delimiter) {
    std::ostringstream stream;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) {
            stream << delimiter;
        }
        stream << values[index];
    }
    return stream.str();
}

template <typename Integer>
std::string join_integers(const std::vector<Integer>& values, const std::string& delimiter) {
    std::ostringstream stream;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) {
            stream << delimiter;
        }
        stream << values[index];
    }
    return stream.str();
}

std::string bigint_string(const BigInt& value) {
    return value.convert_to<std::string>();
}

std::string join_bigints(const std::vector<BigInt>& values, const std::string& delimiter) {
    std::vector<std::string> encoded;
    encoded.reserve(values.size());
    for (const BigInt& value : values) {
        encoded.push_back(bigint_string(value));
    }
    return join_strings(encoded, delimiter);
}

void require_plain_tsv(const fs::path& path, const std::string& label) {
    if (ends_with(path.string(), ".gz")) {
        fail(label + " must be plain TSV; run `gzip -cd INPUT.tsv.gz > INPUT.tsv` first");
    }
    if (!fs::is_regular_file(path)) {
        fail(label + " is not a regular file: " + path.string());
    }
}

void validate_text_field(const std::string& value, const std::string& label) {
    if (value.empty()) {
        fail(label + " must not be empty");
    }
    if (value.find_first_of("\t\r\n") != std::string::npos) {
        fail(label + " contains a forbidden control character");
    }
}

std::uint64_t parse_u64(const std::string& text, const std::string& label, bool positive) {
    if (text.empty() ||
        std::any_of(text.begin(), text.end(), [](unsigned char value) { return value < '0' || value > '9'; })) {
        fail(label + " must be an unsigned decimal integer");
    }
    std::uint64_t value = 0;
    const char* begin = text.data();
    const char* end = begin + text.size();
    const auto parsed = std::from_chars(begin, end, value);
    if (parsed.ec != std::errc() || parsed.ptr != end) {
        fail(label + " is outside uint64 range");
    }
    if (positive && value == 0) {
        fail(label + " must be positive");
    }
    return value;
}

std::int64_t parse_position(const std::string& text, const std::string& label) {
    if (text.empty()) {
        fail(label + " must not be empty");
    }
    std::int64_t value = 0;
    const char* begin = text.data();
    const char* end = begin + text.size();
    const auto parsed = std::from_chars(begin, end, value);
    if (parsed.ec != std::errc() || parsed.ptr != end) {
        fail(label + " must be an int64 coordinate");
    }
    if (value < 1) {
        fail(label + " must be a positive one-based coordinate");
    }
    return value;
}

BigInt parse_nonnegative_bigint(const std::string& text, const std::string& label) {
    if (text.empty() ||
        std::any_of(text.begin(), text.end(), [](unsigned char value) { return value < '0' || value > '9'; })) {
        fail(label + " must be a non-negative decimal integer");
    }
    BigInt value = 0;
    for (const char digit : text) {
        value *= 10;
        value += static_cast<unsigned int>(digit - '0');
    }
    return value;
}

std::vector<std::int64_t> parse_positions(const std::string& text, const std::string& label) {
    const std::vector<std::string> encoded = split_preserving_empty(text, ',');
    if (encoded.empty() || (encoded.size() == 1 && encoded.front().empty())) {
        fail(label + " must not be empty");
    }
    std::vector<std::int64_t> positions;
    positions.reserve(encoded.size());
    for (std::size_t index = 0; index < encoded.size(); ++index) {
        positions.push_back(parse_position(encoded[index], label + " item " + std::to_string(index)));
    }
    for (std::size_t index = 1; index < positions.size(); ++index) {
        if (positions[index] <= positions[index - 1]) {
            fail(label + " must be strictly increasing and unique");
        }
    }
    return positions;
}

std::vector<std::vector<std::string>> read_tsv(
    const fs::path& path,
    const std::vector<std::string>& expected_header,
    const std::string& label
) {
    std::ifstream input(path);
    if (!input) {
        fail("cannot open " + label + ": " + path.string());
    }
    std::string line;
    if (!std::getline(input, line)) {
        fail(label + " is empty");
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    const std::vector<std::string> observed_header = split_preserving_empty(line, '\t');
    if (observed_header != expected_header) {
        fail(label + " header mismatch; expected: " + join_strings(expected_header, "\\t"));
    }
    std::vector<std::vector<std::string>> rows;
    std::size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            fail(label + " contains a blank row at line " + std::to_string(line_number));
        }
        std::vector<std::string> fields = split_preserving_empty(line, '\t');
        if (fields.size() != expected_header.size()) {
            fail(label + " line " + std::to_string(line_number) + " has " +
                 std::to_string(fields.size()) + " fields; expected " +
                 std::to_string(expected_header.size()));
        }
        rows.push_back(std::move(fields));
    }
    if (!input.eof()) {
        fail("I/O error while reading " + label);
    }
    return rows;
}

std::vector<Unit> load_units(const fs::path& path) {
    const auto rows = read_tsv(path, kUnitHeader, "units TSV");
    std::vector<Unit> units;
    std::set<std::string> unit_ids;
    std::set<std::string> component_ids;
    std::set<std::tuple<std::string, std::string, std::string, std::string, std::int64_t>>
        route_site_memberships;
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const auto& row = rows[row_index];
        const std::string prefix = "units row " + std::to_string(row_index + 2) + " ";
        Unit unit;
        unit.dataset = row[0];
        unit.chrom = row[1];
        unit.unit_id = row[2];
        unit.component_id = row[3];
        unit.linkage_basis = row[4];
        unit.hp_family = row[5];
        unit.phase_set = row[6];
        validate_text_field(unit.dataset, prefix + "dataset");
        validate_text_field(unit.chrom, prefix + "chrom");
        validate_text_field(unit.unit_id, prefix + "unit_id");
        validate_text_field(unit.component_id, prefix + "component_id");
        if (!unit_ids.insert(unit.unit_id).second) {
            fail("duplicate unit_id: " + unit.unit_id);
        }
        if (!component_ids.insert(unit.component_id).second) {
            fail("duplicate component_id: " + unit.component_id);
        }
        if (unit.hp_family != "1" && unit.hp_family != "2") {
            fail(prefix + "hp_family must be 1 or 2");
        }
        if (unit.phase_set.empty() || unit.phase_set == ".") {
            fail(prefix + "phase_set must be exact and non-missing");
        }
        validate_text_field(unit.phase_set, prefix + "phase_set");
        const std::string expected_basis = "PS_HP" + unit.hp_family;
        if (unit.linkage_basis != expected_basis) {
            fail(prefix + "linkage_basis does not match hp_family");
        }
        unit.threshold = parse_u64(row[7], prefix + "threshold", true);
        if (unit.threshold != 3) {
            fail(prefix + "threshold must equal 3 under the normalized CLI contract");
        }
        const std::uint64_t declared_k = parse_u64(row[8], prefix + "k", true);
        unit.positions = parse_positions(row[9], prefix + "positions1");
        if (declared_k != unit.positions.size()) {
            fail(prefix + "k does not match positions1 length");
        }
        for (const std::int64_t position : unit.positions) {
            const auto key =
                std::make_tuple(unit.dataset, unit.chrom, unit.hp_family, unit.phase_set, position);
            if (!route_site_memberships.insert(key).second) {
                fail(prefix + "duplicates an exact PS x HP route-site membership at position " +
                     std::to_string(position));
            }
        }
        units.push_back(std::move(unit));
    }
    std::sort(units.begin(), units.end(), [](const Unit& left, const Unit& right) {
        return std::tie(left.dataset, left.chrom, left.unit_id) <
               std::tie(right.dataset, right.chrom, right.unit_id);
    });
    return units;
}

void load_constraints(const fs::path& path, std::vector<Unit>* units) {
    const auto rows = read_tsv(path, kConstraintHeader, "constraints TSV");
    std::map<std::string, std::size_t> unit_index;
    for (std::size_t index = 0; index < units->size(); ++index) {
        unit_index.emplace((*units)[index].unit_id, index);
    }
    std::set<std::string> constraint_ids;
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const auto& row = rows[row_index];
        const std::string prefix = "constraints row " + std::to_string(row_index + 2) + " ";
        Constraint constraint;
        constraint.dataset = row[0];
        constraint.chrom = row[1];
        constraint.unit_id = row[2];
        constraint.constraint_id = row[3];
        constraint.hp_family = row[4];
        constraint.phase_set = row[5];
        validate_text_field(constraint.constraint_id, prefix + "constraint_id");
        if (!constraint_ids.insert(constraint.constraint_id).second) {
            fail("duplicate constraint_id: " + constraint.constraint_id);
        }
        const auto found_unit = unit_index.find(constraint.unit_id);
        if (found_unit == unit_index.end()) {
            fail(prefix + "references unknown unit_id: " + constraint.unit_id);
        }
        Unit& unit = (*units)[found_unit->second];
        if (constraint.dataset != unit.dataset) {
            fail(prefix + "dataset mismatch for unit " + unit.unit_id);
        }
        if (constraint.chrom != unit.chrom) {
            fail(prefix + "chrom mismatch for unit " + unit.unit_id);
        }
        if (constraint.hp_family != unit.hp_family) {
            fail(prefix + "hp_family mismatch for unit " + unit.unit_id);
        }
        if (constraint.phase_set != unit.phase_set) {
            fail(prefix + "phase_set mismatch for unit " + unit.unit_id);
        }
        constraint.positions = parse_positions(row[6], prefix + "positions1");
        constraint.call_codes = row[7];
        if (constraint.call_codes.size() != constraint.positions.size()) {
            fail(prefix + "call_codes length does not match positions1 length");
        }
        if (std::any_of(constraint.call_codes.begin(), constraint.call_codes.end(),
                        [](char value) { return value != 'R' && value != 'A'; })) {
            fail(prefix + "call_codes must contain only fixed R/A calls");
        }
        constraint.molecule_weight = parse_nonnegative_bigint(row[8], prefix + "molecule_weight");
        constraint.pattern_count = parse_nonnegative_bigint(row[9], prefix + "pattern_count");
        if (constraint.pattern_count != 1) {
            fail(prefix + "pattern_count must equal 1 under the normalized CLI contract");
        }
        const std::map<std::int64_t, std::size_t> positions_by_index = [&unit]() {
            std::map<std::int64_t, std::size_t> result;
            for (std::size_t index = 0; index < unit.positions.size(); ++index) {
                result.emplace(unit.positions[index], index);
            }
            return result;
        }();
        std::vector<std::size_t> indices;
        indices.reserve(constraint.positions.size());
        for (const std::int64_t position : constraint.positions) {
            const auto found_position = positions_by_index.find(position);
            if (found_position == positions_by_index.end()) {
                fail(prefix + "references unknown unit position: " + std::to_string(position));
            }
            indices.push_back(found_position->second);
        }
        constraint.lo = indices.front();
        constraint.hi = indices.back();
        unit.constraints.push_back(std::move(constraint));
    }
    for (Unit& unit : *units) {
        std::sort(unit.constraints.begin(), unit.constraints.end(),
                  [](const Constraint& left, const Constraint& right) {
                      return left.constraint_id < right.constraint_id;
                  });
    }
}

bool state_is_better(const State& candidate, const std::optional<State>& incumbent) {
    if (!incumbent.has_value()) {
        return true;
    }
    if (candidate.retained_weight != incumbent->retained_weight) {
        return candidate.retained_weight > incumbent->retained_weight;
    }
    if (candidate.retained_pattern_count != incumbent->retained_pattern_count) {
        return candidate.retained_pattern_count > incumbent->retained_pattern_count;
    }
    if (candidate.block_count != incumbent->block_count) {
        return candidate.block_count < incumbent->block_count;
    }
    if (candidate.cut_gap_sum != incumbent->cut_gap_sum) {
        return candidate.cut_gap_sum > incumbent->cut_gap_sum;
    }
    return std::lexicographical_compare(candidate.cuts.begin(), candidate.cuts.end(),
                                        incumbent->cuts.begin(), incumbent->cuts.end());
}

PartitionResult partition_unit(const Unit& unit, std::size_t max_block_size) {
    const std::size_t n_sites = unit.positions.size();
    std::map<std::pair<std::size_t, std::size_t>, Metric> interval_metrics;
    for (const Constraint& constraint : unit.constraints) {
        const std::size_t span_sites = constraint.hi - constraint.lo + 1;
        if (span_sites <= max_block_size) {
            Metric& metric = interval_metrics[{constraint.lo, constraint.hi}];
            metric.weight += constraint.molecule_weight;
            metric.pattern_count += constraint.pattern_count;
        }
    }

    std::vector<std::vector<Metric>> block_gain(
        n_sites, std::vector<Metric>(max_block_size + 1));
    for (std::size_t start = 0; start < n_sites; ++start) {
        Metric running;
        const std::size_t last = std::min(n_sites, start + max_block_size);
        for (std::size_t end_inclusive = start; end_inclusive < last; ++end_inclusive) {
            for (std::size_t lo = start; lo <= end_inclusive; ++lo) {
                const auto found = interval_metrics.find({lo, end_inclusive});
                if (found != interval_metrics.end()) {
                    running.weight += found->second.weight;
                    running.pattern_count += found->second.pattern_count;
                }
            }
            block_gain[start][end_inclusive - start + 1] = running;
        }
    }

    std::vector<std::optional<State>> states(n_sites + 1);
    states[0] = State{};
    for (std::size_t end_exclusive = 1; end_exclusive <= n_sites; ++end_exclusive) {
        std::optional<State> best;
        const std::size_t first_start = end_exclusive > max_block_size
                                            ? end_exclusive - max_block_size
                                            : 0;
        for (std::size_t start = first_start; start < end_exclusive; ++start) {
            if (!states[start].has_value()) {
                fail("internal DP prefix state is missing for unit " + unit.unit_id);
            }
            const Metric& gain = block_gain[start][end_exclusive - start];
            State candidate = *states[start];
            candidate.retained_weight += gain.weight;
            candidate.retained_pattern_count += gain.pattern_count;
            ++candidate.block_count;
            if (start != 0) {
                candidate.cuts.push_back(start);
                candidate.cut_gap_sum += BigInt(unit.positions[start]) - BigInt(unit.positions[start - 1]);
            }
            if (state_is_better(candidate, best)) {
                best = std::move(candidate);
            }
        }
        states[end_exclusive] = std::move(best);
    }
    if (!states[n_sites].has_value()) {
        fail("internal DP produced no terminal state for unit " + unit.unit_id);
    }
    const State& final_state = *states[n_sites];

    std::vector<std::size_t> boundaries;
    boundaries.push_back(0);
    boundaries.insert(boundaries.end(), final_state.cuts.begin(), final_state.cuts.end());
    boundaries.push_back(n_sites);
    std::vector<Block> blocks;
    blocks.reserve(boundaries.size() - 1);
    for (std::size_t index = 0; index + 1 < boundaries.size(); ++index) {
        Block block;
        block.block_index = index;
        block.start_index = boundaries[index];
        block.end_index_exclusive = boundaries[index + 1];
        if (block.end_index_exclusive <= block.start_index ||
            block.end_index_exclusive - block.start_index > max_block_size) {
            fail("internal block-size contract violation for unit " + unit.unit_id);
        }
        blocks.push_back(std::move(block));
    }

    PartitionResult result;
    result.unit = &unit;
    result.cut_indices = final_state.cuts;
    result.cut_gap_sum = final_state.cut_gap_sum;
    result.blocks = std::move(blocks);
    for (const std::size_t cut : result.cut_indices) {
        result.cut_gaps.push_back(BigInt(unit.positions[cut]) - BigInt(unit.positions[cut - 1]));
    }

    for (const Constraint& constraint : unit.constraints) {
        const auto after_hi = std::upper_bound(result.cut_indices.begin(), result.cut_indices.end(), constraint.hi);
        const auto after_lo = std::upper_bound(result.cut_indices.begin(), result.cut_indices.end(), constraint.lo);
        const std::size_t crossed = static_cast<std::size_t>(std::distance(after_lo, after_hi));
        const std::size_t candidate_block = static_cast<std::size_t>(
            std::distance(result.cut_indices.begin(), after_lo));
        if (candidate_block >= result.blocks.size()) {
            fail("internal block lookup overflow for constraint " + constraint.constraint_id);
        }
        const Block& block = result.blocks[candidate_block];
        const bool retained = block.start_index <= constraint.lo &&
                              constraint.hi < block.end_index_exclusive;
        Disposition disposition;
        disposition.constraint = &constraint;
        disposition.span_sites = constraint.hi - constraint.lo + 1;
        disposition.crossed_cut_count = crossed;
        if (retained) {
            disposition.status = kRetained;
            disposition.block_index = candidate_block;
            Block& mutable_block = result.blocks[candidate_block];
            mutable_block.retained_constraint_ids.push_back(constraint.constraint_id);
            mutable_block.retained_weight += constraint.molecule_weight;
            mutable_block.retained_pattern_count += constraint.pattern_count;
            result.retained_weight += constraint.molecule_weight;
            result.retained_pattern_count += constraint.pattern_count;
        } else {
            if (crossed == 0) {
                fail("internal lost constraint crosses no cut: " + constraint.constraint_id);
            }
            disposition.status = disposition.span_sites > max_block_size ? kUnavoidable : kCut;
        }
        result.total_weight += constraint.molecule_weight;
        result.total_pattern_count += constraint.pattern_count;
        result.dispositions.push_back(std::move(disposition));
    }
    result.lost_weight = result.total_weight - result.retained_weight;
    result.lost_pattern_count = result.total_pattern_count - result.retained_pattern_count;
    if (result.retained_weight != final_state.retained_weight ||
        result.retained_pattern_count != final_state.retained_pattern_count) {
        fail("DP objective does not match dispositions for unit " + unit.unit_id);
    }
    for (Block& block : result.blocks) {
        std::sort(block.retained_constraint_ids.begin(), block.retained_constraint_ids.end());
    }
    return result;
}

std::string block_id(const Unit& unit, std::size_t block_index) {
    std::ostringstream stream;
    stream << unit.unit_id << ":B" << std::setw(4) << std::setfill('0') << (block_index + 1);
    return stream.str();
}

std::string block_positions(const Unit& unit, const Block& block) {
    std::vector<std::int64_t> positions(
        unit.positions.begin() + static_cast<std::ptrdiff_t>(block.start_index),
        unit.positions.begin() + static_cast<std::ptrdiff_t>(block.end_index_exclusive));
    return join_integers(positions, ",");
}

void ensure_new_output_dir(const fs::path& output_dir) {
    std::error_code error;
    if (fs::exists(output_dir, error)) {
        fail("output directory already exists; refusing to overwrite: " + output_dir.string());
    }
    if (error) {
        fail("cannot inspect output directory: " + error.message());
    }
    if (!fs::create_directories(output_dir, error) || error) {
        fail("cannot create output directory: " + output_dir.string() + ": " + error.message());
    }
}

std::ofstream open_output(const fs::path& path) {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output) {
        fail("cannot create output file: " + path.string());
    }
    return output;
}

void verify_output_stream(const std::ofstream& output, const fs::path& path) {
    if (!output) {
        fail("I/O error while writing output file: " + path.string());
    }
}

AggregateSummary summarize(const std::vector<PartitionResult>& results) {
    AggregateSummary summary;
    summary.units = results.size();
    for (const PartitionResult& result : results) {
        summary.sites += result.unit->positions.size();
        summary.constraints += result.dispositions.size();
        summary.blocks += result.blocks.size();
        summary.cuts += result.cut_indices.size();
        summary.total_weight += result.total_weight;
        summary.retained_weight += result.retained_weight;
        summary.lost_weight += result.lost_weight;
        summary.total_pattern_count += result.total_pattern_count;
        summary.retained_pattern_count += result.retained_pattern_count;
        summary.lost_pattern_count += result.lost_pattern_count;
        for (const Disposition& disposition : result.dispositions) {
            if (disposition.status == kRetained) {
                ++summary.retained_constraints;
            } else if (disposition.status == kCut) {
                ++summary.cut_constraints;
            } else if (disposition.status == kUnavoidable) {
                ++summary.unavoidable_constraints;
            } else {
                fail("internal unknown disposition status");
            }
        }
    }
    if (summary.constraints != summary.retained_constraints + summary.cut_constraints +
                                   summary.unavoidable_constraints) {
        fail("aggregate constraint disposition conservation failed");
    }
    if (summary.total_weight != summary.retained_weight + summary.lost_weight ||
        summary.total_pattern_count != summary.retained_pattern_count + summary.lost_pattern_count) {
        fail("aggregate objective conservation failed");
    }
    return summary;
}

void write_blocks(const fs::path& path, const std::vector<PartitionResult>& results) {
    std::ofstream output = open_output(path);
    output << "dataset\tchrom\tunit_id\tcomponent_id\tlinkage_basis\thp_family\tphase_set\tthreshold"
              "\tblock_id\tblock_index\tstart1\tend1\tk\tpositions1"
              "\tretained_molecule_weight\tretained_pattern_count"
              "\tstart_index_zero_based\tend_index_exclusive_zero_based\tretained_constraint_count"
              "\tunit_cut_indices_zero_based\tunit_cut_gaps\tunit_cut_gap_sum\tunit_retained_molecule_weight"
              "\tunit_retained_pattern_count\tunit_total_molecule_weight\tunit_total_pattern_count\n";
    for (const PartitionResult& result : results) {
        const Unit& unit = *result.unit;
        for (const Block& block : result.blocks) {
            output << unit.dataset << '\t' << unit.chrom << '\t' << unit.unit_id << '\t'
                   << unit.component_id << '\t' << unit.linkage_basis << '\t' << unit.hp_family << '\t'
                   << unit.phase_set << '\t' << unit.threshold << '\t'
                   << block_id(unit, block.block_index) << '\t' << (block.block_index + 1) << '\t'
                   << unit.positions[block.start_index] << '\t'
                   << unit.positions[block.end_index_exclusive - 1] << '\t'
                   << (block.end_index_exclusive - block.start_index) << '\t'
                   << block_positions(unit, block) << '\t'
                   << block.retained_weight << '\t' << block.retained_pattern_count << '\t'
                   << block.start_index << '\t' << block.end_index_exclusive << '\t'
                   << block.retained_constraint_ids.size() << '\t'
                   << join_integers(result.cut_indices, ",") << '\t'
                   << join_bigints(result.cut_gaps, ",") << '\t' << result.cut_gap_sum << '\t'
                   << result.retained_weight << '\t' << result.retained_pattern_count << '\t'
                   << result.total_weight << '\t' << result.total_pattern_count << '\n';
        }
    }
    output.flush();
    verify_output_stream(output, path);
}

void write_membership(const fs::path& path, const std::vector<PartitionResult>& results) {
    std::ofstream output = open_output(path);
    output << "dataset\tchrom\tunit_id\tcomponent_id\tlinkage_basis\thp_family\tphase_set\tthreshold"
              "\tunit_local_index\tpos1\tblock_id\tblock_index"
              "\tsite_index_in_unit_zero_based\tsite_index_in_block_zero_based\n";
    for (const PartitionResult& result : results) {
        const Unit& unit = *result.unit;
        for (const Block& block : result.blocks) {
            for (std::size_t site_index = block.start_index;
                 site_index < block.end_index_exclusive; ++site_index) {
                output << unit.dataset << '\t' << unit.chrom << '\t' << unit.unit_id << '\t'
                       << unit.component_id << '\t' << unit.linkage_basis << '\t' << unit.hp_family << '\t'
                       << unit.phase_set << '\t' << unit.threshold << '\t' << (site_index + 1) << '\t'
                       << unit.positions[site_index] << '\t' << block_id(unit, block.block_index) << '\t'
                       << (block.block_index + 1) << '\t' << site_index << '\t'
                       << (site_index - block.start_index) << '\n';
            }
        }
    }
    output.flush();
    verify_output_stream(output, path);
}

void write_dispositions(const fs::path& path, const std::vector<PartitionResult>& results) {
    std::ofstream output = open_output(path);
    output << "dataset\tchrom\tunit_id\tconstraint_id\thp_family\tphase_set\tpositions1\tcall_codes"
              "\tmolecule_weight\tpattern_count\tdisposition\tspan_sites\tcrossed_cut_count"
              "\tretained_block_index\tretained_block_id\n";
    for (const PartitionResult& result : results) {
        const Unit& unit = *result.unit;
        for (const Disposition& disposition : result.dispositions) {
            const Constraint& constraint = *disposition.constraint;
            output << constraint.dataset << '\t' << constraint.chrom << '\t' << constraint.unit_id << '\t'
                   << constraint.constraint_id << '\t' << constraint.hp_family << '\t'
                   << constraint.phase_set << '\t' << join_integers(constraint.positions, ",") << '\t'
                   << constraint.call_codes << '\t' << constraint.molecule_weight << '\t'
                   << constraint.pattern_count << '\t' << disposition.status << '\t'
                   << disposition.span_sites << '\t' << disposition.crossed_cut_count << '\t';
            if (disposition.block_index.has_value()) {
                output << (*disposition.block_index + 1) << '\t'
                       << block_id(unit, *disposition.block_index);
            } else {
                output << '\t';
            }
            output << '\n';
        }
    }
    output.flush();
    verify_output_stream(output, path);
}

void write_summary(
    const fs::path& path,
    const AggregateSummary& summary,
    std::size_t max_block_size
) {
    std::ofstream output = open_output(path);
    output << "{\n"
           << "  \"schema_name\": \"intersubmod.exact_ps_k12_cpp_partition_summary\",\n"
           << "  \"schema_version\": \"1.0.0\",\n"
           << "  \"all_pass\": true,\n"
           << "  \"parameters\": {\n"
           << "    \"max_block_size\": " << max_block_size << ",\n"
           << "    \"phase_set_policy\": \"exact_nonmissing_fail_closed\",\n"
           << "    \"analysis_unit\": \"exact_ps_x_hp_x_read_linked_component\",\n"
           << "    \"partition_type\": \"contiguous_nonoverlapping\",\n"
           << "    \"molecule_weight_type\": \"arbitrary_precision_nonnegative_integer\",\n"
           << "    \"objective_order\": [\"max_retained_molecule_weight\", \"max_retained_pattern_count\", "
              "\"min_blocks\", \"max_cut_gap_sum\", \"lexicographically_smaller_cut_tuple\"]\n"
           << "  },\n"
           << "  \"counts\": {\n"
           << "    \"units\": " << summary.units << ",\n"
           << "    \"sites\": " << summary.sites << ",\n"
           << "    \"constraints_total\": " << summary.constraints << ",\n"
           << "    \"blocks\": " << summary.blocks << ",\n"
           << "    \"cuts\": " << summary.cuts << ",\n"
           << "    \"constraints_retained\": " << summary.retained_constraints << ",\n"
           << "    \"constraints_cut\": " << summary.cut_constraints << ",\n"
           << "    \"constraints_unavoidable\": " << summary.unavoidable_constraints << "\n"
           << "  },\n"
           << "  \"weights\": {\n"
           << "    \"total\": \"" << summary.total_weight << "\",\n"
           << "    \"retained\": \"" << summary.retained_weight << "\",\n"
           << "    \"lost\": \"" << summary.lost_weight << "\"\n"
           << "  },\n"
           << "  \"pattern_counts\": {\n"
           << "    \"total\": \"" << summary.total_pattern_count << "\",\n"
           << "    \"retained\": \"" << summary.retained_pattern_count << "\",\n"
           << "    \"lost\": \"" << summary.lost_pattern_count << "\"\n"
           << "  },\n"
           << "  \"checks\": {\n"
           << "    \"exact_nonmissing_ps_hp_units\": true,\n"
           << "    \"cross_ps_violations_zero\": true,\n"
           << "    \"cross_hp_violations_zero\": true,\n"
           << "    \"strictly_increasing_unique_positions\": true,\n"
           << "    \"constraint_positions_within_unit\": true,\n"
           << "    \"unique_unit_and_constraint_ids\": true,\n"
           << "    \"max_block_size_respected\": true,\n"
           << "    \"site_membership_conserved\": true,\n"
           << "    \"constraint_dispositions_mutually_exclusive_and_conserved\": true,\n"
           << "    \"objective_matches_dispositions\": true\n"
           << "  }\n"
           << "}\n";
    output.flush();
    verify_output_stream(output, path);
}

CliOptions parse_cli(int argc, char** argv) {
    CliOptions options;
    bool have_units = false;
    bool have_constraints = false;
    bool have_output = false;
    bool have_max = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            std::cout << "Usage: exact_ps_k12_partition --units units.tsv --constraints constraints.tsv "
                         "--output-dir DIR [--max-block-size 12]\n"
                         "Inputs must be plain TSV. For .gz producer outputs, first run:\n"
                         "  gzip -cd units.tsv.gz > units.tsv\n"
                         "  gzip -cd constraints.tsv.gz > constraints.tsv\n";
            std::exit(0);
        }
        if (index + 1 >= argc) {
            fail("missing value after CLI option: " + argument);
        }
        const std::string value = argv[++index];
        if (argument == "--units") {
            if (have_units) {
                fail("duplicate --units option");
            }
            options.units_path = value;
            have_units = true;
        } else if (argument == "--constraints") {
            if (have_constraints) {
                fail("duplicate --constraints option");
            }
            options.constraints_path = value;
            have_constraints = true;
        } else if (argument == "--output-dir") {
            if (have_output) {
                fail("duplicate --output-dir option");
            }
            options.output_dir = value;
            have_output = true;
        } else if (argument == "--max-block-size") {
            if (have_max) {
                fail("duplicate --max-block-size option");
            }
            const std::uint64_t parsed = parse_u64(value, "--max-block-size", true);
            if (parsed > 12) {
                fail("--max-block-size must not exceed the exact kernel limit 12");
            }
            options.max_block_size = static_cast<std::size_t>(parsed);
            have_max = true;
        } else {
            fail("unknown CLI option: " + argument);
        }
    }
    if (!have_units || !have_constraints || !have_output) {
        fail("--units, --constraints, and --output-dir are required");
    }
    require_plain_tsv(options.units_path, "--units input");
    require_plain_tsv(options.constraints_path, "--constraints input");
    return options;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliOptions options = parse_cli(argc, argv);
        std::vector<Unit> units = load_units(options.units_path);
        load_constraints(options.constraints_path, &units);

        std::vector<PartitionResult> results;
        results.reserve(units.size());
        for (const Unit& unit : units) {
            results.push_back(partition_unit(unit, options.max_block_size));
        }
        const AggregateSummary summary = summarize(results);

        ensure_new_output_dir(options.output_dir);
        write_blocks(options.output_dir / "blocks.tsv", results);
        write_membership(options.output_dir / "membership.tsv", results);
        write_dispositions(options.output_dir / "dispositions.tsv", results);
        write_summary(options.output_dir / "summary.json", summary, options.max_block_size);

        std::cout << "{\"all_pass\":true,\"units\":" << summary.units
                  << ",\"constraints\":" << summary.constraints
                  << ",\"blocks\":" << summary.blocks
                  << ",\"output_dir\":\"" << options.output_dir.string() << "\"}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
