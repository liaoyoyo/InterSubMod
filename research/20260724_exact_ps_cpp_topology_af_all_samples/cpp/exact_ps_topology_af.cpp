// SPDX-License-Identifier: GPL-3.0-only

#include <jansson.h>
#include <openssl/evp.h>

#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
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

#include <sys/types.h>
#include <unistd.h>

#include "longlineage/solver/obligation_bnb.hpp"
#include "longlineage/solver/parent_mapping.hpp"

namespace {

using boost::multiprecision::cpp_int;
using longlineage::solver::ExactFamilyState;
using longlineage::solver::ExactLegalParentChoices;
using longlineage::solver::ExactObjectiveState;
using longlineage::solver::ExactParentMappingSummary;
using longlineage::solver::ExactStructuralCandidate;
using longlineage::solver::ExactStructuralResult;
using longlineage::solver::ExactTopologyProblem;
using longlineage::solver::HypercubeVertex;
using longlineage::solver::ObligationBnbOptions;

constexpr const char* kInputSchemaName = "intersubmod.exact_ps_layered_tree_input";
constexpr const char* kInputSchemaVersion = "1.0.0";
constexpr const char* kUnitSchemaName = "intersubmod.exact_ps_cpp_topology_af.unit";
constexpr const char* kUnitSchemaVersion = "1.0.0";
constexpr const char* kReceiptSchemaName = "intersubmod.exact_ps_cpp_topology_af.receipt";
constexpr const char* kReceiptSchemaVersion = "1.0.0";

class JsonOwner {
  public:
    explicit JsonOwner(json_t* value = nullptr) : value_(value) {}
    ~JsonOwner() {
        if (value_ != nullptr) {
            json_decref(value_);
        }
    }
    JsonOwner(const JsonOwner&) = delete;
    JsonOwner& operator=(const JsonOwner&) = delete;
    JsonOwner(JsonOwner&& other) noexcept : value_(other.value_) { other.value_ = nullptr; }
    JsonOwner& operator=(JsonOwner&& other) noexcept {
        if (this != &other) {
            if (value_ != nullptr) {
                json_decref(value_);
            }
            value_ = other.value_;
            other.value_ = nullptr;
        }
        return *this;
    }
    json_t* get() const { return value_; }
    json_t* release() {
        json_t* value = value_;
        value_ = nullptr;
        return value;
    }

  private:
    json_t* value_;
};

cpp_int abs_integer(cpp_int value) { return value < 0 ? -value : value; }

cpp_int integer_gcd(cpp_int left, cpp_int right) {
    left = abs_integer(left);
    right = abs_integer(right);
    while (right != 0) {
        cpp_int remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

struct Rational {
    cpp_int numerator = 0;
    cpp_int denominator = 1;

    Rational() = default;
    Rational(cpp_int new_numerator, cpp_int new_denominator)
        : numerator(std::move(new_numerator)), denominator(std::move(new_denominator)) {
        normalize();
    }

    void normalize() {
        if (denominator == 0) {
            throw std::runtime_error("rational denominator is zero");
        }
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
        const cpp_int divisor = integer_gcd(numerator, denominator);
        if (divisor != 0) {
            numerator /= divisor;
            denominator /= divisor;
        }
    }

    std::string str() const {
        std::ostringstream output;
        output << numerator << '/' << denominator;
        return output.str();
    }
};

Rational operator+(const Rational& left, const Rational& right) {
    return Rational(left.numerator * right.denominator + right.numerator * left.denominator,
                    left.denominator * right.denominator);
}

Rational operator-(const Rational& left, const Rational& right) {
    return Rational(left.numerator * right.denominator - right.numerator * left.denominator,
                    left.denominator * right.denominator);
}

bool operator==(const Rational& left, const Rational& right) {
    return left.numerator == right.numerator && left.denominator == right.denominator;
}

bool operator<(const Rational& left, const Rational& right) {
    return left.numerator * right.denominator < right.numerator * left.denominator;
}

bool operator>(const Rational& left, const Rational& right) { return right < left; }

cpp_int parse_nonnegative_integer(const std::string& value, const std::string& label) {
    if (value.empty()) {
        throw std::runtime_error(label + " is empty");
    }
    cpp_int result = 0;
    for (char character : value) {
        if (character < '0' || character > '9') {
            throw std::runtime_error(label + " is not a non-negative decimal integer");
        }
        result *= 10;
        result += character - '0';
    }
    return result;
}

std::string decimal(const cpp_int& value) {
    std::ostringstream output;
    output << value;
    return output.str();
}

json_t* json_string_value(const std::string& value) {
    json_t* result = json_stringn(value.data(), value.size());
    if (result == nullptr) {
        throw std::runtime_error("failed to allocate JSON string");
    }
    return result;
}

void set_new(json_t* object, const char* key, json_t* value) {
    if (value == nullptr) {
        throw std::runtime_error(std::string("cannot set null JSON allocation: ") + key);
    }
    // Jansson's *_set_new API steals the value reference, including on failure.
    if (json_object_set_new(object, key, value) != 0) {
        throw std::runtime_error(std::string("failed to set JSON field: ") + key);
    }
}

void set_string(json_t* object, const char* key, const std::string& value) {
    set_new(object, key, json_string_value(value));
}

void set_optional_string(json_t* object, const char* key, const std::optional<std::string>& value) {
    if (value.has_value()) {
        set_string(object, key, *value);
    } else {
        set_new(object, key, json_null());
    }
}

std::string require_string(json_t* object, const char* key, const std::string& label) {
    json_t* value = json_object_get(object, key);
    if (!json_is_string(value)) {
        throw std::runtime_error(label + "." + key + " must be a string");
    }
    return json_string_value(value);
}

json_int_t require_integer(json_t* object, const char* key, const std::string& label) {
    json_t* value = json_object_get(object, key);
    if (!json_is_integer(value)) {
        throw std::runtime_error(label + "." + key + " must be an integer");
    }
    return json_integer_value(value);
}

json_t* optional_object(json_t* object, const char* key, const std::string& label) {
    json_t* value = json_object_get(object, key);
    if (value == nullptr || json_is_null(value)) {
        return nullptr;
    }
    if (!json_is_object(value)) {
        throw std::runtime_error(label + "." + key + " must be an object when present");
    }
    return value;
}

std::string optional_string(json_t* object, const char* key) {
    json_t* value = json_object_get(object, key);
    return json_is_string(value) ? json_string_value(value) : std::string();
}

std::string sha256_bytes(const unsigned char* data, std::size_t size) {
    EVP_MD_CTX* raw_context = EVP_MD_CTX_new();
    if (raw_context == nullptr) {
        throw std::runtime_error("EVP_MD_CTX_new failed");
    }
    struct ContextGuard {
        EVP_MD_CTX* value;
        ~ContextGuard() { EVP_MD_CTX_free(value); }
    } guard{raw_context};

    if (EVP_DigestInit_ex(raw_context, EVP_sha256(), nullptr) != 1 ||
        (size != 0 && EVP_DigestUpdate(raw_context, data, size) != 1)) {
        throw std::runtime_error("SHA-256 initialization/update failed");
    }
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned int digest_size = 0;
    if (EVP_DigestFinal_ex(raw_context, digest, &digest_size) != 1) {
        throw std::runtime_error("SHA-256 finalization failed");
    }
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (unsigned int index = 0; index < digest_size; ++index) {
        output << std::setw(2) << static_cast<unsigned int>(digest[index]);
    }
    return output.str();
}

std::string sha256_text(const std::string& value) {
    return sha256_bytes(reinterpret_cast<const unsigned char*>(value.data()), value.size());
}

std::string sha256_file(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open file for SHA-256: " + path.string());
    }
    EVP_MD_CTX* raw_context = EVP_MD_CTX_new();
    if (raw_context == nullptr) {
        throw std::runtime_error("EVP_MD_CTX_new failed");
    }
    struct ContextGuard {
        EVP_MD_CTX* value;
        ~ContextGuard() { EVP_MD_CTX_free(value); }
    } guard{raw_context};
    if (EVP_DigestInit_ex(raw_context, EVP_sha256(), nullptr) != 1) {
        throw std::runtime_error("SHA-256 initialization failed");
    }
    std::vector<char> buffer(1024 * 1024);
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize count = input.gcount();
        if (count > 0 &&
            EVP_DigestUpdate(raw_context, buffer.data(), static_cast<std::size_t>(count)) != 1) {
            throw std::runtime_error("SHA-256 update failed");
        }
    }
    if (!input.eof()) {
        throw std::runtime_error("failed while reading file for SHA-256: " + path.string());
    }
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned int digest_size = 0;
    if (EVP_DigestFinal_ex(raw_context, digest, &digest_size) != 1) {
        throw std::runtime_error("SHA-256 finalization failed");
    }
    return [&]() {
        std::ostringstream output;
        output << std::hex << std::setfill('0');
        for (unsigned int index = 0; index < digest_size; ++index) {
            output << std::setw(2) << static_cast<unsigned int>(digest[index]);
        }
        return output.str();
    }();
}

struct FileIdentity {
    std::string path;
    std::uintmax_t size_bytes = 0;
    std::string sha256;
};

FileIdentity identify_file(const std::filesystem::path& path) {
    FileIdentity result;
    result.path = std::filesystem::absolute(path).lexically_normal().string();
    result.size_bytes = std::filesystem::file_size(path);
    result.sha256 = sha256_file(path);
    return result;
}

json_t* file_identity_json(const FileIdentity& identity) {
    json_t* object = json_object();
    set_string(object, "path", identity.path);
    set_new(object, "size_bytes", json_integer(static_cast<json_int_t>(identity.size_bytes)));
    set_string(object, "sha256", identity.sha256);
    return object;
}

struct Arguments {
    std::filesystem::path input;
    std::filesystem::path output;
    std::filesystem::path receipt;
    std::size_t maximum_family_size = 0;
    std::uint64_t maximum_search_nodes = 0;
};

std::uint64_t parse_limit(const std::string& value, const std::string& option) {
    if (value.empty() || value.front() == '-') {
        throw std::runtime_error(option + " requires a non-negative integer");
    }
    std::size_t parsed = 0;
    unsigned long long result = 0;
    try {
        result = std::stoull(value, &parsed);
    } catch (const std::exception&) {
        throw std::runtime_error(option + " requires a non-negative integer");
    }
    if (parsed != value.size()) {
        throw std::runtime_error(option + " requires a non-negative integer");
    }
    return static_cast<std::uint64_t>(result);
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments arguments;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        auto next = [&]() -> std::string {
            if (index + 1 >= argc) {
                throw std::runtime_error("missing value for " + option);
            }
            return argv[++index];
        };
        if (option == "--input") {
            arguments.input = next();
        } else if (option == "--output") {
            arguments.output = next();
        } else if (option == "--receipt") {
            arguments.receipt = next();
        } else if (option == "--max-family-size") {
            const std::uint64_t value = parse_limit(next(), option);
            if (value > std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error(option + " exceeds size_t");
            }
            arguments.maximum_family_size = static_cast<std::size_t>(value);
        } else if (option == "--max-search-nodes") {
            arguments.maximum_search_nodes = parse_limit(next(), option);
        } else if (option == "--help" || option == "-h") {
            std::cout
                << "Usage: exact_ps_topology_af --input MLHP.json --output topology.jsonl "
                   "--receipt topology.receipt.json [--max-family-size N] "
                   "[--max-search-nodes N]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("unknown option: " + option);
        }
    }
    if (arguments.input.empty() || arguments.output.empty() || arguments.receipt.empty()) {
        throw std::runtime_error("--input, --output and --receipt are required");
    }
    if (arguments.output == arguments.receipt) {
        throw std::runtime_error("--output and --receipt must be different paths");
    }
    return arguments;
}

struct CoverageRow {
    std::size_t original_index = 0;
    std::int64_t position = 0;
    std::int64_t reference_count = 0;
    std::int64_t alternate_count = 0;
    std::optional<Rational> fraction;
};

enum class CoverageState {
    kComplete,
    kMissing,
    kZeroDenominator,
};

struct ParsedUnit {
    std::string sample;
    std::string region_id;
    std::string unit_id;
    std::string block_id;
    std::string chromosome;
    std::string hp_family;
    std::string phase_set;
    std::string cn;
    bool input_vaf_eligible = false;
    std::vector<std::int64_t> positions;
    std::vector<std::size_t> active_original_indices;
    std::set<HypercubeVertex> full_observed_vertices;
    ExactTopologyProblem problem;
    CoverageState coverage_state = CoverageState::kComplete;
    std::string coverage_message;
    std::vector<CoverageRow> coverage_rows;
    std::vector<Rational> active_read_af;
};

void validate_pattern(const std::string& pattern, std::size_t k, bool partial,
                      const std::string& label) {
    if (pattern.size() != k) {
        throw std::runtime_error(label + " pattern length differs from positions");
    }
    bool has_x = false;
    for (char character : pattern) {
        const bool valid = character == 'R' || character == 'A' || (partial && character == 'X');
        if (!valid) {
            throw std::runtime_error(label + " pattern contains a character outside R/A" +
                                     std::string(partial ? "/X" : ""));
        }
        has_x = has_x || character == 'X';
    }
    if (partial && !has_x) {
        throw std::runtime_error(label + " partial pattern must contain at least one X");
    }
}

std::map<std::string, std::int64_t> parse_pattern_counts(json_t* by_hp, const std::string& family,
                                                        std::size_t k, bool partial,
                                                        const std::string& label) {
    std::map<std::string, std::int64_t> result;
    if (by_hp == nullptr) {
        return result;
    }
    json_t* patterns = json_object_get(by_hp, family.c_str());
    if (patterns == nullptr) {
        return result;
    }
    if (!json_is_object(patterns)) {
        throw std::runtime_error(label + "." + family + " must be an object");
    }
    const char* key = nullptr;
    json_t* value = nullptr;
    json_object_foreach(patterns, key, value) {
        const std::string pattern = key;
        validate_pattern(pattern, k, partial, label);
        if (!json_is_integer(value) || json_integer_value(value) <= 0) {
            throw std::runtime_error(label + "." + family + "." + pattern +
                                     " must be a positive integer count");
        }
        result.emplace(pattern, static_cast<std::int64_t>(json_integer_value(value)));
    }
    return result;
}

HypercubeVertex compact_mask(const std::string& pattern,
                             const std::vector<std::size_t>& active_original_indices) {
    HypercubeVertex result = 0;
    for (std::size_t compact_index = 0; compact_index < active_original_indices.size(); ++compact_index) {
        if (pattern[active_original_indices[compact_index]] == 'A') {
            result = static_cast<HypercubeVertex>(
                result | static_cast<HypercubeVertex>(HypercubeVertex{1} << compact_index));
        }
    }
    return result;
}

std::vector<HypercubeVertex> partial_group(
    const std::string& pattern, const std::vector<std::size_t>& active_original_indices) {
    HypercubeVertex fixed = 0;
    std::vector<std::size_t> unknown;
    for (std::size_t compact_index = 0; compact_index < active_original_indices.size(); ++compact_index) {
        const char character = pattern[active_original_indices[compact_index]];
        if (character == 'A') {
            fixed = static_cast<HypercubeVertex>(
                fixed | static_cast<HypercubeVertex>(HypercubeVertex{1} << compact_index));
        } else if (character == 'X') {
            unknown.push_back(compact_index);
        }
    }
    const std::size_t completions = std::size_t{1} << unknown.size();
    std::vector<HypercubeVertex> result;
    result.reserve(completions);
    for (std::size_t subset = 0; subset < completions; ++subset) {
        HypercubeVertex vertex = fixed;
        for (std::size_t offset = 0; offset < unknown.size(); ++offset) {
            if ((subset & (std::size_t{1} << offset)) != 0) {
                vertex = static_cast<HypercubeVertex>(
                    vertex | static_cast<HypercubeVertex>(HypercubeVertex{1} << unknown[offset]));
            }
        }
        result.push_back(vertex);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

ParsedUnit parse_unit(json_t* group, std::size_t group_index, const std::string& document_sample) {
    const std::string label = "groups[" + std::to_string(group_index) + "]";
    if (!json_is_object(group)) {
        throw std::runtime_error(label + " must be an object");
    }
    ParsedUnit result;
    result.sample = require_string(group, "dataset", label);
    if (result.sample != document_sample) {
        throw std::runtime_error(label + ".dataset differs from document sample");
    }
    result.region_id = require_string(group, "region_id", label);
    result.unit_id = require_string(group, "unit_id", label);
    result.block_id = require_string(group, "block_id", label);
    result.chromosome = require_string(group, "chrom", label);
    result.hp_family = require_string(group, "hp_family", label);
    result.phase_set = require_string(group, "phase_set", label);
    result.cn = optional_string(group, "cn");
    if (result.hp_family != "1" && result.hp_family != "2") {
        throw std::runtime_error(label + ".hp_family must be 1 or 2");
    }
    json_t* vaf_eligible = json_object_get(group, "vaf_eligible");
    if (vaf_eligible != nullptr && !json_is_boolean(vaf_eligible)) {
        throw std::runtime_error(label + ".vaf_eligible must be boolean");
    }
    result.input_vaf_eligible = json_is_true(vaf_eligible);

    json_t* positions = json_object_get(group, "positions");
    if (!json_is_array(positions)) {
        throw std::runtime_error(label + ".positions must be an array");
    }
    const std::size_t k = json_array_size(positions);
    if (k == 0 || k > longlineage::solver::kMaximumExactHypercubeBits) {
        throw std::runtime_error(label + ".positions length must be in [1,12]");
    }
    result.positions.reserve(k);
    for (std::size_t index = 0; index < k; ++index) {
        json_t* value = json_array_get(positions, index);
        if (!json_is_integer(value) || json_integer_value(value) <= 0) {
            throw std::runtime_error(label + ".positions must contain positive integers");
        }
        const std::int64_t position = static_cast<std::int64_t>(json_integer_value(value));
        if (!result.positions.empty() && result.positions.back() >= position) {
            throw std::runtime_error(label + ".positions must be strictly increasing");
        }
        result.positions.push_back(position);
    }
    const json_int_t declared_k = require_integer(group, "n_sSNV", label);
    if (declared_k != static_cast<json_int_t>(k)) {
        throw std::runtime_error(label + ".n_sSNV differs from positions length");
    }

    const auto full = parse_pattern_counts(optional_object(group, "populations_by_hp", label),
                                           result.hp_family, k, false,
                                           label + ".populations_by_hp");
    const auto partial = parse_pattern_counts(optional_object(group, "subread_groups_by_hp", label),
                                              result.hp_family, k, true,
                                              label + ".subread_groups_by_hp");
    if (full.empty() && partial.empty()) {
        throw std::runtime_error(label + " has no supported full or partial pattern for its HP family");
    }

    std::vector<bool> active(k, false);
    for (const auto& [pattern, count] : full) {
        (void)count;
        for (std::size_t index = 0; index < k; ++index) {
            active[index] = active[index] || pattern[index] == 'A';
        }
    }
    for (const auto& [pattern, count] : partial) {
        (void)count;
        for (std::size_t index = 0; index < k; ++index) {
            active[index] = active[index] || pattern[index] == 'A';
        }
    }
    for (std::size_t index = 0; index < k; ++index) {
        if (active[index]) {
            result.active_original_indices.push_back(index);
        }
    }
    result.problem.bit_count = result.active_original_indices.size();
    for (const auto& [pattern, count] : full) {
        (void)count;
        const HypercubeVertex vertex =
            compact_mask(pattern, result.active_original_indices);
        result.problem.mandatory_vertices.push_back(vertex);
        result.full_observed_vertices.insert(vertex);
    }
    std::sort(result.problem.mandatory_vertices.begin(), result.problem.mandatory_vertices.end());
    result.problem.mandatory_vertices.erase(
        std::unique(result.problem.mandatory_vertices.begin(), result.problem.mandatory_vertices.end()),
        result.problem.mandatory_vertices.end());
    for (const auto& [pattern, count] : partial) {
        (void)count;
        result.problem.terminal_groups.push_back(partial_group(pattern, result.active_original_indices));
    }

    json_t* coverage_by_hp = optional_object(group, "col_coverage_by_hp", label);
    json_t* family_coverage =
        coverage_by_hp == nullptr ? nullptr : json_object_get(coverage_by_hp, result.hp_family.c_str());
    if (family_coverage != nullptr && !json_is_object(family_coverage)) {
        throw std::runtime_error(label + ".col_coverage_by_hp." + result.hp_family +
                                 " must be an object");
    }
    result.coverage_rows.reserve(k);
    std::vector<std::optional<Rational>> all_fractions(k);
    for (std::size_t index = 0; index < k; ++index) {
        CoverageRow row;
        row.original_index = index;
        row.position = result.positions[index];
        json_t* counts = family_coverage == nullptr
                             ? nullptr
                             : json_object_get(family_coverage, std::to_string(row.position).c_str());
        if (counts == nullptr) {
            if (result.coverage_state == CoverageState::kComplete) {
                result.coverage_state = CoverageState::kMissing;
                result.coverage_message = "missing coverage at " + std::to_string(row.position);
            }
            result.coverage_rows.push_back(row);
            continue;
        }
        if (!json_is_array(counts) || json_array_size(counts) < 2 ||
            !json_is_integer(json_array_get(counts, 0)) ||
            !json_is_integer(json_array_get(counts, 1))) {
            throw std::runtime_error(label + ".col_coverage_by_hp." + result.hp_family + "." +
                                     std::to_string(row.position) +
                                     " must be an array with integer [REF,ALT]");
        }
        row.reference_count = static_cast<std::int64_t>(json_integer_value(json_array_get(counts, 0)));
        row.alternate_count = static_cast<std::int64_t>(json_integer_value(json_array_get(counts, 1)));
        if (row.reference_count < 0 || row.alternate_count < 0) {
            throw std::runtime_error(label + " coverage counts must be non-negative");
        }
        const std::int64_t denominator = row.reference_count + row.alternate_count;
        if (denominator <= 0) {
            if (result.coverage_state == CoverageState::kComplete) {
                result.coverage_state = CoverageState::kZeroDenominator;
                result.coverage_message = "zero denominator at " + std::to_string(row.position);
            }
        } else {
            row.fraction = Rational(row.alternate_count, denominator);
            all_fractions[index] = row.fraction;
        }
        result.coverage_rows.push_back(row);
    }
    if (result.coverage_state == CoverageState::kComplete) {
        result.active_read_af.reserve(result.active_original_indices.size());
        for (std::size_t index : result.active_original_indices) {
            if (!all_fractions[index].has_value()) {
                throw std::runtime_error(label + " internal coverage-state mismatch");
            }
            result.active_read_af.push_back(*all_fractions[index]);
        }
    }
    return result;
}

bool rooted_three_gamete_compatible(const std::vector<HypercubeVertex>& vertices,
                                    std::size_t bit_count) {
    for (std::size_t left = 0; left < bit_count; ++left) {
        for (std::size_t right = left + 1; right < bit_count; ++right) {
            bool observed_11 = false;
            bool observed_10 = false;
            bool observed_01 = false;
            const HypercubeVertex left_mask =
                static_cast<HypercubeVertex>(HypercubeVertex{1} << left);
            const HypercubeVertex right_mask =
                static_cast<HypercubeVertex>(HypercubeVertex{1} << right);
            for (HypercubeVertex vertex : vertices) {
                const bool left_present = (vertex & left_mask) != 0;
                const bool right_present = (vertex & right_mask) != 0;
                observed_11 = observed_11 || (left_present && right_present);
                observed_10 = observed_10 || (left_present && !right_present);
                observed_01 = observed_01 || (!left_present && right_present);
            }
            if (observed_11 && observed_10 && observed_01) {
                return false;
            }
        }
    }
    return true;
}

std::size_t one_bit_index(HypercubeVertex value) {
    if (value == 0 || (value & (value - 1)) != 0) {
        throw std::runtime_error("edge is not an upward Hamming-1 transition");
    }
    std::size_t result = 0;
    while ((value & 1U) == 0) {
        value = static_cast<HypercubeVertex>(value >> 1U);
        ++result;
    }
    return result;
}

Rational edge_score(HypercubeVertex parent, HypercubeVertex child,
                    const std::vector<Rational>& read_af) {
    if ((parent & child) != parent) {
        throw std::runtime_error("parent is not a subset of child");
    }
    const std::size_t acquired = one_bit_index(static_cast<HypercubeVertex>(child ^ parent));
    if (acquired >= read_af.size()) {
        throw std::runtime_error("acquired bit lies outside read-AF vector");
    }
    Rational result;
    for (std::size_t bit = 0; bit < read_af.size(); ++bit) {
        if ((parent & static_cast<HypercubeVertex>(HypercubeVertex{1} << bit)) != 0) {
            result = result + (read_af[bit] - read_af[acquired]);
        }
    }
    return result;
}

struct RankedCandidate {
    Rational score;
    cpp_int tie_count = 0;
    std::vector<std::pair<HypercubeVertex, HypercubeVertex>> representative_edges;
    std::vector<HypercubeVertex> vertices;
};

RankedCandidate rank_candidate(const ExactStructuralCandidate& candidate, std::size_t bit_count,
                               const std::vector<Rational>& read_af) {
    const ExactParentMappingSummary summary =
        longlineage::solver::summarize_exact_parent_mappings(bit_count, candidate.vertices);
    if (!summary.valid || !candidate.parent_mapping.valid ||
        summary.tree_count != candidate.parent_mapping.tree_count ||
        summary.legal_parent_count != candidate.parent_mapping.legal_parent_count) {
        throw std::runtime_error("pinned parent factorization disagrees with structural candidate");
    }
    RankedCandidate result;
    result.tie_count = 1;
    result.vertices = candidate.vertices;
    result.representative_edges.reserve(summary.legal_parents.size());
    for (const ExactLegalParentChoices& choices : summary.legal_parents) {
        std::optional<Rational> best;
        std::vector<HypercubeVertex> best_parents;
        for (HypercubeVertex parent : choices.parents) {
            const Rational score = edge_score(parent, choices.vertex, read_af);
            if (!best.has_value() || score > *best) {
                best = score;
                best_parents = {parent};
            } else if (score == *best) {
                best_parents.push_back(parent);
            }
        }
        if (!best.has_value() || best_parents.empty()) {
            throw std::runtime_error("selected vertex lacks a legal parent");
        }
        result.score = result.score + *best;
        result.tie_count *= best_parents.size();
        result.representative_edges.emplace_back(best_parents.front(), choices.vertex);
    }
    return result;
}

std::string original_label(HypercubeVertex vertex, const ParsedUnit& unit) {
    if (vertex == 0) {
        return "ROOT";
    }
    std::string result(unit.positions.size(), 'R');
    for (std::size_t compact_index = 0;
         compact_index < unit.active_original_indices.size(); ++compact_index) {
        if ((vertex & static_cast<HypercubeVertex>(HypercubeVertex{1} << compact_index)) != 0) {
            result[unit.active_original_indices[compact_index]] = 'A';
        }
    }
    return unit.full_observed_vertices.count(vertex) != 0 ? result : "H_" + result;
}

std::string canonical_family_text(const ExactStructuralResult& structural,
                                  std::size_t active_bit_count) {
    std::vector<std::vector<HypercubeVertex>> canonical_family;
    canonical_family.reserve(structural.minimum_family.size());
    for (const ExactStructuralCandidate& candidate : structural.minimum_family) {
        canonical_family.push_back(candidate.vertices);
    }
    std::sort(canonical_family.begin(), canonical_family.end());
    std::ostringstream output;
    output << "active_bit_count=" << active_bit_count << '\n';
    for (const std::vector<HypercubeVertex>& vertices : canonical_family) {
        for (std::size_t index = 0; index < vertices.size(); ++index) {
            if (index != 0) {
                output << ',';
            }
            output << vertices[index];
        }
        output << '\n';
    }
    return output.str();
}

json_t* coverage_rows_json(const ParsedUnit& unit) {
    json_t* rows = json_array();
    for (const CoverageRow& row : unit.coverage_rows) {
        json_t* value = json_object();
        set_new(value, "site_index", json_integer(static_cast<json_int_t>(row.original_index + 1)));
        set_new(value, "position", json_integer(static_cast<json_int_t>(row.position)));
        set_new(value, "ref_reads", json_integer(static_cast<json_int_t>(row.reference_count)));
        set_new(value, "alt_reads", json_integer(static_cast<json_int_t>(row.alternate_count)));
        set_optional_string(value, "fraction",
                            row.fraction.has_value()
                                ? std::optional<std::string>(row.fraction->str())
                                : std::nullopt);
        if (json_array_append_new(rows, value) != 0) {
            json_decref(rows);
            throw std::runtime_error("failed to append coverage row");
        }
    }
    return rows;
}

json_t* active_indices_json(const ParsedUnit& unit) {
    json_t* values = json_array();
    for (std::size_t index : unit.active_original_indices) {
        if (json_array_append_new(values, json_integer(static_cast<json_int_t>(index))) != 0) {
            json_decref(values);
            throw std::runtime_error("failed to append active index");
        }
    }
    return values;
}

json_t* active_positions_json(const ParsedUnit& unit) {
    json_t* values = json_array();
    for (std::size_t index : unit.active_original_indices) {
        if (json_array_append_new(values,
                                  json_integer(static_cast<json_int_t>(unit.positions[index]))) != 0) {
            json_decref(values);
            throw std::runtime_error("failed to append active position");
        }
    }
    return values;
}

json_t* edges_json(const std::vector<std::pair<HypercubeVertex, HypercubeVertex>>& edges,
                   const ParsedUnit& unit) {
    json_t* result = json_array();
    for (const auto& [parent, child] : edges) {
        json_t* edge = json_object();
        set_new(edge, "parent_vertex", json_integer(parent));
        set_new(edge, "child_vertex", json_integer(child));
        set_string(edge, "parent_label", original_label(parent, unit));
        set_string(edge, "child_label", original_label(child, unit));
        const std::size_t acquired =
            one_bit_index(static_cast<HypercubeVertex>(parent ^ child));
        set_new(edge, "acquired_active_bit", json_integer(static_cast<json_int_t>(acquired)));
        const std::size_t original_index = unit.active_original_indices[acquired];
        set_new(edge, "acquired_site_index",
                json_integer(static_cast<json_int_t>(original_index + 1)));
        set_new(edge, "acquired_position",
                json_integer(static_cast<json_int_t>(unit.positions[original_index])));
        set_string(edge, "edge_score_fraction",
                   edge_score(parent, child, unit.active_read_af).str());
        if (json_array_append_new(result, edge) != 0) {
            json_decref(result);
            throw std::runtime_error("failed to append representative edge");
        }
    }
    return result;
}

std::size_t vertex_popcount(HypercubeVertex vertex) {
    std::size_t result = 0;
    while (vertex != 0) {
        vertex = static_cast<HypercubeVertex>(vertex & (vertex - 1));
        ++result;
    }
    return result;
}

std::string unlabeled_shape_at(
    HypercubeVertex vertex,
    const std::map<HypercubeVertex, std::vector<HypercubeVertex>>& children) {
    std::vector<std::string> child_shapes;
    const auto iterator = children.find(vertex);
    if (iterator != children.end()) {
        child_shapes.reserve(iterator->second.size());
        for (HypercubeVertex child : iterator->second) {
            child_shapes.push_back(unlabeled_shape_at(child, children));
        }
    }
    std::sort(child_shapes.begin(), child_shapes.end());
    std::string result = "(";
    for (const std::string& child_shape : child_shapes) {
        result += child_shape;
    }
    result += ")";
    return result;
}

json_t* representative_morphology_json(
    const RankedCandidate& representative) {
    std::map<HypercubeVertex, std::vector<HypercubeVertex>> children;
    bool has_sister = false;
    bool has_direct = false;
    for (HypercubeVertex vertex : representative.vertices) {
        has_direct = has_direct || vertex_popcount(vertex) >= 2;
    }
    for (const auto& [parent, child] : representative.representative_edges) {
        children[parent].push_back(child);
    }
    for (auto& [parent, child_values] : children) {
        (void)parent;
        std::sort(child_values.begin(), child_values.end());
        has_sister = has_sister || child_values.size() >= 2;
    }
    const std::string signature = unlabeled_shape_at(0, children);
    json_t* result = json_object();
    set_new(result, "has_direct", json_boolean(has_direct));
    set_new(result, "has_sister", json_boolean(has_sister));
    set_string(result, "shape_signature", signature);
    set_string(result, "shape_sha256", sha256_text(signature));
    set_string(result, "scope",
               "one deterministic globally-best representative; not an exact top-shape tie census");
    return result;
}

json_t* base_unit_json(json_t* group, std::size_t group_index, const std::string& sample) {
    json_t* output = json_object();
    set_string(output, "schema_name", kUnitSchemaName);
    set_string(output, "schema_version", kUnitSchemaVersion);
    set_string(output, "sample", sample);
    set_new(output, "group_index", json_integer(static_cast<json_int_t>(group_index)));
    for (const char* key : {"region_id", "unit_id", "block_id", "chrom", "hp_family", "phase_set"}) {
        json_t* value = json_is_object(group) ? json_object_get(group, key) : nullptr;
        if (json_is_string(value)) {
            set_new(output, key, json_incref(value));
        } else {
            set_new(output, key, json_null());
        }
    }
    return output;
}

struct UnitOutcome {
    JsonOwner record;
    std::string unit_status;
    std::string objective_status;
    std::string family_status;
    std::string read_af_status;
    bool mutation_bearing = false;
    bool objective_certified = false;
    bool family_complete = false;
    bool ranking_complete = false;
};

UnitOutcome malformed_unit(json_t* group, std::size_t group_index, const std::string& sample,
                           const std::string& message) {
    JsonOwner record(base_unit_json(group, group_index, sample));
    set_string(record.get(), "unit_status", "malformed_unit");
    set_string(record.get(), "objective_status", "ABSTAIN_NOT_IDENTIFIABLE");
    set_string(record.get(), "family_status", "ABSTAIN_NOT_IDENTIFIABLE");
    set_string(record.get(), "reason", "MALFORMED_UNIT");
    set_string(record.get(), "message", message);
    set_new(record.get(), "objective_h", json_null());
    set_new(record.get(), "active_bit_count", json_null());
    set_new(record.get(), "minimum_vertex_set_count", json_null());
    set_new(record.get(), "minimum_family_sha256", json_null());
    set_new(record.get(), "total_tree_count", json_null());
    set_new(record.get(), "has_recfree", json_null());
    set_new(record.get(), "recurrence_required", json_null());
    set_string(record.get(), "read_af_status", "not_evaluable_malformed");
    set_new(record.get(), "best_score_fraction", json_null());
    set_new(record.get(), "best_tree_tie_count", json_null());
    set_new(record.get(), "best_tree_unique", json_null());
    set_new(record.get(), "representative_best_edges", json_array());
    set_new(record.get(), "search_nodes", json_integer(0));
    set_new(record.get(), "solver_elapsed_microseconds", json_integer(0));
    return {std::move(record), "malformed_unit", "ABSTAIN_NOT_IDENTIFIABLE",
            "ABSTAIN_NOT_IDENTIFIABLE", "not_evaluable_malformed", false, false, false, false};
}

UnitOutcome analyze_unit(json_t* group, std::size_t group_index, const std::string& sample,
                         const Arguments& arguments) {
    ParsedUnit unit;
    try {
        unit = parse_unit(group, group_index, sample);
    } catch (const std::exception& error) {
        return malformed_unit(group, group_index, sample, error.what());
    }

    JsonOwner record(base_unit_json(group, group_index, sample));
    set_new(record.get(), "original_bit_count",
            json_integer(static_cast<json_int_t>(unit.positions.size())));
    set_new(record.get(), "active_bit_count",
            json_integer(static_cast<json_int_t>(unit.active_original_indices.size())));
    set_new(record.get(), "active_original_indices", active_indices_json(unit));
    set_new(record.get(), "active_positions", active_positions_json(unit));
    set_new(record.get(), "input_vaf_eligible", json_boolean(unit.input_vaf_eligible));
    set_string(record.get(), "af_basis",
               "family-specific ALT/(REF+ALT) from exact-PS projected col_coverage_by_hp");
    set_new(record.get(), "af_coverage", coverage_rows_json(unit));
    set_string(record.get(), "model",
               "recurrence-allowed unit-cost rooted directed Boolean-hypercube group arborescence");

    ObligationBnbOptions options;
    options.maximum_complete_family_size = arguments.maximum_family_size;
    options.maximum_search_nodes = arguments.maximum_search_nodes;
    const auto solver_start = std::chrono::steady_clock::now();
    const ExactStructuralResult structural =
        longlineage::solver::solve_obligation_bnb(unit.problem, options);
    const auto solver_finish = std::chrono::steady_clock::now();
    const auto solver_microseconds =
        std::chrono::duration_cast<std::chrono::microseconds>(solver_finish - solver_start)
            .count();
    const std::string objective_status = longlineage::solver::to_string(structural.objective_state);
    const std::string family_status = longlineage::solver::to_string(structural.family_state);
    const std::string reason = longlineage::solver::to_string(structural.reason);
    const bool objective_certified =
        structural.objective_state == ExactObjectiveState::kObjectiveCertified;
    const bool family_complete = structural.family_state == ExactFamilyState::kFamilyComplete;
    const bool mutation_bearing = !unit.active_original_indices.empty();

    set_string(record.get(), "objective_status", objective_status);
    set_string(record.get(), "family_status", family_status);
    set_string(record.get(), "reason", reason);
    set_string(record.get(), "message", structural.message);
    if (structural.objective_h.has_value()) {
        set_new(record.get(), "objective_h",
                json_integer(static_cast<json_int_t>(*structural.objective_h)));
    } else {
        set_new(record.get(), "objective_h", json_null());
    }
    set_new(record.get(), "search_nodes",
            json_integer(static_cast<json_int_t>(structural.search_nodes)));
    set_new(record.get(), "solver_elapsed_microseconds",
            json_integer(static_cast<json_int_t>(solver_microseconds)));
    set_new(record.get(), "objective_search_exhausted",
            json_boolean(structural.objective_search_exhausted));
    set_new(record.get(), "family_enumeration_exhausted",
            json_boolean(structural.family_enumeration_exhausted));

    if (!family_complete) {
        set_string(record.get(), "unit_status", "family_incomplete");
        set_new(record.get(), "minimum_vertex_set_count", json_null());
        set_new(record.get(), "minimum_family_sha256", json_null());
        set_new(record.get(), "total_tree_count", json_null());
        set_new(record.get(), "has_recfree", json_null());
        set_new(record.get(), "recurrence_required", json_null());
        set_string(record.get(), "read_af_status", "not_evaluable_family_incomplete");
        set_new(record.get(), "best_score_fraction", json_null());
        set_new(record.get(), "best_tree_tie_count", json_null());
        set_new(record.get(), "best_tree_unique", json_null());
        set_new(record.get(), "representative_best_edges", json_array());
        return {std::move(record), "family_incomplete", objective_status, family_status,
                "not_evaluable_family_incomplete", mutation_bearing, objective_certified,
                false, false};
    }

    set_new(record.get(), "minimum_vertex_set_count",
            json_integer(static_cast<json_int_t>(structural.minimum_family.size())));
    const std::string family_text =
        canonical_family_text(structural, unit.active_original_indices.size());
    set_string(record.get(), "minimum_family_sha256", sha256_text(family_text));
    set_string(record.get(), "minimum_family_digest_format",
               "SHA256 of 'active_bit_count=k\\n' followed by one sorted decimal vertex set per line");

    cpp_int total_tree_count = 0;
    std::size_t compatible_count = 0;
    std::size_t incompatible_count = 0;
    for (const ExactStructuralCandidate& candidate : structural.minimum_family) {
        const ExactParentMappingSummary summary =
            longlineage::solver::summarize_exact_parent_mappings(
                unit.active_original_indices.size(), candidate.vertices);
        if (!summary.valid || summary.tree_count != candidate.parent_mapping.tree_count) {
            throw std::runtime_error("parent-factorization differential mismatch");
        }
        total_tree_count += parse_nonnegative_integer(summary.tree_count, "tree_count");
        if (rooted_three_gamete_compatible(candidate.vertices,
                                           unit.active_original_indices.size())) {
            ++compatible_count;
        } else {
            ++incompatible_count;
        }
    }
    const bool has_recfree = compatible_count > 0;
    const bool recurrence_required = !structural.minimum_family.empty() && !has_recfree;
    set_string(record.get(), "total_tree_count", decimal(total_tree_count));
    set_new(record.get(), "recurrence_compatible_vertex_set_count",
            json_integer(static_cast<json_int_t>(compatible_count)));
    set_new(record.get(), "recurrence_incompatible_vertex_set_count",
            json_integer(static_cast<json_int_t>(incompatible_count)));
    set_new(record.get(), "has_recfree", json_boolean(has_recfree));
    set_new(record.get(), "recurrence_required", json_boolean(recurrence_required));

    if (!mutation_bearing) {
        set_string(record.get(), "unit_status", "no_active_alt");
        set_string(record.get(), "read_af_status", "not_applicable_no_active_alt");
        set_new(record.get(), "best_score_fraction", json_null());
        set_new(record.get(), "best_tree_tie_count", json_null());
        set_new(record.get(), "best_tree_unique", json_null());
        set_new(record.get(), "representative_best_edges", json_array());
        return {std::move(record), "no_active_alt", objective_status, family_status,
                "not_applicable_no_active_alt", false, objective_certified, true, false};
    }
    if (recurrence_required) {
        set_string(record.get(), "unit_status", "recurrence_required");
        set_string(record.get(), "read_af_status", "recurrence_not_evaluable");
        set_new(record.get(), "best_score_fraction", json_null());
        set_new(record.get(), "best_tree_tie_count", json_null());
        set_new(record.get(), "best_tree_unique", json_null());
        set_new(record.get(), "representative_best_edges", json_array());
        return {std::move(record), "recurrence_required", objective_status, family_status,
                "recurrence_not_evaluable", true, objective_certified, true, false};
    }
    if (unit.coverage_state != CoverageState::kComplete) {
        const std::string coverage_status =
            unit.coverage_state == CoverageState::kMissing ? "missing_read_af" : "zero_denominator";
        set_string(record.get(), "unit_status", coverage_status);
        set_string(record.get(), "read_af_status", coverage_status);
        set_string(record.get(), "read_af_message", unit.coverage_message);
        set_new(record.get(), "best_score_fraction", json_null());
        set_new(record.get(), "best_tree_tie_count", json_null());
        set_new(record.get(), "best_tree_unique", json_null());
        set_new(record.get(), "representative_best_edges", json_array());
        return {std::move(record), coverage_status, objective_status, family_status,
                coverage_status, true, objective_certified, true, false};
    }

    std::optional<Rational> global_best;
    cpp_int global_tie_count = 0;
    std::size_t best_vertex_set_count = 0;
    std::optional<RankedCandidate> representative;
    for (const ExactStructuralCandidate& candidate : structural.minimum_family) {
        RankedCandidate ranked =
            rank_candidate(candidate, unit.active_original_indices.size(), unit.active_read_af);
        if (!global_best.has_value() || ranked.score > *global_best) {
            global_best = ranked.score;
            global_tie_count = ranked.tie_count;
            best_vertex_set_count = 1;
            representative = std::move(ranked);
        } else if (ranked.score == *global_best) {
            global_tie_count += ranked.tie_count;
            ++best_vertex_set_count;
        }
    }
    if (!global_best.has_value() || !representative.has_value()) {
        throw std::runtime_error("complete structural family produced no rankable candidate");
    }
    set_string(record.get(), "unit_status", "ranked");
    set_string(record.get(), "read_af_status", "ranked_complete");
    set_string(record.get(), "best_score_fraction", global_best->str());
    set_string(record.get(), "best_tree_tie_count", decimal(global_tie_count));
    set_new(record.get(), "best_tree_unique", json_boolean(global_tie_count == 1));
    set_new(record.get(), "best_vertex_set_count",
            json_integer(static_cast<json_int_t>(best_vertex_set_count)));
    json_t* representative_vertices = json_array();
    for (HypercubeVertex vertex : representative->vertices) {
        json_t* value = json_object();
        set_new(value, "vertex", json_integer(vertex));
        set_string(value, "label", original_label(vertex, unit));
        if (json_array_append_new(representative_vertices, value) != 0) {
            json_decref(representative_vertices);
            throw std::runtime_error("failed to append representative vertex");
        }
    }
    set_new(record.get(), "representative_best_vertices", representative_vertices);
    set_new(record.get(), "representative_best_edges",
            edges_json(representative->representative_edges, unit));
    set_new(record.get(), "representative_best_morphology",
            representative_morphology_json(*representative));
    return {std::move(record), "ranked", objective_status, family_status, "ranked_complete",
            true, objective_certified, true, true};
}

std::string utc_now() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t time = std::chrono::system_clock::to_time_t(now);
    std::tm value{};
    gmtime_r(&time, &value);
    std::ostringstream output;
    output << std::put_time(&value, "%Y-%m-%dT%H:%M:%SZ");
    return output.str();
}

std::filesystem::path temporary_sibling(const std::filesystem::path& target) {
    return target.parent_path() /
           ("." + target.filename().string() + ".pending." + std::to_string(getpid()));
}

void write_json_line(std::ofstream& output, json_t* value) {
    char* payload = json_dumps(value, JSON_COMPACT | JSON_SORT_KEYS);
    if (payload == nullptr) {
        throw std::runtime_error("failed to serialize unit JSON");
    }
    output << payload << '\n';
    std::free(payload);
    if (!output) {
        throw std::runtime_error("failed while writing JSONL output");
    }
}

void write_json_file(const std::filesystem::path& path, json_t* value) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("cannot open JSON output: " + path.string());
    }
    char* payload = json_dumps(value, JSON_COMPACT | JSON_SORT_KEYS);
    if (payload == nullptr) {
        throw std::runtime_error("failed to serialize JSON");
    }
    output << payload << '\n';
    std::free(payload);
    output.flush();
    if (!output) {
        throw std::runtime_error("failed while writing JSON output: " + path.string());
    }
}

json_t* census_json(const std::map<std::string, std::uint64_t>& census) {
    json_t* result = json_object();
    for (const auto& [key, value] : census) {
        set_new(result, key.c_str(), json_integer(static_cast<json_int_t>(value)));
    }
    return result;
}

struct RunCounters {
    std::uint64_t total_units = 0;
    std::uint64_t malformed_units = 0;
    std::uint64_t mutation_bearing_units = 0;
    std::uint64_t objective_certified_units = 0;
    std::uint64_t family_complete_mutation_bearing_units = 0;
    std::uint64_t ranking_complete_count = 0;
    std::map<std::string, std::uint64_t> unit_status;
    std::map<std::string, std::uint64_t> objective_status;
    std::map<std::string, std::uint64_t> family_status;
    std::map<std::string, std::uint64_t> read_af_status;
};

json_t* build_receipt(const Arguments& arguments, const FileIdentity& input_identity,
                      const FileIdentity& output_identity, const RunCounters& counters,
                      const std::string& started_at, const std::string& finished_at,
                      double elapsed_seconds) {
    JsonOwner receipt(json_object());
    set_string(receipt.get(), "schema_name", kReceiptSchemaName);
    set_string(receipt.get(), "schema_version", kReceiptSchemaVersion);
    const bool all_units_well_formed = counters.malformed_units == 0;
    set_new(receipt.get(), "all_pass", json_boolean(all_units_well_formed));
    set_new(receipt.get(), "input", file_identity_json(input_identity));
    set_new(receipt.get(), "output", file_identity_json(output_identity));
    set_new(receipt.get(), "all_units_well_formed",
            json_boolean(all_units_well_formed));
    set_new(receipt.get(), "all_units_objective_certified",
            json_boolean(counters.objective_certified_units == counters.total_units));
    set_new(receipt.get(), "all_mutation_bearing_families_complete",
            json_boolean(all_units_well_formed &&
                         counters.family_complete_mutation_bearing_units ==
                             counters.mutation_bearing_units));
    set_new(receipt.get(), "ranking_complete_count",
            json_integer(static_cast<json_int_t>(counters.ranking_complete_count)));

    json_t* counts = json_object();
    set_new(counts, "total_units", json_integer(static_cast<json_int_t>(counters.total_units)));
    set_new(counts, "malformed_units",
            json_integer(static_cast<json_int_t>(counters.malformed_units)));
    set_new(counts, "mutation_bearing_units",
            json_integer(static_cast<json_int_t>(counters.mutation_bearing_units)));
    set_new(counts, "objective_certified_units",
            json_integer(static_cast<json_int_t>(counters.objective_certified_units)));
    set_new(counts, "family_complete_mutation_bearing_units",
            json_integer(static_cast<json_int_t>(
                counters.family_complete_mutation_bearing_units)));
    set_new(counts, "ranking_complete_count",
            json_integer(static_cast<json_int_t>(counters.ranking_complete_count)));
    set_new(receipt.get(), "counts", counts);

    json_t* status = json_object();
    set_new(status, "unit_status", census_json(counters.unit_status));
    set_new(status, "objective_status", census_json(counters.objective_status));
    set_new(status, "family_status", census_json(counters.family_status));
    set_new(status, "read_af_status", census_json(counters.read_af_status));
    set_new(receipt.get(), "status_census", status);

    json_t* parameters = json_object();
    set_new(parameters, "max_family_size",
            json_integer(static_cast<json_int_t>(arguments.maximum_family_size)));
    set_new(parameters, "max_search_nodes",
            json_integer(static_cast<json_int_t>(arguments.maximum_search_nodes)));
    set_string(parameters, "model",
               "recurrence-allowed exact obligation B&B, k_active<=12");
    set_string(parameters, "af_score",
               "sum over edges and parent-present bits of AF(ancestor)-AF(acquired)");
    set_string(parameters, "parent_optimization",
               "independent per child for a fixed vertex set; exact cpp_int tie counts");
    set_string(parameters, "mixed_recurrence_scope",
               "rank every recurrence-allowed minimum candidate when at least one "
               "rooted-three-gamete-compatible minimum vertex set exists; abstain only "
               "when every minimum vertex set requires recurrence");
    set_new(receipt.get(), "parameters", parameters);

    json_t* runtime = json_object();
    set_string(runtime, "started_at_utc", started_at);
    set_string(runtime, "finished_at_utc", finished_at);
    set_new(runtime, "elapsed_seconds", json_real(elapsed_seconds));
    set_new(receipt.get(), "runtime", runtime);
    return receipt.release();
}

int run(const Arguments& arguments) {
    if (!std::filesystem::is_regular_file(arguments.input)) {
        throw std::runtime_error("input is not a regular file: " + arguments.input.string());
    }
    if (std::filesystem::exists(arguments.output) ||
        std::filesystem::exists(arguments.receipt)) {
        throw std::runtime_error("output and receipt targets must not already exist");
    }
    std::filesystem::create_directories(arguments.output.parent_path());
    std::filesystem::create_directories(arguments.receipt.parent_path());

    json_error_t error{};
    JsonOwner document(json_load_file(arguments.input.c_str(), JSON_REJECT_DUPLICATES, &error));
    if (document.get() == nullptr) {
        throw std::runtime_error("top-level JSON parse failed at line " +
                                 std::to_string(error.line) + ": " + error.text);
    }
    if (!json_is_object(document.get())) {
        throw std::runtime_error("top-level JSON must be an object");
    }
    const std::string schema_name =
        require_string(document.get(), "schema_name", "document");
    const std::string schema_version =
        require_string(document.get(), "schema_version", "document");
    if (schema_name != kInputSchemaName || schema_version != kInputSchemaVersion) {
        throw std::runtime_error("unsupported input schema " + schema_name + " " +
                                 schema_version);
    }
    const std::string sample = require_string(document.get(), "sample", "document");
    json_t* groups = json_object_get(document.get(), "groups");
    if (!json_is_array(groups)) {
        throw std::runtime_error("document.groups must be an array");
    }

    const std::string started_at = utc_now();
    const auto start = std::chrono::steady_clock::now();
    const std::filesystem::path output_temp = temporary_sibling(arguments.output);
    const std::filesystem::path receipt_temp = temporary_sibling(arguments.receipt);
    if (std::filesystem::exists(output_temp) || std::filesystem::exists(receipt_temp)) {
        throw std::runtime_error("atomic temporary target already exists");
    }

    RunCounters counters;
    try {
        std::ofstream output(output_temp, std::ios::binary | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("cannot open output temporary file: " +
                                     output_temp.string());
        }
        for (std::size_t index = 0; index < json_array_size(groups); ++index) {
            UnitOutcome outcome =
                analyze_unit(json_array_get(groups, index), index, sample, arguments);
            write_json_line(output, outcome.record.get());
            ++counters.total_units;
            counters.unit_status[outcome.unit_status]++;
            counters.objective_status[outcome.objective_status]++;
            counters.family_status[outcome.family_status]++;
            counters.read_af_status[outcome.read_af_status]++;
            if (outcome.unit_status == "malformed_unit") {
                ++counters.malformed_units;
            }
            if (outcome.mutation_bearing) {
                ++counters.mutation_bearing_units;
                if (outcome.family_complete) {
                    ++counters.family_complete_mutation_bearing_units;
                }
            }
            if (outcome.objective_certified) {
                ++counters.objective_certified_units;
            }
            if (outcome.ranking_complete) {
                ++counters.ranking_complete_count;
            }
        }
        output.flush();
        if (!output) {
            throw std::runtime_error("failed to flush output temporary file");
        }
        output.close();

        const FileIdentity input_identity = identify_file(arguments.input);
        FileIdentity output_identity = identify_file(output_temp);
        output_identity.path =
            std::filesystem::absolute(arguments.output).lexically_normal().string();
        const auto finish = std::chrono::steady_clock::now();
        const std::string finished_at = utc_now();
        const double elapsed =
            std::chrono::duration<double>(finish - start).count();
        JsonOwner receipt(build_receipt(arguments, input_identity, output_identity, counters,
                                        started_at, finished_at, elapsed));
        write_json_file(receipt_temp, receipt.get());

        std::filesystem::rename(output_temp, arguments.output);
        std::filesystem::rename(receipt_temp, arguments.receipt);
        std::cerr << "sample=" << sample << " units=" << counters.total_units
                  << " objective_certified=" << counters.objective_certified_units
                  << " mutation_family_complete="
                  << counters.family_complete_mutation_bearing_units << '/'
                  << counters.mutation_bearing_units
                  << " ranked=" << counters.ranking_complete_count
                  << " elapsed_seconds=" << std::fixed << std::setprecision(6)
                  << elapsed << '\n';
    } catch (...) {
        std::cerr << "atomic publication aborted; diagnostic temporary files, if any, were retained: "
                  << output_temp << " " << receipt_temp << '\n';
        throw;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(parse_arguments(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
