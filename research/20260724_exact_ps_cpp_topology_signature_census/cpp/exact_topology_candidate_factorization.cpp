// SPDX-License-Identifier: GPL-3.0-only
//
// Exact candidate-space factorization companion for the frozen 2026-07-24
// exact-PS topology run. The frozen implementation is included directly so
// parsing, structural search, and exact read-AF scoring cannot drift.

#define main exact_ps_topology_af_frozen_main
#include "../../20260724_exact_ps_cpp_topology_af_all_samples/cpp/exact_ps_topology_af.cpp"
#undef main

#include <cerrno>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <functional>
#include <linux/fs.h>
#include <limits>
#include <sys/syscall.h>
#include <unistd.h>

namespace {

constexpr const char* kFactorizationSchemaName =
    "intersubmod.exact_ps_topology_candidate_factorization.unit";
constexpr const char* kFactorizationSchemaVersion = "1.0.0";
constexpr const char* kCanonicalSchemaName =
    "intersubmod.exact_ps_cpp_topology_af.unit";
constexpr const char* kCanonicalSchemaVersion = "1.0.0";
constexpr const char* kCensusSchemaName =
    "intersubmod.exact_ps_cpp_topology_signature_census.unit";
constexpr const char* kCensusSchemaVersion = "1.0.0";

struct FactorizationArguments {
    std::filesystem::path input;
    std::filesystem::path canonical;
    std::filesystem::path census;
    std::filesystem::path output;
    std::size_t maximum_family_size = 100000;
    std::uint64_t maximum_search_nodes = 1000;
};

FactorizationArguments parse_factorization_arguments(int argc, char** argv) {
    FactorizationArguments arguments;
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
        } else if (option == "--canonical") {
            arguments.canonical = next();
        } else if (option == "--census") {
            arguments.census = next();
        } else if (option == "--output") {
            arguments.output = next();
        } else if (option == "--max-family-size") {
            const std::uint64_t parsed = parse_limit(next(), option);
            if (parsed > std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error(option + " exceeds size_t range");
            }
            arguments.maximum_family_size = static_cast<std::size_t>(parsed);
        } else if (option == "--max-search-nodes") {
            arguments.maximum_search_nodes = parse_limit(next(), option);
        } else {
            throw std::runtime_error("unknown option: " + option);
        }
    }
    if (arguments.input.empty() || arguments.canonical.empty() ||
        arguments.census.empty() || arguments.output.empty()) {
        throw std::runtime_error(
            "--input, --canonical, --census and --output are required");
    }
    return arguments;
}

void append_new(json_t* array, json_t* value, const std::string& label) {
    if (value == nullptr) {
        throw std::runtime_error("cannot append null JSON allocation: " +
                                 label);
    }
    // Jansson's *_append_new API steals the value reference on all paths.
    if (json_array_append_new(array, value) != 0) {
        throw std::runtime_error("failed to append " + label);
    }
}

json_t* require_array(json_t* object, const char* key,
                      const std::string& label) {
    json_t* value = json_object_get(object, key);
    if (!json_is_array(value)) {
        throw std::runtime_error(label + "." + key + " must be an array");
    }
    return value;
}

json_t* require_object(json_t* object, const char* key,
                       const std::string& label) {
    json_t* value = json_object_get(object, key);
    if (!json_is_object(value)) {
        throw std::runtime_error(label + "." + key + " must be an object");
    }
    return value;
}

bool require_boolean(json_t* object, const char* key,
                     const std::string& label) {
    json_t* value = json_object_get(object, key);
    if (!json_is_boolean(value)) {
        throw std::runtime_error(label + "." + key + " must be a boolean");
    }
    return json_is_true(value);
}

JsonOwner load_factorization_json_line(const std::string& line,
                                       std::size_t line_number,
                                       const std::string& label) {
    json_error_t error{};
    JsonOwner row(json_loads(line.c_str(), JSON_REJECT_DUPLICATES, &error));
    if (!json_is_object(row.get())) {
        throw std::runtime_error(
            label + " JSONL line " + std::to_string(line_number) +
            " is not a valid object: " + error.text);
    }
    return row;
}

json_t* integer_array(const std::vector<HypercubeVertex>& values) {
    json_t* result = json_array();
    if (result == nullptr) {
        throw std::runtime_error("failed to allocate integer array");
    }
    for (HypercubeVertex value : values) {
        append_new(result, json_integer(static_cast<json_int_t>(value)),
                   "integer array value");
    }
    return result;
}

bool json_integer_array_equals(json_t* value,
                               const std::vector<std::int64_t>& expected) {
    if (!json_is_array(value) || json_array_size(value) != expected.size()) {
        return false;
    }
    for (std::size_t index = 0; index < expected.size(); ++index) {
        json_t* item = json_array_get(value, index);
        if (!json_is_integer(item) ||
            json_integer_value(item) != expected[index]) {
            return false;
        }
    }
    return true;
}

std::vector<std::int64_t> active_indices_as_integers(
    const ParsedUnit& unit) {
    std::vector<std::int64_t> result;
    result.reserve(unit.active_original_indices.size());
    for (std::size_t value : unit.active_original_indices) {
        result.push_back(static_cast<std::int64_t>(value));
    }
    return result;
}

std::vector<std::int64_t> active_positions_as_integers(
    const ParsedUnit& unit) {
    std::vector<std::int64_t> result;
    result.reserve(unit.active_original_indices.size());
    for (std::size_t index : unit.active_original_indices) {
        result.push_back(unit.positions[index]);
    }
    return result;
}

void verify_schema(json_t* row, const std::string& expected_name,
                   const std::string& expected_version,
                   const std::string& label) {
    if (require_string(row, "schema_name", label) != expected_name ||
        require_string(row, "schema_version", label) != expected_version) {
        throw std::runtime_error(label + " has an unsupported schema");
    }
}

void verify_identity(json_t* row, const ParsedUnit& unit,
                     std::size_t group_index, const std::string& label) {
    if (require_integer(row, "group_index", label) !=
            static_cast<json_int_t>(group_index) ||
        require_string(row, "sample", label) != unit.sample ||
        require_string(row, "region_id", label) != unit.region_id ||
        require_string(row, "unit_id", label) != unit.unit_id ||
        require_string(row, "block_id", label) != unit.block_id ||
        require_string(row, "chrom", label) != unit.chromosome ||
        require_string(row, "hp_family", label) != unit.hp_family ||
        require_string(row, "phase_set", label) != unit.phase_set) {
        throw std::runtime_error(label + " identity disagrees with input group");
    }
}

std::string factorization_coarse_class(bool has_sister, bool has_direct) {
    if (!has_sister && !has_direct) {
        return "Single-only";
    }
    if (has_sister && !has_direct) {
        return "Sister-only";
    }
    if (!has_sister && has_direct) {
        return "Direct-only";
    }
    return "Sister+direct";
}

struct FactorizedTreeSignature {
    std::string signature;
    std::string coarse;
};

FactorizedTreeSignature factorization_signature_for_tree(
    const std::vector<HypercubeVertex>& vertices,
    const std::vector<std::pair<HypercubeVertex, HypercubeVertex>>& edges) {
    std::map<HypercubeVertex, std::vector<HypercubeVertex>> children;
    bool has_sister = false;
    bool has_direct = false;
    for (HypercubeVertex vertex : vertices) {
        has_direct = has_direct || vertex_popcount(vertex) >= 2;
    }
    for (const auto& [parent, child] : edges) {
        children[parent].push_back(child);
    }
    for (auto& entry : children) {
        std::vector<HypercubeVertex>& child_values = entry.second;
        std::sort(child_values.begin(), child_values.end());
        has_sister = has_sister || child_values.size() >= 2;
    }
    return {unlabeled_shape_at(0, children),
            factorization_coarse_class(has_sister, has_direct)};
}

struct ChildFactorization {
    HypercubeVertex child = 0;
    std::vector<HypercubeVertex> legal_parents;
    std::vector<HypercubeVertex> best_parents;
    std::vector<std::pair<HypercubeVertex, Rational>> edge_scores;
    Rational best_score;
};

struct VertexSetFactorization {
    std::vector<HypercubeVertex> vertices;
    std::vector<ChildFactorization> children;
    Rational best_score;
    cpp_int tree_count = 1;
    cpp_int best_tree_count = 1;
    bool is_global_best = false;
};

VertexSetFactorization factorize_vertex_set(
    const ExactStructuralCandidate& candidate, std::size_t bit_count,
    const std::vector<Rational>& read_af) {
    const ExactParentMappingSummary summary =
        longlineage::solver::summarize_exact_parent_mappings(
            bit_count, candidate.vertices);
    if (!summary.valid || !candidate.parent_mapping.valid ||
        summary.tree_count != candidate.parent_mapping.tree_count ||
        summary.legal_parent_count !=
            candidate.parent_mapping.legal_parent_count) {
        throw std::runtime_error(
            "pinned parent factorization disagrees with structural candidate");
    }

    VertexSetFactorization result;
    result.vertices = candidate.vertices;
    std::sort(result.vertices.begin(), result.vertices.end());
    result.tree_count =
        parse_nonnegative_integer(summary.tree_count, "summary.tree_count");
    cpp_int reproduced_tree_count = 1;

    for (const ExactLegalParentChoices& choices : summary.legal_parents) {
        if (choices.parents.empty()) {
            throw std::runtime_error(
                "valid parent factorization contains no legal parent");
        }
        ChildFactorization child;
        child.child = choices.vertex;
        child.legal_parents = choices.parents;
        std::sort(child.legal_parents.begin(), child.legal_parents.end());
        child.legal_parents.erase(
            std::unique(child.legal_parents.begin(),
                        child.legal_parents.end()),
            child.legal_parents.end());
        if (child.legal_parents.size() != choices.parents.size()) {
            throw std::runtime_error(
                "parent factorization contains duplicate legal parents");
        }
        std::optional<Rational> best;
        for (HypercubeVertex parent : child.legal_parents) {
            const Rational score = edge_score(parent, child.child, read_af);
            child.edge_scores.emplace_back(parent, score);
            if (!best.has_value() || score > *best) {
                best = score;
                child.best_parents = {parent};
            } else if (score == *best) {
                child.best_parents.push_back(parent);
            }
        }
        if (!best.has_value() || child.best_parents.empty()) {
            throw std::runtime_error(
                "valid parent factorization has no best parent");
        }
        child.best_score = *best;
        result.best_score = result.best_score + child.best_score;
        result.best_tree_count *=
            static_cast<unsigned long long>(child.best_parents.size());
        reproduced_tree_count *=
            static_cast<unsigned long long>(child.legal_parents.size());
        result.children.push_back(std::move(child));
    }
    std::sort(result.children.begin(), result.children.end(),
              [](const ChildFactorization& left,
                 const ChildFactorization& right) {
                  return left.child < right.child;
              });
    if (reproduced_tree_count != result.tree_count) {
        throw std::runtime_error(
            "legal-parent product does not reproduce candidate tree_count");
    }
    return result;
}

using Edge = std::pair<HypercubeVertex, HypercubeVertex>;
using EdgeIncidence = std::map<Edge, cpp_int>;

void publish_no_replace(const std::filesystem::path& pending,
                        const std::filesystem::path& target) {
    const long result = syscall(
        SYS_renameat2, AT_FDCWD, pending.c_str(), AT_FDCWD,
        target.c_str(), RENAME_NOREPLACE);
    if (result != 0) {
        throw std::runtime_error(
            "atomic no-replace publish failed for " + target.string() +
            ": " + std::strerror(errno));
    }
}

std::filesystem::path create_exclusive_pending(
    const std::filesystem::path& target) {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch()
                           .count();
    const std::filesystem::path pending =
        target.parent_path() /
        (target.filename().string() + ".pending." +
         std::to_string(static_cast<long long>(getpid())) + "." +
         std::to_string(stamp));
    const int descriptor = open(
        pending.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0664);
    if (descriptor < 0) {
        throw std::runtime_error(
            "cannot exclusively create temporary output " +
            pending.string() + ": " + std::strerror(errno));
    }
    if (close(descriptor) != 0) {
        throw std::runtime_error(
            "cannot close newly created temporary output " +
            pending.string() + ": " + std::strerror(errno));
    }
    return pending;
}

void sync_file(const std::filesystem::path& path) {
    const int descriptor = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (descriptor < 0) {
        throw std::runtime_error(
            "cannot open output for fsync " + path.string() + ": " +
            std::strerror(errno));
    }
    if (fsync(descriptor) != 0) {
        const std::string message =
            "cannot fsync output " + path.string() + ": " +
            std::strerror(errno);
        close(descriptor);
        throw std::runtime_error(message);
    }
    if (close(descriptor) != 0) {
        throw std::runtime_error(
            "cannot close fsynced output " + path.string() + ": " +
            std::strerror(errno));
    }
}

class PendingOutputGuard {
  public:
    explicit PendingOutputGuard(std::filesystem::path path)
        : path_(std::move(path)) {}
    ~PendingOutputGuard() {
        if (!active_) {
            return;
        }
        const auto stamp = std::chrono::steady_clock::now()
                               .time_since_epoch()
                               .count();
        const std::filesystem::path failed =
            path_.parent_path() /
            (path_.filename().string() + ".failed." +
             std::to_string(static_cast<long long>(getpid())) + "." +
             std::to_string(stamp));
        try {
            publish_no_replace(path_, failed);
        } catch (...) {
            // Preserve the original pending artifact if no-replace archival
            // cannot be completed; destructors must not throw.
        }
    }
    PendingOutputGuard(const PendingOutputGuard&) = delete;
    PendingOutputGuard& operator=(const PendingOutputGuard&) = delete;
    void release() { active_ = false; }

  private:
    std::filesystem::path path_;
    bool active_ = true;
};

void add_edge_incidence(const VertexSetFactorization& factorization,
                        bool best_only, EdgeIncidence& incidence) {
    const cpp_int product =
        best_only ? factorization.best_tree_count : factorization.tree_count;
    for (const ChildFactorization& child : factorization.children) {
        const std::vector<HypercubeVertex>& parents =
            best_only ? child.best_parents : child.legal_parents;
        if (parents.empty()) {
            throw std::runtime_error(
                "cannot calculate incidence with no parent choices");
        }
        const cpp_int divisor =
            static_cast<unsigned long long>(parents.size());
        if (product % divisor != 0) {
            throw std::runtime_error(
                "tree product is not divisible by parent-choice count");
        }
        const cpp_int count = product / divisor;
        for (HypercubeVertex parent : parents) {
            incidence[{parent, child.child}] += count;
        }
    }
}

struct SignatureAggregate {
    std::string coarse;
    cpp_int tree_count = 0;
    std::vector<HypercubeVertex> exemplar_vertices;
    std::vector<Edge> exemplar_edges;
};

using SignatureAggregates = std::map<std::string, SignatureAggregate>;

void enumerate_global_best_trees(
    const VertexSetFactorization& factorization, std::size_t index,
    std::vector<Edge>& edges, SignatureAggregates& signatures,
    std::map<std::string, cpp_int>& coarse_counts,
    cpp_int& enumerated_count) {
    if (index == factorization.children.size()) {
        const FactorizedTreeSignature signature =
            factorization_signature_for_tree(factorization.vertices, edges);
        auto [iterator, inserted] = signatures.emplace(
            signature.signature, SignatureAggregate{});
        SignatureAggregate& aggregate = iterator->second;
        if (inserted) {
            aggregate.coarse = signature.coarse;
            aggregate.exemplar_vertices = factorization.vertices;
            aggregate.exemplar_edges = edges;
            std::sort(aggregate.exemplar_edges.begin(),
                      aggregate.exemplar_edges.end());
        } else if (aggregate.coarse != signature.coarse) {
            throw std::runtime_error(
                "one rooted-unlabeled signature mapped to multiple coarse "
                "classes");
        }
        aggregate.tree_count += 1;
        coarse_counts[signature.coarse] += 1;
        enumerated_count += 1;
        return;
    }
    const ChildFactorization& child = factorization.children[index];
    for (HypercubeVertex parent : child.best_parents) {
        edges.emplace_back(parent, child.child);
        enumerate_global_best_trees(factorization, index + 1, edges,
                                    signatures, coarse_counts,
                                    enumerated_count);
        edges.pop_back();
    }
}

json_t* compact_edge_array(const std::vector<Edge>& edges) {
    json_t* result = json_array();
    if (result == nullptr) {
        throw std::runtime_error("failed to allocate compact edge array");
    }
    for (const auto& [parent, child] : edges) {
        json_t* edge = json_array();
        append_new(edge, json_integer(static_cast<json_int_t>(parent)),
                   "edge parent");
        append_new(edge, json_integer(static_cast<json_int_t>(child)),
                   "edge child");
        append_new(result, edge, "compact edge");
    }
    return result;
}

json_t* edge_incidence_json(const EdgeIncidence& incidence) {
    json_t* result = json_array();
    if (result == nullptr) {
        throw std::runtime_error("failed to allocate edge incidence array");
    }
    for (const auto& [edge, count] : incidence) {
        json_t* row = json_array();
        append_new(row,
                   json_integer(static_cast<json_int_t>(edge.first)),
                   "incidence parent");
        append_new(row,
                   json_integer(static_cast<json_int_t>(edge.second)),
                   "incidence child");
        append_new(row, json_string(decimal(count).c_str()),
                   "incidence count");
        append_new(result, row, "edge incidence row");
    }
    return result;
}

json_t* parent_factorization_json(
    const VertexSetFactorization& factorization) {
    json_t* result = json_array();
    if (result == nullptr) {
        throw std::runtime_error(
            "failed to allocate parent factorization array");
    }
    for (const ChildFactorization& child : factorization.children) {
        json_t* row = json_array();
        append_new(row, json_integer(static_cast<json_int_t>(child.child)),
                   "factorization child");
        append_new(row, integer_array(child.legal_parents),
                   "factorization legal parents");
        append_new(row, integer_array(child.best_parents),
                   "factorization best parents");
        append_new(result, row, "parent factorization row");
    }
    return result;
}

json_t* legal_edge_scores_json(
    const VertexSetFactorization& factorization) {
    std::vector<std::tuple<HypercubeVertex, HypercubeVertex, Rational>>
        scores;
    for (const ChildFactorization& child : factorization.children) {
        for (const auto& [parent, score] : child.edge_scores) {
            scores.emplace_back(parent, child.child, score);
        }
    }
    std::sort(scores.begin(), scores.end(),
              [](const auto& left, const auto& right) {
                  if (std::get<0>(left) != std::get<0>(right)) {
                      return std::get<0>(left) < std::get<0>(right);
                  }
                  return std::get<1>(left) < std::get<1>(right);
              });
    json_t* result = json_array();
    if (result == nullptr) {
        throw std::runtime_error("failed to allocate edge score array");
    }
    for (const auto& [parent, child, score] : scores) {
        json_t* row = json_array();
        append_new(row, json_integer(static_cast<json_int_t>(parent)),
                   "edge score parent");
        append_new(row, json_integer(static_cast<json_int_t>(child)),
                   "edge score child");
        append_new(row, json_string(score.str().c_str()),
                   "edge score fraction");
        append_new(result, row, "legal edge score row");
    }
    return result;
}

json_t* minimum_vertex_sets_json(
    const std::vector<VertexSetFactorization>& factorizations) {
    json_t* result = json_array();
    if (result == nullptr) {
        throw std::runtime_error(
            "failed to allocate minimum vertex set array");
    }
    for (const VertexSetFactorization& factorization : factorizations) {
        json_t* row = json_object();
        set_new(row, "vertices", integer_array(factorization.vertices));
        set_string(row, "best_score_fraction",
                   factorization.best_score.str());
        set_string(row, "total_tree_count",
                   decimal(factorization.tree_count));
        set_string(row, "best_tree_count",
                   decimal(factorization.best_tree_count));
        set_new(row, "is_global_best",
                json_boolean(factorization.is_global_best));
        set_new(row, "parent_factorization",
                parent_factorization_json(factorization));
        set_new(row, "legal_edge_score_fractions",
                legal_edge_scores_json(factorization));
        append_new(result, row, "minimum vertex set");
    }
    return result;
}

json_t* signature_aggregates_json(
    const SignatureAggregates& signatures) {
    json_t* result = json_array();
    if (result == nullptr) {
        throw std::runtime_error("failed to allocate signature array");
    }
    for (const auto& [signature, aggregate] : signatures) {
        json_t* row = json_object();
        set_string(row, "shape_signature", signature);
        set_string(row, "shape_sha256", sha256_text(signature));
        set_string(row, "coarse_class", aggregate.coarse);
        set_string(row, "tree_count", decimal(aggregate.tree_count));
        set_new(row, "exemplar_vertices",
                integer_array(aggregate.exemplar_vertices));
        set_new(row, "exemplar_edges",
                compact_edge_array(aggregate.exemplar_edges));
        append_new(result, row, "global-best signature");
    }
    return result;
}

json_t* string_integer_map_json(
    const std::map<std::string, cpp_int>& values) {
    json_t* result = json_object();
    if (result == nullptr) {
        throw std::runtime_error("failed to allocate string-count object");
    }
    for (const auto& [key, value] : values) {
        set_string(result, key.c_str(), decimal(value));
    }
    return result;
}

struct ExpectedSignature {
    std::string shape_sha256;
    std::string coarse;
    cpp_int tree_count = 0;
};

void verify_census_signatures(
    json_t* census, const SignatureAggregates& signatures,
    const std::map<std::string, cpp_int>& coarse_counts) {
    const json_int_t expected_signature_count =
        require_integer(census, "topology_signature_count", "census");
    if (expected_signature_count < 0 ||
        static_cast<std::size_t>(expected_signature_count) !=
            signatures.size()) {
        throw std::runtime_error(
            "census topology_signature_count did not reproduce");
    }
    std::map<std::string, ExpectedSignature> expected;
    json_t* rows =
        require_array(census, "topology_signatures", "census");
    std::size_t index = 0;
    json_t* row = nullptr;
    json_array_foreach(rows, index, row) {
        if (!json_is_object(row)) {
            throw std::runtime_error(
                "census.topology_signatures contains a non-object");
        }
        const std::string signature =
            require_string(row, "shape_signature", "census signature");
        ExpectedSignature value;
        value.shape_sha256 =
            require_string(row, "shape_sha256", "census signature");
        value.coarse =
            require_string(row, "coarse_class", "census signature");
        value.tree_count = parse_nonnegative_integer(
            require_string(row, "tree_count", "census signature"),
            "census signature.tree_count");
        if (!expected.emplace(signature, std::move(value)).second) {
            throw std::runtime_error(
                "census contains a duplicate topology signature");
        }
    }
    if (expected.size() != signatures.size()) {
        throw std::runtime_error(
            "census topology signature set size did not reproduce");
    }
    for (const auto& [signature, aggregate] : signatures) {
        const auto iterator = expected.find(signature);
        if (iterator == expected.end() ||
            iterator->second.shape_sha256 != sha256_text(signature) ||
            iterator->second.coarse != aggregate.coarse ||
            iterator->second.tree_count != aggregate.tree_count) {
            throw std::runtime_error(
                "census topology signature counts did not reproduce");
        }
    }

    json_t* expected_coarse =
        require_object(census, "coarse_class_tree_counts", "census");
    if (json_object_size(expected_coarse) != coarse_counts.size()) {
        throw std::runtime_error(
            "census coarse-class set size did not reproduce");
    }
    for (const auto& [coarse, count] : coarse_counts) {
        json_t* value = json_object_get(expected_coarse, coarse.c_str());
        if (!json_is_string(value) ||
            parse_nonnegative_integer(
                json_string_value(value),
                "census.coarse_class_tree_counts") != count) {
            throw std::runtime_error(
                "census coarse-class tree counts did not reproduce");
        }
    }
}

JsonOwner analyze_ranked_factorization(
    json_t* group, std::size_t group_index, const std::string& sample,
    json_t* canonical, json_t* census,
    const FactorizationArguments& arguments) {
    ParsedUnit unit = parse_unit(group, group_index, sample);
    if (unit.coverage_state != CoverageState::kComplete ||
        unit.active_original_indices.empty()) {
        throw std::runtime_error(
            "ranked row does not have complete active read-AF input");
    }
    verify_schema(canonical, kCanonicalSchemaName,
                  kCanonicalSchemaVersion, "canonical");
    verify_schema(census, kCensusSchemaName, kCensusSchemaVersion,
                  "census");
    verify_identity(canonical, unit, group_index, "canonical");
    verify_identity(census, unit, group_index, "census");
    if (!require_boolean(census, "canonical_reproduction_pass",
                         "census")) {
        throw std::runtime_error(
            "census does not certify canonical reproduction");
    }

    const std::vector<std::int64_t> active_indices =
        active_indices_as_integers(unit);
    const std::vector<std::int64_t> active_positions =
        active_positions_as_integers(unit);
    if (!json_integer_array_equals(
            json_object_get(canonical, "active_original_indices"),
            active_indices) ||
        !json_integer_array_equals(
            json_object_get(canonical, "active_positions"),
            active_positions) ||
        require_integer(canonical, "active_bit_count", "canonical") !=
            static_cast<json_int_t>(active_indices.size()) ||
        require_integer(census, "active_bit_count", "census") !=
            static_cast<json_int_t>(active_indices.size())) {
        throw std::runtime_error(
            "active indices/positions did not reproduce canonical inputs");
    }

    ObligationBnbOptions options;
    options.maximum_complete_family_size = arguments.maximum_family_size;
    options.maximum_search_nodes = arguments.maximum_search_nodes;
    const ExactStructuralResult structural =
        longlineage::solver::solve_obligation_bnb(unit.problem, options);
    if (structural.family_state != ExactFamilyState::kFamilyComplete) {
        throw std::runtime_error(
            "ranked row did not reproduce a complete minimum family");
    }
    const std::string family_sha256 = sha256_text(
        canonical_family_text(structural, active_indices.size()));
    if (family_sha256 !=
        require_string(canonical, "minimum_family_sha256", "canonical")) {
        throw std::runtime_error("minimum-family SHA-256 mismatch");
    }

    std::vector<VertexSetFactorization> factorizations;
    factorizations.reserve(structural.minimum_family.size());
    for (const ExactStructuralCandidate& candidate :
         structural.minimum_family) {
        factorizations.push_back(factorize_vertex_set(
            candidate, active_indices.size(), unit.active_read_af));
    }
    std::sort(factorizations.begin(), factorizations.end(),
              [](const VertexSetFactorization& left,
                 const VertexSetFactorization& right) {
                  return left.vertices < right.vertices;
              });
    for (std::size_t index = 1; index < factorizations.size(); ++index) {
        if (factorizations[index - 1].vertices ==
            factorizations[index].vertices) {
            throw std::runtime_error(
                "minimum family contains a duplicate vertex set");
        }
    }
    if (factorizations.empty()) {
        throw std::runtime_error(
            "complete minimum family unexpectedly contains no candidates");
    }

    cpp_int total_tree_count = 0;
    std::optional<Rational> global_best;
    for (VertexSetFactorization& factorization : factorizations) {
        total_tree_count += factorization.tree_count;
        if (!global_best.has_value() ||
            factorization.best_score > *global_best) {
            global_best = factorization.best_score;
        }
    }
    if (!global_best.has_value()) {
        throw std::runtime_error("minimum family has no best score");
    }

    std::size_t best_vertex_set_count = 0;
    cpp_int best_tree_tie_count = 0;
    EdgeIncidence all_minimum_incidence;
    EdgeIncidence global_best_incidence;
    SignatureAggregates signatures;
    std::map<std::string, cpp_int> coarse_counts;
    cpp_int enumerated_best_tree_count = 0;
    std::vector<Edge> edges;
    for (VertexSetFactorization& factorization : factorizations) {
        add_edge_incidence(factorization, false,
                           all_minimum_incidence);
        factorization.is_global_best =
            factorization.best_score == *global_best;
        if (!factorization.is_global_best) {
            continue;
        }
        ++best_vertex_set_count;
        best_tree_tie_count += factorization.best_tree_count;
        add_edge_incidence(factorization, true,
                           global_best_incidence);
        edges.clear();
        edges.reserve(factorization.children.size());
        enumerate_global_best_trees(
            factorization, 0, edges, signatures, coarse_counts,
            enumerated_best_tree_count);
    }
    if (best_vertex_set_count == 0 ||
        enumerated_best_tree_count != best_tree_tie_count) {
        throw std::runtime_error(
            "global-best tree enumeration disagrees with factorization");
    }

    const json_int_t canonical_minimum_sets =
        require_integer(canonical, "minimum_vertex_set_count",
                        "canonical");
    const json_int_t census_minimum_sets =
        require_integer(census, "minimum_vertex_set_count", "census");
    const json_int_t canonical_best_sets =
        require_integer(canonical, "best_vertex_set_count", "canonical");
    const json_int_t census_best_sets =
        require_integer(census, "best_vertex_set_count", "census");
    if (canonical_minimum_sets < 0 || census_minimum_sets < 0 ||
        canonical_best_sets < 0 || census_best_sets < 0 ||
        static_cast<std::size_t>(canonical_minimum_sets) !=
            factorizations.size() ||
        static_cast<std::size_t>(census_minimum_sets) !=
            factorizations.size() ||
        static_cast<std::size_t>(canonical_best_sets) !=
            best_vertex_set_count ||
        static_cast<std::size_t>(census_best_sets) !=
            best_vertex_set_count) {
        throw std::runtime_error(
            "minimum/global-best vertex-set counts did not reproduce");
    }
    if (require_string(canonical, "total_tree_count", "canonical") !=
            decimal(total_tree_count) ||
        require_string(canonical, "best_score_fraction", "canonical") !=
            global_best->str() ||
        require_string(census, "best_score_fraction", "census") !=
            global_best->str() ||
        require_string(canonical, "best_tree_tie_count", "canonical") !=
            decimal(best_tree_tie_count) ||
        require_string(census, "best_tree_tie_count", "census") !=
            decimal(best_tree_tie_count) ||
        require_string(census, "enumerated_best_tree_count", "census") !=
            decimal(enumerated_best_tree_count) ||
        require_boolean(canonical, "best_tree_unique", "canonical") !=
            (best_tree_tie_count == 1)) {
        throw std::runtime_error(
            "tree total/best score/tie did not reproduce canonical/census");
    }

    verify_census_signatures(census, signatures, coarse_counts);
    const std::string resolution =
        best_tree_tie_count == 1
            ? "UNIQUE_TREE"
            : (signatures.size() == 1 ? "TIED_SAME_TOPOLOGY"
                                      : "TIED_CROSS_TOPOLOGY");
    if (require_string(census, "resolution_class", "census") !=
        resolution) {
        throw std::runtime_error(
            "census resolution class did not reproduce");
    }

    JsonOwner output(json_object());
    set_string(output.get(), "schema_name", kFactorizationSchemaName);
    set_string(output.get(), "schema_version",
               kFactorizationSchemaVersion);
    set_string(output.get(), "sample", unit.sample);
    set_new(output.get(), "group_index",
            json_integer(static_cast<json_int_t>(group_index)));
    set_string(output.get(), "region_id", unit.region_id);
    set_string(output.get(), "unit_id", unit.unit_id);
    set_string(output.get(), "block_id", unit.block_id);
    set_string(output.get(), "chrom", unit.chromosome);
    set_string(output.get(), "hp_family", unit.hp_family);
    set_string(output.get(), "phase_set", unit.phase_set);
    set_new(output.get(), "active_original_indices",
            active_indices_json(unit));
    set_new(output.get(), "active_positions",
            active_positions_json(unit));
    set_new(output.get(), "active_bit_count",
            json_integer(static_cast<json_int_t>(active_indices.size())));
    std::vector<HypercubeVertex> observed_vertices(
        unit.full_observed_vertices.begin(),
        unit.full_observed_vertices.end());
    set_new(output.get(), "observed_vertices",
            integer_array(observed_vertices));
    set_string(
        output.get(), "vertex_bit_encoding",
        "bit i corresponds to active_positions[i] in ascending active-index "
        "order; 1=ALT, 0=REF; vertex 0 is ROOT; vertices absent from "
        "observed_vertices are solver-required hidden ancestors");
    set_string(output.get(), "minimum_family_sha256", family_sha256);
    set_new(output.get(), "minimum_vertex_set_count",
            json_integer(static_cast<json_int_t>(
                factorizations.size())));
    set_string(output.get(), "total_tree_count",
               decimal(total_tree_count));
    set_string(output.get(), "best_score_fraction",
               global_best->str());
    set_new(output.get(), "best_vertex_set_count",
            json_integer(static_cast<json_int_t>(
                best_vertex_set_count)));
    set_string(output.get(), "best_tree_tie_count",
               decimal(best_tree_tie_count));
    set_string(output.get(), "resolution_class", resolution);
    set_new(output.get(), "minimum_vertex_sets",
            minimum_vertex_sets_json(factorizations));
    set_string(
        output.get(), "all_minimum_tree_edge_incidence_definition",
        "count of exact legal parent mappings containing [parent,child] "
        "across every minimum vertex set");
    set_new(output.get(), "all_minimum_tree_edge_incidence",
            edge_incidence_json(all_minimum_incidence));
    set_string(
        output.get(), "global_best_edge_incidence_definition",
        "count of AF-global-best parent mappings containing "
        "[parent,child] across every global-best vertex set");
    set_new(output.get(), "global_best_edge_incidence",
            edge_incidence_json(global_best_incidence));
    set_new(output.get(), "global_best_signature_count",
            json_integer(static_cast<json_int_t>(signatures.size())));
    set_new(output.get(), "global_best_signatures",
            signature_aggregates_json(signatures));
    set_new(output.get(), "global_best_coarse_class_tree_counts",
            string_integer_map_json(coarse_counts));
    set_new(output.get(), "canonical_reproduction_pass", json_true());
    set_new(output.get(), "census_reproduction_pass", json_true());
    return output;
}

int run_factorization(const FactorizationArguments& arguments) {
    for (const auto& path :
         {arguments.input, arguments.canonical, arguments.census}) {
        if (!std::filesystem::is_regular_file(path)) {
            throw std::runtime_error("input is not a regular file: " +
                                     path.string());
        }
    }
    if (std::filesystem::exists(arguments.output)) {
        throw std::runtime_error("output target already exists");
    }
    const std::filesystem::path output_parent =
        arguments.output.parent_path();
    if (!output_parent.empty()) {
        std::filesystem::create_directories(output_parent);
    }

    json_error_t error{};
    JsonOwner document(json_load_file(arguments.input.c_str(),
                                      JSON_REJECT_DUPLICATES, &error));
    if (!json_is_object(document.get())) {
        throw std::runtime_error("input JSON parse failed: " +
                                 std::string(error.text));
    }
    if (require_string(document.get(), "schema_name", "document") !=
            kInputSchemaName ||
        require_string(document.get(), "schema_version", "document") !=
            kInputSchemaVersion) {
        throw std::runtime_error("unsupported input schema");
    }
    const std::string sample =
        require_string(document.get(), "sample", "document");
    json_t* groups = require_array(document.get(), "groups", "document");

    std::ifstream canonical_input(arguments.canonical);
    if (!canonical_input) {
        throw std::runtime_error("cannot open canonical JSONL");
    }
    std::ifstream census_input(arguments.census);
    if (!census_input) {
        throw std::runtime_error("cannot open census JSONL");
    }

    const std::filesystem::path temporary =
        create_exclusive_pending(arguments.output);
    PendingOutputGuard pending_output(temporary);
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("cannot open temporary output");
    }

    std::uint64_t ranked_count = 0;
    cpp_int all_minimum_tree_count = 0;
    cpp_int global_best_tree_count = 0;
    std::string canonical_line;
    std::string census_line;
    for (std::size_t index = 0; index < json_array_size(groups);
         ++index) {
        if (!std::getline(canonical_input, canonical_line)) {
            throw std::runtime_error(
                "canonical JSONL ended before input groups");
        }
        JsonOwner canonical = load_factorization_json_line(
            canonical_line, index + 1, "canonical");
        if (require_integer(canonical.get(), "group_index",
                            "canonical") !=
            static_cast<json_int_t>(index)) {
            throw std::runtime_error(
                "canonical group_index order mismatch");
        }
        const std::string status =
            require_string(canonical.get(), "read_af_status",
                           "canonical");
        if (status != "ranked_complete") {
            continue;
        }
        if (!std::getline(census_input, census_line)) {
            throw std::runtime_error(
                "census JSONL ended before ranked canonical rows");
        }
        JsonOwner census = load_factorization_json_line(
            census_line, ranked_count + 1, "census");
        JsonOwner factorization = analyze_ranked_factorization(
            json_array_get(groups, index), index, sample,
            canonical.get(), census.get(), arguments);
        all_minimum_tree_count += parse_nonnegative_integer(
            require_string(factorization.get(), "total_tree_count",
                           "factorization"),
            "factorization.total_tree_count");
        global_best_tree_count += parse_nonnegative_integer(
            require_string(factorization.get(), "best_tree_tie_count",
                           "factorization"),
            "factorization.best_tree_tie_count");
        write_json_line(output, factorization.get());
        ++ranked_count;
    }
    if (std::getline(canonical_input, canonical_line)) {
        throw std::runtime_error(
            "canonical JSONL contains more rows than input groups");
    }
    if (std::getline(census_input, census_line)) {
        throw std::runtime_error(
            "census JSONL contains more rows than ranked canonical rows");
    }

    output.flush();
    if (!output) {
        throw std::runtime_error(
            "failed to flush candidate factorization output");
    }
    output.close();
    sync_file(temporary);
    publish_no_replace(temporary, arguments.output);
    pending_output.release();
    std::cerr << "sample=" << sample << " ranked=" << ranked_count
              << " all_minimum_trees="
              << decimal(all_minimum_tree_count)
              << " global_best_trees="
              << decimal(global_best_tree_count) << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_factorization(
            parse_factorization_arguments(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
