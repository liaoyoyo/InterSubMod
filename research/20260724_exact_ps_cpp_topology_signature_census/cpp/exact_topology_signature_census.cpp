// SPDX-License-Identifier: GPL-3.0-only
//
// Read-only companion analysis for the 2026-07-24 exact-PS topology run.
// It intentionally includes the frozen production-research implementation so
// parsing, exact rational AF scoring, and structural B&B semantics stay pinned.

#define main exact_ps_topology_af_frozen_main
#include "../../20260724_exact_ps_cpp_topology_af_all_samples/cpp/exact_ps_topology_af.cpp"
#undef main

#include <functional>

namespace {

constexpr const char* kCensusSchemaName =
    "intersubmod.exact_ps_cpp_topology_signature_census.unit";
constexpr const char* kCensusSchemaVersion = "1.0.0";

struct CensusArguments {
    std::filesystem::path input;
    std::filesystem::path canonical;
    std::filesystem::path output;
    std::size_t maximum_family_size = 100000;
    std::uint64_t maximum_search_nodes = 1000;
};

CensusArguments parse_census_arguments(int argc, char** argv) {
    CensusArguments arguments;
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
        } else if (option == "--output") {
            arguments.output = next();
        } else if (option == "--max-family-size") {
            arguments.maximum_family_size =
                static_cast<std::size_t>(parse_limit(next(), option));
        } else if (option == "--max-search-nodes") {
            arguments.maximum_search_nodes = parse_limit(next(), option);
        } else {
            throw std::runtime_error("unknown option: " + option);
        }
    }
    if (arguments.input.empty() || arguments.canonical.empty() ||
        arguments.output.empty()) {
        throw std::runtime_error("--input, --canonical and --output are required");
    }
    return arguments;
}

std::string coarse_class(bool has_sister, bool has_direct) {
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

struct TreeSignature {
    std::string signature;
    std::string coarse;
};

TreeSignature signature_for_tree(
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
    for (auto& [parent, child_values] : children) {
        (void)parent;
        std::sort(child_values.begin(), child_values.end());
        has_sister = has_sister || child_values.size() >= 2;
    }
    return {unlabeled_shape_at(0, children),
            coarse_class(has_sister, has_direct)};
}

struct BestParentProduct {
    Rational score;
    cpp_int tie_count = 1;
    std::vector<HypercubeVertex> children;
    std::vector<std::vector<HypercubeVertex>> best_parents;
};

BestParentProduct best_parent_product(
    const ExactStructuralCandidate& candidate, std::size_t bit_count,
    const std::vector<Rational>& read_af) {
    const ExactParentMappingSummary summary =
        longlineage::solver::summarize_exact_parent_mappings(
            bit_count, candidate.vertices);
    if (!summary.valid || !candidate.parent_mapping.valid ||
        summary.tree_count != candidate.parent_mapping.tree_count ||
        summary.legal_parent_count != candidate.parent_mapping.legal_parent_count) {
        throw std::runtime_error(
            "pinned parent factorization disagrees with structural candidate");
    }
    BestParentProduct result;
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
            throw std::runtime_error("selected vertex lacks a best legal parent");
        }
        std::sort(best_parents.begin(), best_parents.end());
        result.score = result.score + *best;
        result.tie_count *= best_parents.size();
        result.children.push_back(choices.vertex);
        result.best_parents.push_back(std::move(best_parents));
    }
    return result;
}

void enumerate_best_trees(
    const std::vector<HypercubeVertex>& vertices,
    const BestParentProduct& product, std::size_t index,
    std::vector<std::pair<HypercubeVertex, HypercubeVertex>>& edges,
    std::map<std::string, cpp_int>& signature_counts,
    std::map<std::string, cpp_int>& coarse_counts,
    std::map<std::string, std::string>& signature_coarse,
    cpp_int& enumerated_count) {
    if (index == product.children.size()) {
        const TreeSignature signature = signature_for_tree(vertices, edges);
        signature_counts[signature.signature] += 1;
        coarse_counts[signature.coarse] += 1;
        const auto [iterator, inserted] =
            signature_coarse.emplace(signature.signature, signature.coarse);
        if (!inserted && iterator->second != signature.coarse) {
            throw std::runtime_error(
                "one rooted-unlabeled signature mapped to multiple coarse classes");
        }
        enumerated_count += 1;
        return;
    }
    for (HypercubeVertex parent : product.best_parents[index]) {
        edges.emplace_back(parent, product.children[index]);
        enumerate_best_trees(vertices, product, index + 1, edges,
                             signature_counts, coarse_counts, signature_coarse,
                             enumerated_count);
        edges.pop_back();
    }
}

JsonOwner load_json_line(const std::string& line, std::size_t line_number) {
    json_error_t error{};
    JsonOwner row(json_loads(line.c_str(), JSON_REJECT_DUPLICATES, &error));
    if (!json_is_object(row.get())) {
        throw std::runtime_error("canonical JSONL line " +
                                 std::to_string(line_number) +
                                 " is not a valid object: " + error.text);
    }
    return row;
}

json_t* signature_count_json(
    const std::map<std::string, cpp_int>& signature_counts,
    const std::map<std::string, std::string>& signature_coarse) {
    json_t* result = json_array();
    for (const auto& [signature, count] : signature_counts) {
        json_t* row = json_object();
        set_string(row, "shape_signature", signature);
        set_string(row, "shape_sha256", sha256_text(signature));
        set_string(row, "coarse_class", signature_coarse.at(signature));
        set_string(row, "tree_count", decimal(count));
        if (json_array_append_new(result, row) != 0) {
            json_decref(result);
            throw std::runtime_error("failed to append topology signature");
        }
    }
    return result;
}

json_t* string_integer_map_json(
    const std::map<std::string, cpp_int>& values) {
    json_t* result = json_object();
    for (const auto& [key, value] : values) {
        set_string(result, key.c_str(), decimal(value));
    }
    return result;
}

JsonOwner analyze_ranked_census(
    json_t* group, std::size_t group_index, const std::string& sample,
    json_t* canonical, const CensusArguments& arguments) {
    ParsedUnit unit = parse_unit(group, group_index, sample);
    if (unit.coverage_state != CoverageState::kComplete ||
        unit.active_original_indices.empty()) {
        throw std::runtime_error(
            "canonical ranked row does not have complete active read-AF input");
    }

    ObligationBnbOptions options;
    options.maximum_complete_family_size = arguments.maximum_family_size;
    options.maximum_search_nodes = arguments.maximum_search_nodes;
    const ExactStructuralResult structural =
        longlineage::solver::solve_obligation_bnb(unit.problem, options);
    if (structural.family_state != ExactFamilyState::kFamilyComplete) {
        throw std::runtime_error(
            "canonical ranked row did not reproduce a complete minimum family");
    }
    const std::string reproduced_family_sha256 =
        sha256_text(canonical_family_text(
            structural, unit.active_original_indices.size()));
    const std::string canonical_family_sha256 =
        require_string(canonical, "minimum_family_sha256", "canonical");
    if (reproduced_family_sha256 != canonical_family_sha256) {
        throw std::runtime_error("minimum-family SHA-256 mismatch");
    }

    std::optional<Rational> global_best;
    std::vector<std::pair<const ExactStructuralCandidate*, BestParentProduct>>
        best_products;
    for (const ExactStructuralCandidate& candidate :
         structural.minimum_family) {
        BestParentProduct product =
            best_parent_product(candidate,
                                unit.active_original_indices.size(),
                                unit.active_read_af);
        if (!global_best.has_value() || product.score > *global_best) {
            global_best = product.score;
            best_products.clear();
            best_products.emplace_back(&candidate, std::move(product));
        } else if (product.score == *global_best) {
            best_products.emplace_back(&candidate, std::move(product));
        }
    }
    if (!global_best.has_value() || best_products.empty()) {
        throw std::runtime_error("complete family produced no best products");
    }

    cpp_int factorized_count = 0;
    cpp_int enumerated_count = 0;
    std::map<std::string, cpp_int> signature_counts;
    std::map<std::string, cpp_int> coarse_counts;
    std::map<std::string, std::string> signature_coarse;
    std::vector<std::pair<HypercubeVertex, HypercubeVertex>> edges;
    for (const auto& [candidate, product] : best_products) {
        factorized_count += product.tie_count;
        edges.clear();
        edges.reserve(product.children.size());
        enumerate_best_trees(candidate->vertices, product, 0, edges,
                             signature_counts, coarse_counts,
                             signature_coarse, enumerated_count);
    }

    const cpp_int canonical_tie_count = parse_nonnegative_integer(
        require_string(canonical, "best_tree_tie_count", "canonical"),
        "canonical.best_tree_tie_count");
    const std::string canonical_best_score =
        require_string(canonical, "best_score_fraction", "canonical");
    const json_int_t canonical_best_vertex_sets =
        require_integer(canonical, "best_vertex_set_count", "canonical");
    if (global_best->str() != canonical_best_score ||
        factorized_count != canonical_tie_count ||
        enumerated_count != canonical_tie_count ||
        best_products.size() !=
            static_cast<std::size_t>(canonical_best_vertex_sets)) {
        throw std::runtime_error(
            "global-best score/count did not reproduce canonical output");
    }

    json_t* representative_morphology =
        json_object_get(canonical, "representative_best_morphology");
    if (!json_is_object(representative_morphology)) {
        throw std::runtime_error(
            "canonical representative_best_morphology is missing");
    }
    const std::string representative_signature = require_string(
        representative_morphology, "shape_signature",
        "canonical.representative_best_morphology");
    if (signature_counts.count(representative_signature) == 0) {
        throw std::runtime_error(
            "canonical representative shape is absent from exact census");
    }

    const std::string resolution =
        canonical_tie_count == 1
            ? "UNIQUE_TREE"
            : (signature_counts.size() == 1
                   ? "TIED_SAME_TOPOLOGY"
                   : "TIED_CROSS_TOPOLOGY");
    JsonOwner output(json_object());
    set_string(output.get(), "schema_name", kCensusSchemaName);
    set_string(output.get(), "schema_version", kCensusSchemaVersion);
    set_string(output.get(), "sample", sample);
    set_new(output.get(), "group_index",
            json_integer(static_cast<json_int_t>(group_index)));
    for (const char* key :
         {"region_id", "unit_id", "block_id", "chrom", "hp_family",
          "phase_set"}) {
        json_t* value = json_object_get(canonical, key);
        set_new(output.get(), key,
                json_is_string(value) ? json_incref(value) : json_null());
    }
    set_new(output.get(), "active_bit_count",
            json_integer(static_cast<json_int_t>(
                unit.active_original_indices.size())));
    set_new(output.get(), "minimum_vertex_set_count",
            json_integer(static_cast<json_int_t>(
                structural.minimum_family.size())));
    set_new(output.get(), "best_vertex_set_count",
            json_integer(static_cast<json_int_t>(best_products.size())));
    set_string(output.get(), "best_score_fraction", global_best->str());
    set_string(output.get(), "best_tree_tie_count",
               decimal(canonical_tie_count));
    set_string(output.get(), "enumerated_best_tree_count",
               decimal(enumerated_count));
    set_new(output.get(), "topology_signature_count",
            json_integer(static_cast<json_int_t>(signature_counts.size())));
    set_new(output.get(), "coarse_class_count",
            json_integer(static_cast<json_int_t>(coarse_counts.size())));
    set_new(output.get(), "topology_signatures",
            signature_count_json(signature_counts, signature_coarse));
    set_new(output.get(), "coarse_class_tree_counts",
            string_integer_map_json(coarse_counts));
    set_string(output.get(), "resolution_class", resolution);
    set_string(output.get(), "topology_signature_definition",
               "root-preserving unlabeled rooted-tree canonical "
               "parenthesis signature; sibling signatures sorted");
    set_string(output.get(), "coarse_class_definition",
               "Single-only/Sister-only/Direct-only/Sister+direct from "
               "branching and root-to-node depth>=2");
    set_new(output.get(), "canonical_reproduction_pass", json_true());
    return output;
}

int run_census(const CensusArguments& arguments) {
    for (const auto& path :
         {arguments.input, arguments.canonical}) {
        if (!std::filesystem::is_regular_file(path)) {
            throw std::runtime_error("input is not a regular file: " +
                                     path.string());
        }
    }
    if (std::filesystem::exists(arguments.output)) {
        throw std::runtime_error("output target already exists");
    }
    std::filesystem::create_directories(arguments.output.parent_path());

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
    json_t* groups = json_object_get(document.get(), "groups");
    if (!json_is_array(groups)) {
        throw std::runtime_error("document.groups must be an array");
    }

    std::ifstream canonical_input(arguments.canonical);
    if (!canonical_input) {
        throw std::runtime_error("cannot open canonical JSONL");
    }
    const std::filesystem::path temporary =
        temporary_sibling(arguments.output);
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("cannot open temporary output");
    }

    std::uint64_t ranked_count = 0;
    std::uint64_t unique_count = 0;
    std::uint64_t same_topology_count = 0;
    std::uint64_t cross_topology_count = 0;
    std::string line;
    for (std::size_t index = 0; index < json_array_size(groups); ++index) {
        if (!std::getline(canonical_input, line)) {
            throw std::runtime_error(
                "canonical JSONL ended before input groups");
        }
        JsonOwner canonical = load_json_line(line, index + 1);
        if (require_integer(canonical.get(), "group_index", "canonical") !=
            static_cast<json_int_t>(index)) {
            throw std::runtime_error("canonical group_index order mismatch");
        }
        const std::string status =
            require_string(canonical.get(), "read_af_status", "canonical");
        if (status != "ranked_complete") {
            continue;
        }
        JsonOwner census = analyze_ranked_census(
            json_array_get(groups, index), index, sample, canonical.get(),
            arguments);
        const std::string resolution =
            require_string(census.get(), "resolution_class", "census");
        ++ranked_count;
        unique_count += resolution == "UNIQUE_TREE";
        same_topology_count += resolution == "TIED_SAME_TOPOLOGY";
        cross_topology_count += resolution == "TIED_CROSS_TOPOLOGY";
        write_json_line(output, census.get());
    }
    if (std::getline(canonical_input, line)) {
        throw std::runtime_error(
            "canonical JSONL contains more rows than input groups");
    }
    output.flush();
    if (!output) {
        throw std::runtime_error("failed to flush census output");
    }
    output.close();
    std::filesystem::rename(temporary, arguments.output);
    std::cerr << "sample=" << sample << " ranked=" << ranked_count
              << " unique_tree=" << unique_count
              << " tied_same_topology=" << same_topology_count
              << " tied_cross_topology=" << cross_topology_count << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_census(parse_census_arguments(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
