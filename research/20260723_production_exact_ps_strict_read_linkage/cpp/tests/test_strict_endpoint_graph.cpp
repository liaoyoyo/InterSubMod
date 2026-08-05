#include "../strict_endpoint_graph.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace graph = InterSubMod::StrictEndpointGraph;

namespace {

graph::Molecule molecule(const std::string& id, const std::string& ps, const std::string& hp,
                         const std::vector<std::int64_t>& indices, const std::vector<std::int64_t>& positions,
                         const std::string& codes) {
    assert(indices.size() == positions.size());
    assert(indices.size() == codes.size());
    graph::Molecule result{"TEST", "chr1", id, hp, ps, {}};
    for (std::size_t index = 0; index < indices.size(); ++index) {
        result.observations.push_back({{indices[index], positions[index]}, codes[index]});
    }
    return result;
}

std::vector<graph::Molecule> fixture() {
    return {
        // A-B and B-C each have two distinct molecules, so transitivity yields A-B-C.
        molecule("m01", "100", "HP1", {0, 1}, {100, 200}, "RA"),
        molecule("m02", "100", "1", {0, 1}, {100, 200}, "AA"),
        molecule("m03", "100", "HP1", {1, 2}, {200, 300}, "RA"),
        molecule("m04", "100", "1", {1, 2}, {200, 300}, "AA"),

        // A-X-C supports A-C only; node B must remain isolated.
        molecule("m05", "100", "HP1", {3, 4, 5}, {400, 500, 600}, "AXR"),
        molecule("m06", "100", "1", {3, 4, 5}, {400, 500, 600}, "RXA"),

        // Same HP but different exact PS values must remain separate containers.
        molecule("m07", "200", "HP1", {6, 7}, {700, 800}, "RR"),
        molecule("m08", "200", "1", {6, 7}, {700, 800}, "AA"),
        molecule("m09", "300", "HP1", {7, 8}, {800, 900}, "RA"),
        molecule("m10", "300", "1", {7, 8}, {800, 900}, "AR"),

        // Physical distance is not a linkage veto when two molecules support both endpoints.
        molecule("m11", "400", "HP2", {9, 10}, {1000, 101000}, "RA"),
        molecule("m12", "400", "2", {9, 10}, {1000, 101000}, "AA"),

        // These rows are deliberately outside the primary exact-PS x HP1/2 contract.
        molecule("m13", "500", "HP3", {11, 12}, {1100, 1200}, "AA"),
        molecule("m14", ".", "HP1", {13, 14}, {1300, 1400}, "AA"),
    };
}

const graph::EdgeRecord& find_edge(const graph::BuildResult& result, const std::string& ps, std::int64_t left,
                                   std::int64_t right) {
    for (const auto& edge : result.edges) {
        if (edge.container.phase_set == ps && edge.edge.left.site_index == left &&
            edge.edge.right.site_index == right) {
            return edge;
        }
    }
    throw std::runtime_error("edge not found");
}

bool has_component(const graph::BuildResult& result, const std::string& ps,
                   const std::vector<std::int64_t>& expected_indices) {
    for (const auto& component : result.components) {
        if (component.container.phase_set != ps) continue;
        std::vector<std::int64_t> indices;
        for (const auto& node : component.nodes) indices.push_back(node.site_index);
        if (indices == expected_indices) return true;
    }
    return false;
}

}  // namespace

int main() {
    const std::vector<std::string> missing_tokens = {"", ".", "NA", "N/A", "NaN", "None", " NULL "};
    for (const std::string& token : missing_tokens) {
        assert(graph::normalize_phase_set(token).empty());
    }
    assert(graph::normalize_phase_set(" 123 ") == "123");
    std::vector<graph::Molecule> missing_phase_sets;
    std::size_t missing_index = 0;
    for (const std::string& token : missing_tokens) {
        missing_phase_sets.push_back(
            molecule("missing-" + std::to_string(missing_index++), token, "HP1", {0, 1}, {100, 200}, "AA"));
    }
    const graph::BuildResult missing_result = graph::build_strict_endpoint_graph(missing_phase_sets, 1);
    assert(missing_result.stats.input_molecules == 7);
    assert(missing_result.stats.eligible_molecules == 0);
    assert(missing_result.stats.skipped_missing_phase_set == 7);
    assert(missing_result.stats.containers == 0);
    assert(missing_result.components.empty());

    const std::vector<graph::Molecule> molecules = fixture();
    const graph::BuildResult threshold2 = graph::build_strict_endpoint_graph(molecules, 2);
    assert(threshold2.stats.input_molecules == 14);
    assert(threshold2.stats.eligible_molecules == 12);
    assert(threshold2.stats.skipped_nonprimary_hp == 1);
    assert(threshold2.stats.skipped_missing_phase_set == 1);
    assert(threshold2.stats.containers == 4);
    assert(threshold2.stats.observed_edges == 6);
    assert(threshold2.stats.qualifying_edges == 6);
    assert(threshold2.stats.components == 5);
    assert(threshold2.stats.singleton_components == 0);
    assert(threshold2.stats.multisite_components == 5);

    assert(has_component(threshold2, "100", {0, 1, 2}));
    assert(has_component(threshold2, "100", {3, 5}));
    assert(!has_component(threshold2, "100", {3, 4, 5}));
    assert(!has_component(threshold2, "100", {4}));
    assert(has_component(threshold2, "200", {6, 7}));
    assert(has_component(threshold2, "300", {7, 8}));
    assert(has_component(threshold2, "400", {9, 10}));

    const graph::EdgeRecord& ab = find_edge(threshold2, "100", 0, 1);
    assert(ab.support.total == 2 && ab.support.ra == 1 && ab.support.aa == 1 && ab.qualifies);
    const graph::EdgeRecord& ac = find_edge(threshold2, "100", 3, 5);
    assert(ac.support.total == 2 && ac.support.ra == 1 && ac.support.ar == 1 && ac.qualifies);
    for (const auto& edge : threshold2.edges) {
        assert(!(edge.container.phase_set == "100" &&
                 (edge.edge.left.site_index == 4 || edge.edge.right.site_index == 4)));
    }

    const graph::BuildResult threshold3 = graph::build_strict_endpoint_graph(molecules, 3);
    assert(threshold3.stats.observed_edges == 6);
    assert(threshold3.stats.qualifying_edges == 0);
    assert(threshold3.stats.container_nodes == 11);
    assert(threshold3.stats.components == 11);
    assert(threshold3.stats.singleton_components == 11);
    assert(threshold3.stats.multisite_components == 0);

    std::vector<graph::Molecule> reversed = molecules;
    std::reverse(reversed.begin(), reversed.end());
    const graph::BuildResult reversed_result = graph::build_strict_endpoint_graph(reversed, 2);
    assert(reversed_result.graph_digest_fnv1a64 == threshold2.graph_digest_fnv1a64);

    std::vector<graph::Molecule> duplicate = molecules;
    duplicate.push_back(molecules.front());
    bool duplicate_failed = false;
    try {
        (void)graph::build_strict_endpoint_graph(duplicate, 2);
    } catch (const std::invalid_argument& error) {
        duplicate_failed = std::string(error.what()).find("duplicate molecule_id") != std::string::npos;
    }
    assert(duplicate_failed);

    std::cout << "test_strict_endpoint_graph: PASS\n"
              << "threshold2: containers=4 observed_edges=6 qualifying_edges=6 components=5\n"
              << "threshold3: qualifying_edges=0 components=11\n"
              << "digest=" << threshold2.graph_digest_fnv1a64 << "\n";
    return 0;
}
