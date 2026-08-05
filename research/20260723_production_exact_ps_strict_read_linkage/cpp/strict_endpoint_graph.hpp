#ifndef INTERSUBMOD_STRICT_ENDPOINT_GRAPH_HPP_
#define INTERSUBMOD_STRICT_ENDPOINT_GRAPH_HPP_

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace InterSubMod::StrictEndpointGraph {

struct Node {
    std::int64_t site_index = -1;
    std::int64_t position1 = -1;

    bool operator<(const Node& other) const {
        return std::tie(position1, site_index) < std::tie(other.position1, other.site_index);
    }
    bool operator==(const Node& other) const {
        return site_index == other.site_index && position1 == other.position1;
    }
};

struct ContainerKey {
    std::string dataset;
    std::string chrom;
    std::string phase_set;
    std::string hp_family;

    bool operator<(const ContainerKey& other) const {
        return std::tie(dataset, chrom, phase_set, hp_family) <
               std::tie(other.dataset, other.chrom, other.phase_set, other.hp_family);
    }
    bool operator==(const ContainerKey& other) const {
        return dataset == other.dataset && chrom == other.chrom && phase_set == other.phase_set &&
               hp_family == other.hp_family;
    }
};

struct Observation {
    Node node;
    char call_code = 'X';
};

struct Molecule {
    std::string dataset;
    std::string chrom;
    std::string molecule_id;
    std::string hp_family;
    std::string phase_set;
    std::vector<Observation> observations;
};

struct EdgeKey {
    Node left;
    Node right;

    bool operator<(const EdgeKey& other) const {
        return std::tie(left, right) < std::tie(other.left, other.right);
    }
};

struct EdgeSupport {
    std::uint64_t total = 0;
    std::uint64_t rr = 0;
    std::uint64_t ra = 0;
    std::uint64_t ar = 0;
    std::uint64_t aa = 0;

    void add(char left_call, char right_call) {
        ++total;
        if (left_call == 'R' && right_call == 'R') {
            ++rr;
        } else if (left_call == 'R' && right_call == 'A') {
            ++ra;
        } else if (left_call == 'A' && right_call == 'R') {
            ++ar;
        } else if (left_call == 'A' && right_call == 'A') {
            ++aa;
        } else {
            throw std::logic_error("edge state must have fixed R/A endpoints");
        }
    }
};

struct EdgeRecord {
    ContainerKey container;
    EdgeKey edge;
    EdgeSupport support;
    bool qualifies = false;
};

struct ComponentRecord {
    ContainerKey container;
    std::string component_id;
    std::vector<Node> nodes;
    std::uint64_t qualifying_edge_count = 0;
};

struct BuildStats {
    std::uint64_t input_molecules = 0;
    std::uint64_t eligible_molecules = 0;
    std::uint64_t skipped_nonprimary_hp = 0;
    std::uint64_t skipped_missing_phase_set = 0;
    std::uint64_t nonlinking_calls = 0;
    std::uint64_t fixed_calls = 0;
    std::uint64_t containers = 0;
    std::uint64_t container_nodes = 0;
    std::uint64_t observed_edges = 0;
    std::uint64_t qualifying_edges = 0;
    std::uint64_t components = 0;
    std::uint64_t singleton_components = 0;
    std::uint64_t multisite_components = 0;
};

struct BuildResult {
    std::uint64_t threshold = 0;
    std::vector<EdgeRecord> edges;
    std::vector<ComponentRecord> components;
    BuildStats stats;
    std::string graph_digest_fnv1a64;
};

inline bool is_fixed_call(char code) { return code == 'R' || code == 'A'; }

inline bool is_valid_call(char code) {
    return code == 'R' || code == 'A' || code == 'X' || code == 'L' || code == 'O' || code == 'D' ||
           code == 'S';
}

inline std::string normalize_primary_hp(const std::string& hp) {
    if (hp == "1" || hp == "HP1") return "HP1";
    if (hp == "2" || hp == "HP2") return "HP2";
    return "";
}

inline std::string normalize_phase_set(const std::string& phase_set) {
    const auto first = std::find_if_not(phase_set.begin(), phase_set.end(), [](unsigned char value) {
        return std::isspace(value) != 0;
    });
    const auto last = std::find_if_not(phase_set.rbegin(), phase_set.rend(), [](unsigned char value) {
                          return std::isspace(value) != 0;
                      }).base();
    const std::string trimmed = first < last ? std::string(first, last) : std::string();
    std::string upper;
    upper.reserve(trimmed.size());
    std::transform(trimmed.begin(), trimmed.end(), std::back_inserter(upper), [](unsigned char value) {
        return static_cast<char>(std::toupper(value));
    });
    if (upper.empty() || upper == "." || upper == "NA" || upper == "N/A" || upper == "NAN" ||
        upper == "NONE" || upper == "NULL") {
        return "";
    }
    return trimmed;
}

namespace detail {

class DisjointSet {
  public:
    explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0) {
        for (std::size_t index = 0; index < size; ++index) parent_[index] = index;
    }

    std::size_t find(std::size_t item) {
        if (parent_[item] != item) parent_[item] = find(parent_[item]);
        return parent_[item];
    }

    void unite(std::size_t left, std::size_t right) {
        left = find(left);
        right = find(right);
        if (left == right) return;
        if (rank_[left] < rank_[right]) std::swap(left, right);
        parent_[right] = left;
        if (rank_[left] == rank_[right]) ++rank_[left];
    }

  private:
    std::vector<std::size_t> parent_;
    std::vector<std::uint8_t> rank_;
};

class Fnv1a64 {
  public:
    void update(const std::string& value) {
        for (const unsigned char byte : value) {
            value_ ^= static_cast<std::uint64_t>(byte);
            value_ *= 1099511628211ULL;
        }
    }

    std::string hex() const {
        constexpr char kHex[] = "0123456789abcdef";
        std::string result(16, '0');
        std::uint64_t value = value_;
        for (int index = 15; index >= 0; --index) {
            result[static_cast<std::size_t>(index)] = kHex[value & 0x0FULL];
            value >>= 4U;
        }
        return result;
    }

  private:
    std::uint64_t value_ = 14695981039346656037ULL;
};

inline std::string join_node_indices(const std::vector<Node>& nodes) {
    std::string result;
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        if (index != 0) result += ',';
        result += std::to_string(nodes[index].site_index);
    }
    return result;
}

inline std::string join_node_positions(const std::vector<Node>& nodes) {
    std::string result;
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        if (index != 0) result += ',';
        result += std::to_string(nodes[index].position1);
    }
    return result;
}

}  // namespace detail

inline BuildResult build_strict_endpoint_graph(const std::vector<Molecule>& molecules, std::uint64_t threshold) {
    if (threshold == 0) throw std::invalid_argument("threshold must be >= 1");

    BuildResult result;
    result.threshold = threshold;
    result.stats.input_molecules = molecules.size();

    std::set<std::string> molecule_ids;
    std::map<std::tuple<std::string, std::string, std::int64_t>, std::int64_t> known_positions;
    std::map<ContainerKey, std::set<Node>> nodes_by_container;
    std::map<ContainerKey, std::map<EdgeKey, EdgeSupport>> edges_by_container;

    for (const Molecule& molecule : molecules) {
        if (molecule.molecule_id.empty()) throw std::invalid_argument("molecule_id must not be empty");
        if (!molecule_ids.insert(molecule.molecule_id).second) {
            throw std::invalid_argument("duplicate molecule_id: " + molecule.molecule_id);
        }
        if (molecule.dataset.empty() || molecule.chrom.empty()) {
            throw std::invalid_argument("dataset and chrom must not be empty for molecule " + molecule.molecule_id);
        }

        const std::string hp = normalize_primary_hp(molecule.hp_family);
        if (hp.empty()) {
            ++result.stats.skipped_nonprimary_hp;
            continue;
        }
        const std::string phase_set = normalize_phase_set(molecule.phase_set);
        if (phase_set.empty()) {
            ++result.stats.skipped_missing_phase_set;
            continue;
        }
        ++result.stats.eligible_molecules;

        ContainerKey container{molecule.dataset, molecule.chrom, phase_set, hp};
        std::set<std::int64_t> row_site_indices;
        std::vector<std::pair<Node, char>> fixed;
        for (const Observation& observation : molecule.observations) {
            if (observation.node.site_index < 0 || observation.node.position1 <= 0) {
                throw std::invalid_argument("site_index must be >=0 and position1 >=1 for molecule " +
                                            molecule.molecule_id);
            }
            if (!row_site_indices.insert(observation.node.site_index).second) {
                throw std::invalid_argument("duplicate site_index within molecule " + molecule.molecule_id);
            }
            if (!is_valid_call(observation.call_code)) {
                throw std::invalid_argument("unknown call code in molecule " + molecule.molecule_id);
            }
            const auto site_key = std::make_tuple(molecule.dataset, molecule.chrom, observation.node.site_index);
            const auto known = known_positions.find(site_key);
            if (known == known_positions.end()) {
                known_positions.emplace(site_key, observation.node.position1);
            } else if (known->second != observation.node.position1) {
                throw std::invalid_argument("inconsistent position for site_index " +
                                            std::to_string(observation.node.site_index));
            }
            if (is_fixed_call(observation.call_code)) {
                nodes_by_container[container].insert(observation.node);
                fixed.emplace_back(observation.node, observation.call_code);
                ++result.stats.fixed_calls;
            } else {
                ++result.stats.nonlinking_calls;
            }
        }
        std::sort(fixed.begin(), fixed.end(), [](const auto& left, const auto& right) { return left.first < right.first; });
        for (std::size_t left_index = 0; left_index < fixed.size(); ++left_index) {
            for (std::size_t right_index = left_index + 1; right_index < fixed.size(); ++right_index) {
                const EdgeKey edge{fixed[left_index].first, fixed[right_index].first};
                edges_by_container[container][edge].add(fixed[left_index].second, fixed[right_index].second);
            }
        }
    }

    result.stats.containers = nodes_by_container.size();
    for (const auto& [container, node_set] : nodes_by_container) {
        const std::vector<Node> nodes(node_set.begin(), node_set.end());
        result.stats.container_nodes += nodes.size();
        std::map<Node, std::size_t> node_to_index;
        for (std::size_t index = 0; index < nodes.size(); ++index) node_to_index.emplace(nodes[index], index);
        detail::DisjointSet dsu(nodes.size());

        const auto edge_map_it = edges_by_container.find(container);
        if (edge_map_it != edges_by_container.end()) {
            for (const auto& [edge, support] : edge_map_it->second) {
                const bool qualifies = support.total >= threshold;
                result.edges.push_back(EdgeRecord{container, edge, support, qualifies});
                ++result.stats.observed_edges;
                if (qualifies) {
                    ++result.stats.qualifying_edges;
                    dsu.unite(node_to_index.at(edge.left), node_to_index.at(edge.right));
                }
            }
        }

        std::map<std::size_t, std::vector<Node>> component_nodes;
        for (std::size_t index = 0; index < nodes.size(); ++index) component_nodes[dsu.find(index)].push_back(nodes[index]);
        std::vector<std::vector<Node>> ordered_components;
        for (auto& [unused_root, members] : component_nodes) {
            (void)unused_root;
            std::sort(members.begin(), members.end());
            ordered_components.push_back(std::move(members));
        }
        std::sort(ordered_components.begin(), ordered_components.end());

        for (std::size_t ordinal = 0; ordinal < ordered_components.size(); ++ordinal) {
            const std::vector<Node>& members = ordered_components[ordinal];
            const std::set<Node> member_set(members.begin(), members.end());
            std::uint64_t edge_count = 0;
            if (edge_map_it != edges_by_container.end()) {
                for (const auto& [edge, support] : edge_map_it->second) {
                    if (support.total >= threshold && member_set.count(edge.left) != 0 &&
                        member_set.count(edge.right) != 0) {
                        ++edge_count;
                    }
                }
            }
            std::string component_id = "C";
            const std::string ordinal_text = std::to_string(ordinal + 1);
            component_id.append(6 - std::min<std::size_t>(6, ordinal_text.size()), '0');
            component_id += ordinal_text;
            result.components.push_back(ComponentRecord{container, component_id, members, edge_count});
            ++result.stats.components;
            if (members.size() == 1) {
                ++result.stats.singleton_components;
            } else {
                ++result.stats.multisite_components;
            }
        }
    }

    std::sort(result.edges.begin(), result.edges.end(), [](const EdgeRecord& left, const EdgeRecord& right) {
        return std::tie(left.container, left.edge) < std::tie(right.container, right.edge);
    });
    std::sort(result.components.begin(), result.components.end(),
              [](const ComponentRecord& left, const ComponentRecord& right) {
                  return std::tie(left.container, left.nodes, left.component_id) <
                         std::tie(right.container, right.nodes, right.component_id);
              });

    detail::Fnv1a64 digest;
    digest.update("strict_endpoint_graph_v1\nthreshold=" + std::to_string(threshold) + "\n");
    for (const EdgeRecord& edge : result.edges) {
        digest.update("E\t" + edge.container.dataset + "\t" + edge.container.chrom + "\t" +
                      edge.container.phase_set + "\t" + edge.container.hp_family + "\t" +
                      std::to_string(edge.edge.left.site_index) + "\t" +
                      std::to_string(edge.edge.right.site_index) + "\t" +
                      std::to_string(edge.support.total) + "\t" + std::to_string(edge.support.rr) + "\t" +
                      std::to_string(edge.support.ra) + "\t" + std::to_string(edge.support.ar) + "\t" +
                      std::to_string(edge.support.aa) + "\t" + (edge.qualifies ? "1\n" : "0\n"));
    }
    for (const ComponentRecord& component : result.components) {
        digest.update("C\t" + component.container.dataset + "\t" + component.container.chrom + "\t" +
                      component.container.phase_set + "\t" + component.container.hp_family + "\t" +
                      detail::join_node_indices(component.nodes) + "\t" +
                      detail::join_node_positions(component.nodes) + "\n");
    }
    result.graph_digest_fnv1a64 = digest.hex();
    return result;
}

inline std::string join_site_indices(const std::vector<Node>& nodes) { return detail::join_node_indices(nodes); }
inline std::string join_positions1(const std::vector<Node>& nodes) { return detail::join_node_positions(nodes); }

}  // namespace InterSubMod::StrictEndpointGraph

#endif  // INTERSUBMOD_STRICT_ENDPOINT_GRAPH_HPP_
