#include "strict_endpoint_graph.hpp"

#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace graph = InterSubMod::StrictEndpointGraph;

namespace {

struct Options {
    fs::path input;
    fs::path edges_output;
    fs::path components_output;
    fs::path receipt_output;
    std::uint64_t threshold = 3;
};

std::string usage() {
    return
        "Usage: strict_endpoint_graph_verify --input MOLECULES.tsv --threshold N "
        "--edges-output EDGES.tsv --components-output COMPONENTS.tsv "
        "--receipt-output RECEIPT.json\n\n"
        "Required TSV columns: dataset, chrom, molecule_id, hp_family, phase_set, "
        "site_indices, positions1, call_codes.\n"
        "Use --input - to read an uncompressed TSV stream (for example gzip -dc file.tsv.gz | ...).\n"
        "Only exact non-missing PS and HP1/HP2 rows enter the graph. Only fixed "
        "R/A endpoint pairs from a unique molecule add one unit of edge support.\n";
}

std::uint64_t parse_uint64(const std::string& text, const std::string& label) {
    if (text.empty()) throw std::invalid_argument(label + " must not be empty");
    std::size_t consumed = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(label + " is not an unsigned integer: " + text);
    }
    if (consumed != text.size()) throw std::invalid_argument(label + " has trailing characters: " + text);
    return static_cast<std::uint64_t>(value);
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            std::cout << usage();
            std::exit(0);
        }
        if (index + 1 >= argc) throw std::invalid_argument("missing value after " + argument);
        const std::string value = argv[++index];
        if (argument == "--input") {
            options.input = value;
        } else if (argument == "--edges-output") {
            options.edges_output = value;
        } else if (argument == "--components-output") {
            options.components_output = value;
        } else if (argument == "--receipt-output") {
            options.receipt_output = value;
        } else if (argument == "--threshold") {
            options.threshold = parse_uint64(value, "threshold");
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }
    if (options.input.empty() || options.edges_output.empty() || options.components_output.empty() ||
        options.receipt_output.empty()) {
        throw std::invalid_argument("--input, --edges-output, --components-output, and --receipt-output are required");
    }
    if (options.threshold == 0) throw std::invalid_argument("threshold must be >= 1");
    return options;
}

std::vector<std::string> split_preserving_empty(const std::string& text, char delimiter) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        const std::size_t end = text.find(delimiter, start);
        if (end == std::string::npos) {
            parts.push_back(text.substr(start));
            break;
        }
        parts.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    return parts;
}

std::int64_t parse_int64(const std::string& text, const std::string& label, std::size_t line_number) {
    if (text.empty()) {
        throw std::invalid_argument("line " + std::to_string(line_number) + ": empty " + label);
    }
    std::size_t consumed = 0;
    long long value = 0;
    try {
        value = std::stoll(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("line " + std::to_string(line_number) + ": invalid " + label + ": " + text);
    }
    if (consumed != text.size()) {
        throw std::invalid_argument("line " + std::to_string(line_number) + ": trailing characters in " + label);
    }
    return static_cast<std::int64_t>(value);
}

std::vector<std::int64_t> parse_csv_int64(const std::string& text, const std::string& label,
                                          std::size_t line_number) {
    if (text.empty()) return {};
    std::vector<std::int64_t> result;
    for (const std::string& item : split_preserving_empty(text, ',')) {
        result.push_back(parse_int64(item, label, line_number));
    }
    return result;
}

std::vector<graph::Molecule> read_molecules(const fs::path& input) {
    std::ifstream file_stream;
    std::istream* stream = nullptr;
    if (input == fs::path("-")) {
        stream = &std::cin;
    } else {
        file_stream.open(input);
        if (!file_stream) throw std::runtime_error("cannot open input: " + input.string());
        stream = &file_stream;
    }

    std::string header_line;
    if (!std::getline(*stream, header_line)) throw std::runtime_error("input TSV is empty: " + input.string());
    if (!header_line.empty() && header_line.back() == '\r') header_line.pop_back();
    const std::vector<std::string> header = split_preserving_empty(header_line, '\t');
    std::map<std::string, std::size_t> column;
    for (std::size_t index = 0; index < header.size(); ++index) {
        if (header[index].empty()) throw std::runtime_error("input TSV has an empty header field");
        if (!column.emplace(header[index], index).second) {
            throw std::runtime_error("input TSV has duplicate header: " + header[index]);
        }
    }
    const std::vector<std::string> required = {"dataset",      "chrom",      "molecule_id", "hp_family",
                                               "phase_set",    "site_indices", "positions1",  "call_codes"};
    for (const std::string& field : required) {
        if (column.count(field) == 0) throw std::runtime_error("input TSV missing required column: " + field);
    }

    std::vector<graph::Molecule> molecules;
    std::string line;
    std::size_t line_number = 1;
    while (std::getline(*stream, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) throw std::runtime_error("line " + std::to_string(line_number) + ": blank rows are not allowed");
        const std::vector<std::string> fields = split_preserving_empty(line, '\t');
        if (fields.size() != header.size()) {
            throw std::runtime_error("line " + std::to_string(line_number) + ": expected " +
                                     std::to_string(header.size()) + " TSV fields, found " +
                                     std::to_string(fields.size()));
        }
        const auto get = [&](const std::string& name) -> const std::string& { return fields[column.at(name)]; };
        const std::vector<std::int64_t> indices = parse_csv_int64(get("site_indices"), "site_indices", line_number);
        const std::vector<std::int64_t> positions = parse_csv_int64(get("positions1"), "positions1", line_number);
        const std::string& codes = get("call_codes");
        if (indices.size() != positions.size() || indices.size() != codes.size()) {
            throw std::runtime_error("line " + std::to_string(line_number) +
                                     ": site_indices/positions1/call_codes length mismatch");
        }
        graph::Molecule molecule;
        molecule.dataset = get("dataset");
        molecule.chrom = get("chrom");
        molecule.molecule_id = get("molecule_id");
        molecule.hp_family = get("hp_family");
        molecule.phase_set = get("phase_set");
        for (std::size_t index = 0; index < indices.size(); ++index) {
            molecule.observations.push_back(graph::Observation{{indices[index], positions[index]}, codes[index]});
        }
        molecules.push_back(std::move(molecule));
    }
    if (!stream->eof()) throw std::runtime_error("I/O error while reading input: " + input.string());
    return molecules;
}

void ensure_parent_directory(const fs::path& path) {
    if (path.has_parent_path()) fs::create_directories(path.parent_path());
}

std::string json_escape(const std::string& value) {
    std::string result;
    for (const unsigned char byte : value) {
        switch (byte) {
            case '"': result += "\\\""; break;
            case '\\': result += "\\\\"; break;
            case '\b': result += "\\b"; break;
            case '\f': result += "\\f"; break;
            case '\n': result += "\\n"; break;
            case '\r': result += "\\r"; break;
            case '\t': result += "\\t"; break;
            default:
                if (byte < 0x20U) {
                    constexpr char kHex[] = "0123456789abcdef";
                    result += "\\u00";
                    result += kHex[(byte >> 4U) & 0x0FU];
                    result += kHex[byte & 0x0FU];
                } else {
                    result += static_cast<char>(byte);
                }
        }
    }
    return result;
}

void write_edges(const fs::path& output, const graph::BuildResult& result) {
    ensure_parent_directory(output);
    std::ofstream stream(output);
    if (!stream) throw std::runtime_error("cannot open edge output: " + output.string());
    stream << "dataset\tchrom\tphase_set\thp_family\tleft_site_index\tright_site_index\tleft_pos1\tright_pos1"
              "\tsupport_total\tsupport_RR\tsupport_RA\tsupport_AR\tsupport_AA\tthreshold\tqualifies\n";
    for (const graph::EdgeRecord& edge : result.edges) {
        stream << edge.container.dataset << '\t' << edge.container.chrom << '\t' << edge.container.phase_set << '\t'
               << edge.container.hp_family << '\t' << edge.edge.left.site_index << '\t'
               << edge.edge.right.site_index << '\t' << edge.edge.left.position1 << '\t'
               << edge.edge.right.position1 << '\t' << edge.support.total << '\t' << edge.support.rr << '\t'
               << edge.support.ra << '\t' << edge.support.ar << '\t' << edge.support.aa << '\t' << result.threshold
               << '\t' << (edge.qualifies ? 1 : 0) << '\n';
    }
    if (!stream) throw std::runtime_error("failed while writing edge output: " + output.string());
}

void write_components(const fs::path& output, const graph::BuildResult& result) {
    ensure_parent_directory(output);
    std::ofstream stream(output);
    if (!stream) throw std::runtime_error("cannot open component output: " + output.string());
    stream << "dataset\tchrom\tphase_set\thp_family\tcomponent_id\tthreshold\tk\tsite_indices\tpositions1"
              "\tqualifying_edge_count\tcomponent_class\ttree_eligible\tinference_role\n";
    for (const graph::ComponentRecord& component : result.components) {
        stream << component.container.dataset << '\t' << component.container.chrom << '\t'
               << component.container.phase_set << '\t' << component.container.hp_family << '\t'
               << component.component_id << '\t' << result.threshold << '\t' << component.nodes.size() << '\t'
               << graph::join_site_indices(component.nodes) << '\t' << graph::join_positions1(component.nodes) << '\t'
               << component.qualifying_edge_count << '\t'
               << (component.nodes.size() > 1 ? "READ_LINKED_MULTISITE_REGION" : "UNLINKED_SINGLETON_COMPONENT")
               << '\t' << (component.nodes.size() > 1 ? "true" : "false") << '\t'
               << (component.nodes.size() > 1 ? "PRIMARY_PS_AWARE" : "ABSTAIN_SINGLETON_UNLINKED") << '\n';
    }
    if (!stream) throw std::runtime_error("failed while writing component output: " + output.string());
}

void write_receipt(const fs::path& output, const Options& options, const graph::BuildResult& result) {
    ensure_parent_directory(output);
    std::ofstream stream(output);
    if (!stream) throw std::runtime_error("cannot open receipt output: " + output.string());
    const graph::BuildStats& stats = result.stats;
    stream << "{\n"
           << "  \"schema_name\": \"intersubmod.strict_endpoint_graph_receipt\",\n"
           << "  \"schema_version\": \"1.1.0\",\n"
           << "  \"linkage_rule\": \"distinct_canonical_molecule_fixed_RA_endpoint_pair\",\n"
           << "  \"container_rule\": \"dataset_x_chromosome_x_exact_nonmissing_PS_x_HP1_or_HP2\",\n"
           << "  \"input\": \"" << json_escape(options.input.string()) << "\",\n"
           << "  \"edges_output\": \"" << json_escape(options.edges_output.string()) << "\",\n"
           << "  \"components_output\": \"" << json_escape(options.components_output.string()) << "\",\n"
           << "  \"threshold\": " << result.threshold << ",\n"
           << "  \"graph_digest_algorithm\": \"fnv1a64_canonical_records\",\n"
           << "  \"graph_digest\": \"" << result.graph_digest_fnv1a64 << "\",\n"
           << "  \"counts\": {\n"
           << "    \"input_molecules\": " << stats.input_molecules << ",\n"
           << "    \"eligible_molecules\": " << stats.eligible_molecules << ",\n"
           << "    \"skipped_nonprimary_hp\": " << stats.skipped_nonprimary_hp << ",\n"
           << "    \"skipped_missing_phase_set\": " << stats.skipped_missing_phase_set << ",\n"
           << "    \"fixed_calls\": " << stats.fixed_calls << ",\n"
           << "    \"nonlinking_calls\": " << stats.nonlinking_calls << ",\n"
           << "    \"containers\": " << stats.containers << ",\n"
           << "    \"container_nodes\": " << stats.container_nodes << ",\n"
           << "    \"observed_edges\": " << stats.observed_edges << ",\n"
           << "    \"qualifying_edges\": " << stats.qualifying_edges << ",\n"
           << "    \"components\": " << stats.components << ",\n"
           << "    \"singleton_components\": " << stats.singleton_components << ",\n"
           << "    \"multisite_components\": " << stats.multisite_components << "\n"
           << "  },\n"
           << "  \"invariants\": {\n"
           << "    \"molecule_id_global_unique\": true,\n"
           << "    \"only_exact_nonmissing_PS_HP1_HP2\": true,\n"
           << "    \"only_fixed_RA_endpoints_add_support\": true,\n"
           << "    \"edge_support_is_distinct_molecule_count\": true,\n"
           << "    \"components_use_only_qualifying_edges\": true,\n"
           << "    \"singleton_components_are_abstain_not_tree_regions\": true\n"
           << "  }\n"
           << "}\n";
    if (!stream) throw std::runtime_error("failed while writing receipt output: " + output.string());
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        const std::vector<graph::Molecule> molecules = read_molecules(options.input);
        const graph::BuildResult result = graph::build_strict_endpoint_graph(molecules, options.threshold);
        write_edges(options.edges_output, result);
        write_components(options.components_output, result);
        write_receipt(options.receipt_output, options, result);
        std::cout << "strict_endpoint_graph_verify: PASS\n"
                  << "input=" << options.input << "\n"
                  << "threshold=" << result.threshold << "\n"
                  << "eligible_molecules=" << result.stats.eligible_molecules << "\n"
                  << "containers=" << result.stats.containers << "\n"
                  << "qualifying_edges=" << result.stats.qualifying_edges << "\n"
                  << "components=" << result.stats.components << "\n"
                  << "graph_digest_fnv1a64=" << result.graph_digest_fnv1a64 << "\n"
                  << "edges_output=" << options.edges_output << "\n"
                  << "components_output=" << options.components_output << "\n"
                  << "receipt_output=" << options.receipt_output << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "strict_endpoint_graph_verify: FAIL: " << error.what() << '\n';
        return 2;
    }
}
