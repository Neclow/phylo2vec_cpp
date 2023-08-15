#include <iostream>

#include "cxxopts.hpp"
#include "phylo2vec.hpp"

cxxopts::Options get_options() {
    // Parse options with CXXOpts
    cxxopts::Options options("phylo2vec",
                             "Convert Newick strings to integer vectors and vice-versa");

    // clang-format off
    options.add_options()
        ("h,help", "Show help")
        ("toNewick", "Convert to Newick format. Example input: 0 1 4",cxxopts::value<std::vector<int>>())
        ("toVector", "Convert to integer vector. Example input: \"(((2,1)4,0)5,3)6;\"", cxxopts::value<std::string>())
        ("with_mapping", "For Newicks that do not only contain digits.", cxxopts::value<bool>()->default_value("false"))
        ("num_leaves", "Number of leaves (optional, but recommended when using toVector)", cxxopts::value<int>());
    // clang-format on

    options.positional_help("toNewick toVector");

    options.parse_positional({"toNewick", "toVector"});

    options.show_positional_help();

    return options;
}

std::map<int, std::string> convertMapping(std::map<std::string, std::string> mapping) {
    std::map<int, std::string> converted_mapping;
    for (const auto& pair : mapping) {
        int key = std::stoi(pair.first);
        converted_mapping[key] = pair.second;
    }

    return converted_mapping;
}

void doToNewick(std::vector<int> v) {
    std::string newick = toNewick(v);
    std::cout << "Newick string: " << newick << std::endl;
}

void doToVector(std::string newick, int num_leaves, bool with_mapping) {
    std::vector<int> converted_v;

    if (with_mapping) {
        Newick2VResult tmp = newick2vWithMapping(newick, num_leaves);
        converted_v = tmp.v;

        std::map<int, std::string> converted_mapping = convertMapping(tmp.mapping);

        std::cout << "Number of leaves: " << tmp.num_leaves << std::endl;

        std::cout << "Mapping:" << std::endl;
        for (const auto& elem : converted_mapping) {
            std::cout << elem.first << "->" << elem.second << std::endl;
        }
    } else {
        converted_v = newick2v(newick, num_leaves).v;
    }

    std::cout << "Integer vector: ";
    for (int val : converted_v) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    cxxopts::Options options = get_options();

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (result.count("toNewick")) {
        std::vector<int> v = result["toNewick"].as<std::vector<int>>();
        doToNewick(v);
    } else if (result.count("toVector")) {
        std::string newick = result["toVector"].as<std::string>();
        int num_leaves = result.count("num_leaves") ? result["num_leaves"].as<int>() : -1;

        bool with_mapping = result["with_mapping"].as<bool>();
        doToVector(newick, num_leaves, with_mapping);
    } else {
        std::cerr << "Invalid arguments. Use --help for usage information." << std::endl;
        return 1;
    }

    return 0;
}