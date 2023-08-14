#include <array>
#include <iostream>
#include <vector>

#include "cxxopts.hpp"
#include "phylo2vec.hpp"

int main(int argc, char* argv[]) {
    cxxopts::Options options("phylo2vec",
                             "Convert Newick strings to integer vectors and vice-versa");

    // clang-format off
    options.add_options()
        ("h,help", "Show help")
        ("toNewick", "Convert to Newick format. Example input: 0 1 4",cxxopts::value<std::vector<int>>())
        ("toVector", "Convert to integer vector. Example input: \"(((2,1)4,0)5,3)6;\"", cxxopts::value<std::string>())
        ("num_leaves", "Number of leaves (optional, but recommended when using toVector)", cxxopts::value<int>());
    // clang-format on

    options.positional_help("toNewick toVector");

    options.parse_positional({"toNewick", "toVector"});

    options.show_positional_help();

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (result.count("toNewick")) {
        auto v = result["toNewick"].as<std::vector<int>>();
        std::string newick = toNewick(v);
        std::cout << "Newick string: " << newick << std::endl;
    } else if (result.count("toVector")) {
        std::string newick = result["toVector"].as<std::string>();
        int num_leaves = result.count("num_leaves") ? result["num_leaves"].as<int>() : -1;
        std::vector<int> converted_v = newick2v(newick, num_leaves);
        std::cout << "Integer vector: ";
        for (int val : converted_v) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    } else {
        std::cerr << "Invalid arguments. Use --help for usage information." << std::endl;
        return 1;
    }

    return 0;
}