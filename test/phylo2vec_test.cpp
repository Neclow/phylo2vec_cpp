#include "../src/phylo2vec.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <unordered_map>

const int MIN_K = 3;
const int NUM_TESTS = 100;

class TestCountTracker : public testing::EmptyTestEventListener {
   public:
    void OnTestStart(const testing::TestInfo& test_info) override {
        std::string test_name = test_info.name();
        test_name = test_name.substr(0, test_name.find("/", 0));
        testCounts[test_info.test_case_name() + std::string(".") + test_name]++;
    }

    void OnTestEnd(const testing::TestInfo& test_info) override {
        if (test_info.result()->Passed()) {
            std::string test_name = test_info.name();
            test_name = test_name.substr(0, test_name.find("/", 0));
            testSuccesses[test_info.test_case_name() + std::string(".") + test_name]++;
        }
    }

    void OnTestProgramEnd(const testing::UnitTest& unit_test) override {
        std::cout << "\nDetails:\n";
        for (const auto& testCount : testCounts) {
            int successes = testSuccesses.find(testCount.first) == testSuccesses.end()
                                ? 0
                                : testSuccesses[testCount.first];
            std::cout << successes << "/" << testCount.second << " passed tests from "
                      << testCount.first << std::endl;
        }
    }

   private:
    std::unordered_map<std::string, int> testCounts;
    std::unordered_map<std::string, int> testSuccesses;
};

class Phylo2VecTest : public ::testing::TestWithParam<int> {
   protected:
};

INSTANTIATE_TEST_SUITE_P(RandomTests, Phylo2VecTest, ::testing::Range(MIN_K, NUM_TESTS + MIN_K));

TEST_P(Phylo2VecTest, TestSamplingRandomV) {
    int k = GetParam();

    std::vector<int> v = sample(k);

    // Check that it's valid
    EXPECT_NO_THROW(check_v(v));
}

TEST_P(Phylo2VecTest, TestVtoIntNewickBacktoV) {
    int k = GetParam();

    std::vector<int> v = sample(k);

    // Convert to Newick
    std::string nw = toNewick(v);

    // Convert back to v
    Newick2VResult res = newick2v(nw, k + 1);

    std::vector<int> converted_v = res.v;

    converted_v.erase(converted_v.begin());

    EXPECT_EQ(v, converted_v);

    // Replace newick integers with tip_0, tip1, ...
    // tmp = newick2vWithMapping(nw, k + 1)
    // converted_v_from_taxon_newick
    // EXPECT_EQ(v, converted_v_from_taxon_newick
}

TEST_P(Phylo2VecTest, TestGetNumLeavesFromNewick) {
    int k = GetParam();

    std::vector<int> v = sample(k);

    // Convert to Newick
    std::string nw = toNewick(v);

    processNewick(nw);

    EXPECT_EQ(getNumLeavesFromNewick(nw), k + 1);
}

TEST(StringNewickTest, TestStringNewickToV) {
    std::ifstream file("../test/100trees.txt");

    // Newick where nodes are strings and not just ints
    // The nodes are all called "tip_1", "tip_2"...
    std::string stringNewick;
    std::regex pattern("tip_");
    int line = 0;
    while (std::getline(file, stringNewick)) {
        std::string intNewick = std::regex_replace(stringNewick, pattern, "");

        processNewick(intNewick);

        int num_leaves = getNumLeavesFromNewick(intNewick);

        std::vector<int> converted_v_from_int_newick = newick2v(intNewick, num_leaves).v;

        std::vector<int> converted_v_from_string_newick =
            newick2vWithMapping(stringNewick, num_leaves).v;

        EXPECT_EQ(converted_v_from_int_newick, converted_v_from_string_newick);

        line++;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Add tracker for each test to get a proper breakdown
    TestCountTracker* tracker = new TestCountTracker;
    ::testing::UnitTest::GetInstance()->listeners().Append(tracker);

    return RUN_ALL_TESTS();
}