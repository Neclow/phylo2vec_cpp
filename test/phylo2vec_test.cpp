#include "../src/phylo2vec.hpp"

#include <gtest/gtest.h>

#include <random>

const int NUM_TESTS = 100;
const int MIN_K = 3;
const int MAX_K = 50;

class Phylo2VecTest : public ::testing::TestWithParam<int> {
   protected:
    std::mt19937 gen;
    std::uniform_int_distribution<int> distrib;

    Phylo2VecTest() : gen(std::random_device{}()), distrib(MIN_K, MAX_K) {}
};

TEST_P(Phylo2VecTest, TestVtoNewickBacktoV) {
    int k = GetParam();

    std::vector<int> v = sample(k);

    // Check that it's valid
    EXPECT_NO_THROW(check_v(v));

    // Convert to Newick
    std::string nw = toNewick(v);

    // Convert back to v
    std::vector<int> converted_v = newick2v(nw, k + 1);

    converted_v.erase(converted_v.begin());

    EXPECT_EQ(v, converted_v);
}

INSTANTIATE_TEST_SUITE_P(RandomTests, Phylo2VecTest, ::testing::Range(1, 101));

// TEST(Phylo2VecTest, V2M2V) {
//     std::random_device rd;  // a seed source for the random number engine
//     std::mt19937 gen(rd());
//     std::uniform_int_distribution<> distrib(MIN_K, MAX_K);
//     for (int i = 0; i < NUM_TESTS; ++i) {
//         // Generate a random v
//         int k = distrib(gen);
//         std::vector<int> v = sample(k);

//         // Check that it's valid
//         EXPECT_NO_THROW(check_v(v));

//         // Convert to Newick
//         std::string nw = toNewick(v);

//         // Convert back to v
//         std::vector<int> converted_v = newick2v(nw, k + 1);

//         converted_v.erase(converted_v.begin());

//         // Check equality
//         EXPECT_EQ(v, converted_v);
//     }
// }

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}