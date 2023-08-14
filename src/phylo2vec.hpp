#ifndef PHYLO2VEC_HPP
#define PHYLO2VEC_HPP

#include <array>
#include <string>
#include <utility>
#include <vector>

/**
 * @brief Sample a random Phylo2Vec v for n_leaves = k + 1
 */
std::vector<int> sample(const int &k);

/**
 * @brief check that Phylo2Vec v is correct
 * i.e. that each 0 <= v[i] <= 2i
 */

void check_v(const std::vector<int> &v);

/**
 * @brief flip the rows of a 2D matrix
 */
template <typename T>
void flipRows(T &vec);

/**
 * @brief flip the columns of a 2D matrix
 */
template <typename T>
void flipColumns(T &vec);

/**
 * @brief flip rows and/or columns of a 2D matrix
 * Equivalent of np.flip
 * @tparam T
 * @param vec
 * @param axis Axis or axes along which to flip over
 * By default: flip all axes
 */
template <typename T>
void flip(T &vec, int axis = -1);

/**
 * @brief initialise the view matrix for getAncestry
 * In Python: np.tril([np.arange(k+1)]*(k+1))
 *
 * @param k number of leaves - 1
 * @return std::vector<std::vector<int>> matrix of size k, k+1
 * The last row is kept apart to clarify the implementation of getAncestry
 */
std::vector<std::vector<int>> initViewMatrix(const int &k);

/**
 * @brief Get ancestry for each node given a v-representation.
 *
 * @param v Phylo2Vec vector
 * @return std::vector<std::array<int, 3>>
 * 1st column: parent
 * 2nd and 3rd columns: children
 */
std::vector<std::array<int, 3>> getAncestry(const std::vector<int> &v);

/**
 * @brief Build a Newick string from an "ancestry" array to describe a tree
 * M is processed such that we iteratively write a Newick string
 * to describe the tree.
 * @param M cf. getAncestry
 * @return std::string Newick-format representation of a tree
 */
std::string buildNewick(std::vector<std::array<int, 3>> M);

/**
 * @brief Wrapper of getAncestry and toNewick
 */
std::string toNewick(const std::vector<int> &v);

/**
 * @brief remove parent nodes from a Newick string
 * @todo remove not only ints but any string
 */
void removeParentAnnotations(std::string &newick);

/**
 * @brief Calculate the number of leaves in a tree from its Newick
 */
int getNumLeavesFromNewick(const std::string &newick);

// Copyright Contributors to the Pystring project.
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/imageworks/pystring/blob/master/LICENSE
/**
 * @brief Split the string around first occurrence of sep
 * @param result a 3-tuple containing the part before the separator, the separator itself, and the
 * part after the separator
 * If the separator is not found, return a 3-tuple containing the string itself, followed by two
 * empty strings
 */
void partition(const std::string &str, const std::string &sep, std::vector<std::string> &result);
inline std::vector<std::string> partition(const std::string &str, const std::string &sep) {
    std::vector<std::string> result;
    partition(str, sep, result);
    return result;
}

// Copyright Contributors to the Pystring project.
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/imageworks/pystring/blob/master/LICENSE
/**
 * @brief Split the string at the last occurrence of sep
 * @param result a 3-tuple containing the part before the separator, the separator itself, and the
 * part after the separator If the separator is not found, return a 3-tuple containing two empty
 * strings, followed by the string itself
 */
void rpartition(const std::string &str, const std::string &sep, std::vector<std::string> &result);
inline std::vector<std::string> rpartition(const std::string &str, const std::string &sep) {
    std::vector<std::string> result;
    rpartition(str, sep, result);
    return result;
}

/**
 * @brief find the most inner left leaf
 *
 * @param nw Newick representation of a tree
 * @param labels @todo
 * @param processed Whether each leaf node has been processed or not
 * @param num_leaves Number of leaves
 * @return std::pair<int, int> the left leaf + at which iteration it was found
 */
std::pair<int, int> findLeftLeaf(std::string nw, const std::vector<int> &labels,
                                 const std::vector<bool> &processed, int num_leaves);

void updateVmin(std::vector<int> &vmin, int right_leaf, int num_leaves,
                const std::vector<bool> &processed);

/**
 * @brief Update the newick string by fusing the left leaf and the right leaf
 *
 * @param nw  Newick representation of a tree
 * @param left_leaf_ind @todo
 * @param left_leaf @todo
 * @param right_leaf @todo
 * @param labels @todo
 */
void updateNewick(std::string &nw, int left_leaf_ind, int left_leaf, int right_leaf,
                  const std::vector<int> &labels);

/**
 * @brief Convert a newick-format tree to its v representation
 *
 * @param newick Newick representation of a tree
 * @param num_leaves Number of leaves
 * Saves some computation time if fed in advance
 * @return std::vector<int> Phylo2Vec representation of nw
 */
std::vector<int> newick2v(std::string newick, int num_leaves = -1);

#endif  // PHYLO2VEC_HPP