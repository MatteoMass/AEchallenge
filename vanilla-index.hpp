#pragma once

#include "piecewise_linear_model.hpp"
#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)

/**
 * An implementation of a one-level PGM-index that uses binary search on the segments.
 */
class MyIndex {
  size_t n;                             ///< The number of elements this index was built on.
  std::vector<pgm::Segment> segments;   ///< The segments composing the index.
  const std::vector<uint64_t> &data;    ///< The data vector this index was built on.
  static constexpr size_t epsilon = 64; ///< The space-time trade-off parameter of the PGM-index.

  /** Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key. */
  auto segment_for_key(uint64_t key) const {
    // We use a simple binary search here
    auto it = std::upper_bound(segments.begin(), segments.end(), key);
    return it == segments.begin() ? it : std::prev(it);
  }

public:
  /** Constructs the index on the given sorted keys. */
  MyIndex(const std::vector<uint64_t> &data) : n(data.size()), segments(), data(data) {
    // We run piecewise linear model construction algorithm with the given epsilon value
    segments = pgm::make_pgm_segments(data, epsilon);

    // We add an extra sentinel segment with intercept equal to n, so that when last segment (*it) outputs a pos > n
    // the expression std::min((*it)(k), std::next(it)->intercept) in nextGEQ will be equal to n.
    segments.emplace_back(data.back() + 1, 0, n);
  }

  /** Returns the element greater than or equal to a given key. */
  uint64_t nextGEQ(uint64_t key) const { // This function is required for the challenge
    auto it = segment_for_key(key);
    auto pos = std::min<size_t>((*it)(key), std::next(it)->intercept);
    auto lo = PGM_SUB_EPS(pos, epsilon);
    auto hi = PGM_ADD_EPS(pos, epsilon, n);
    return *std::lower_bound(data.begin() + lo, data.begin() + hi, key);
  }

  /** Returns the size of the index in bytes. */
  size_t size_in_bytes() const { // This function is required for the challenge
    return sizeof(this) + segments.size() * sizeof(pgm::Segment);
  }
};