#include "helper_hash.h"

#include <vector>

// boost::hash_combine
// https://stackoverflow.com/a/27952689/3446913
std::size_t hash_combine(std::size_t lhs, std::size_t rhs) {
  lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  return lhs;
}
