#ifndef HELPER_HASH_H
#define HELPER_HASH_H

#define MACHINE_EPS 1e-10

#include <vector>

std::size_t hash_combine(std::size_t lhs, std::size_t rhs);

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return hash_combine(x.first, x.second);
  }
};


struct vechash {
public:
  template <typename T>
  std::size_t operator()(const std::vector<T> &x) const
  {
    std::size_t h = 0;
    
    for (std::size_t i = 0; i < x.size(); ++i) {
      h = hash_combine(h, x[i]);
    }
    
    return h;
  }
};

#endif
