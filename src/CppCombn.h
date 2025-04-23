#ifndef CPP_COMBN_H
#define CPP_COMBN_H

#include <vector>
#include <string>
#include <functional>

/**
 * @brief Generate all combinations of m elements from a given vector vec.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The input vector to generate combinations from.
 * @param m The number of elements in each combination.
 * @return std::vector<std::vector<T>> A vector containing all combinations.
 */
template <typename T>
std::vector<std::vector<T>> CppCombn(const std::vector<T>& vec, int m) {
  std::vector<std::vector<T>> result;
  std::vector<T> current;

  std::function<void(int)> combnHelper = [&](int start) {
    if (static_cast<int>(current.size()) == m) {
      result.push_back(current);
      return;
    }
    for (int i = start; i <= static_cast<int>(vec.size()) - (m - current.size()); ++i) {
      current.push_back(vec[i]);
      combnHelper(i + 1);
      current.pop_back();
    }
  };

  combnHelper(0);
  return result;
}

#endif // CPP_COMBN_H
