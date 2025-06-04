#include <vector>
#include <numeric>
#include <cmath>
#include "CppLatticeUtils.h"

std::vector<std::vector<double>> SLMUni4Lattice(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    size_t step,
    double alpha,
    double escape_threshold = 1e10
){
  std::vector<int> lib(vec.size());
  // for(size_t i = 0; i < vec.size(); ++i){
  //   lib[i] = static_cast<int>(i);
  // }
  std::iota(lib.begin(), lib.end(), 0);

  std::vector<std::vector<int>> neighbors = GenLatticeNeighbors(vec, nb, lib, k);

  std::vector<std::vector<double>> res(vec.size(),
                                       std::vector<double>(step + 1,
                                                           std::numeric_limits<double>::quiet_NaN()));

  for(size_t i = 0; i < vec.size(); ++i){
    res[i][0] = vec[i];
  }

  // Simulation loop
  for (size_t s = 1; s <= step; ++s){
    for(size_t currentIndex = 0; currentIndex < vec.size(); ++currentIndex){
      double v_neighbors = 0;
      double valid_neighbors = 0;
      std::vector<int> local_neighbors = neighbors[currentIndex];
      for (size_t i = 0; i < local_neighbors.size(); ++i) {
        v_neighbors += res[local_neighbors[i]][s - 1];
        valid_neighbors += 1;
      }

      // Apply spatial logistic map equation (adjusted for varying degrees)
      double v_next = 1 - alpha * res[currentIndex][s - 1] * v_neighbors / valid_neighbors;

      // Escape check
      if (std::abs(v_next) <= escape_threshold){
        res[currentIndex][s] = v_next;
      }
    }
  }

  return res;
}
