#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <utility>
#include <limits>
#include "DataStruct.h"

PatternProjectionRes PatternProjectionSingle(
    const std::vector<std::vector<double>>& SMy,
    const std::vector<std::vector<double>>& Dx,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    int num_neighbors = 4,
    int zero_tolerance = std::numeric_limits<int>::min()
) {
  const size_t& n_row = SMy.size();
  const size_t& n_col = SMy[0].size();

  if (zero_tolerance == std::numeric_limits<int>::min()) {
    zero_tolerance = static_cast<int>(n_col);
  }

  std::vector<std::vector<double>> predSMy(n_row,
                                           std::vector<double>(n_col,std::numeric_limits<double>::quiet_NaN()));

}
