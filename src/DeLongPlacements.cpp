#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include "spEDMDataStruct.h"

DeLongPlacementsRes CppDeLongPlacements(const std::vector<double>& cases,
                                        const std::vector<double>& controls,
                                        const std::string& direction) {
  // Initialize variables
  size_t m = cases.size();
  size_t n = controls.size();
  size_t L = m + n;

  // Create working copies to handle direction
  std::vector<double> proc_cases = cases;
  std::vector<double> proc_controls = controls;

  // Handle direction parameter
  if (direction == ">") {
    for (auto& val : proc_cases) val = -val;
    for (auto& val : proc_controls) val = -val;
  }

  // Create combined vector with indices and class labels
  std::vector<std::pair<int, double>> Z;
  std::vector<bool> labels(L);

  // Populate case data
  for (size_t i = 0; i < m; ++i) {
    Z.emplace_back(i, proc_cases[i]);
    labels[i] = true;
  }

  // Populate control data
  for (size_t j = 0; j < n; ++j) {
    Z.emplace_back(m + j, proc_controls[j]);
    labels[m + j] = false;
  }

  // Sort combined data by value
  std::sort(Z.begin(), Z.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.second < b.second;
  });

  // Calculate XY placements
  std::vector<double> XY(L, 0.0);
  size_t current_m = 0, current_n = 0;
  size_t i = 0;

  while (i < L) {
    std::vector<int> case_indices, control_indices;
    size_t case_count = 0, control_count = 0;

    // Handle tied values
    while (true) {
      int index = Z[i].first;
      if (labels[index]) {
        ++case_count;
        case_indices.push_back(index);
      } else {
        ++control_count;
        control_indices.push_back(index);
      }

      if (i == L-1 || Z[i].second != Z[i+1].second) break;
      ++i;
    }

    // Update XY values for cases
    for (int idx : case_indices) {
      XY[idx] = current_n + control_count/2.0;
    }

    // Update XY values for controls
    for (int idx : control_indices) {
      XY[idx] = current_m + case_count/2.0;
    }

    // Accumulate counts
    current_m += case_count;
    current_n += control_count;
    ++i;
  }

  // Calculate final X and Y vectors
  DeLongPlacementsRes ret;
  double sum = 0.0;
  const double norm_n = static_cast<double>(n);
  const double norm_m = static_cast<double>(m);

  for (size_t k = 0; k < L; ++k) {
    if (labels[k]) {
      sum += XY[k];
      ret.X.push_back(XY[k] / norm_n);
    } else {
      ret.Y.push_back(1.0 - XY[k] / norm_m);
    }
  }

  // Calculate theta (AUC estimate)
  ret.theta = sum / (norm_m * norm_n);

  return ret;
}
