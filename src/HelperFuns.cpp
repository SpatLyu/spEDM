#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits

// Function to check if a value is NA
bool isNA(double value) {
  return std::isnan(value);
}

// Function to check if a indice of int type is NA
bool checkIntNA(int value) {
  return value == std::numeric_limits<int>::min();
}

// Function to calculate the mean of a vector, ignoring NA values
double CppMean(const std::vector<double>& vec, bool NA_rm = false) {
  double sum = 0.0;
  size_t count = 0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      sum += value;
      ++count;
    }
  }
  return count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate the sum of a vector, ignoring NA values if NA_rm is true
double CppSum(const std::vector<double>& vec,
              bool NA_rm = false) {
  double sum = 0.0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      sum += value;
    }
  }
  return sum;
}

// Function to calculate the absolute difference between two vectors
std::vector<double> CppAbs(const std::vector<double>& vec1,
                           const std::vector<double>& vec2) {
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  std::vector<double> result(vec1.size());
  for (size_t i = 0; i < vec1.size(); ++i) {
    result[i] = std::abs(vec1[i] - vec2[i]);
  }
  return result;
}

// Function to normalize a vector by dividing each element by the sum of all elements
std::vector<double> CppSumNormalize(const std::vector<double>& vec,
                                    bool NA_rm = false) {
  double sum = CppSum(vec, NA_rm);
  if (sum == 0.0) {
    throw std::invalid_argument("Sum of vector elements is zero, cannot normalize.");
  }

  std::vector<double> normalizedVec(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    if (!isNA(vec[i])) {
      normalizedVec[i] = vec[i] / sum;
    } else {
      normalizedVec[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  return normalizedVec;
}
