#ifndef CppUtils_H
#define CppUtils_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <limits>
#include "HelperFuns.h"

std::vector<std::vector<int>> CppLaggedIndices(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nbmat,
                                               int lagNum);

std::vector<std::vector<double>> GenEmbeddings(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nbmat,
                                               int E);

std::vector<std::vector<double>> CppDist(const std::vector<std::vector<double>>& matrix);

std::vector<std::vector<int>> CppClosestIndices(const std::vector<std::vector<double>>& distmat,
                                                int libsize);

std::vector<std::vector<double>> CCMWeight(const std::vector<std::vector<double>>& distmat,
                                           const std::vector<std::vector<int>>& closestIndices,
                                           int libsize);

#endif // CppUtils_H
