#ifndef CppLatticeUtils_H
#define CppLatticeUtils_H

#include <iostream>
#include <vector>
#include <numeric>   // for std::accumulate
#include <algorithm> // for std::sort, std::unique, std::accumulate
#include <unordered_set> // for std::unordered_set
#include <limits> // for std::numeric_limits
#include <cmath> // For std::isnan

std::vector<std::vector<int>> CppLaggedVar4Lattice(std::vector<std::vector<int>> spNeighbor,
                                                   int lagNum);

std::vector<std::vector<double>> GenLatticeEmbeddings(const std::vector<double>& vec,
                                                      const std::vector<std::vector<int>>& nb,
                                                      int E);

#endif // CppLatticeUtils_H
