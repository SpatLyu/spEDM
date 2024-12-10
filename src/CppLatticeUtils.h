#ifndef CppLatticeUtils_H
#define CppLatticeUtils_H

#include <iostream>
#include <vector>
#include <algorithm> // for std::sort and std::unique
#include <numeric>   // for std::accumulate
#include <Rcpp.h>

std::vector<std::vector<int>> nb2vec(Rcpp::List nb);

std::vector<std::vector<int>> CppLaggedVar4Lattice(std::vector<std::vector<int>> spNeighbor,
                                                   int lagNum);

std::vector<std::vector<double>> GenEmbeddings(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nb,
                                               int E);

#endif // CppLatticeUtils_H
