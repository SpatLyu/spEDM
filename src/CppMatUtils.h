#ifndef CppMatUtils_H
#define CppMatUtils_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

std::vector<std::vector<double>> CppLaggedVar4Mat(
    std::vector<std::vector<double>> mat,
    int lagNum);

std::vector<std::vector<double>> GenEmbeddings(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nb,
                                               int E);

#endif // CppMatUtils_H
