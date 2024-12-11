#ifndef CppGridUtils_H
#define CppGridUtils_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

std::vector<std::vector<double>> CppLaggedVar4Grid(
    std::vector<std::vector<double>> mat,
    int lagNum);

std::vector<std::vector<double>> GenGridEmbeddings(const std::vector<double>& vec,
                                                   const std::vector<std::vector<int>>& nb,
                                                   int E);

#endif // CppGridUtils_H
