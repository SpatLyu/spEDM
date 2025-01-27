#ifndef CppGridUtils_H
#define CppGridUtils_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

std::vector<std::vector<double>> CppLaggedVar4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum
);

std::vector<std::vector<double>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat,
    int E,
    bool includeself);

#endif // CppGridUtils_H
