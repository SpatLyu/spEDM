#ifndef CppGridUtils_H
#define CppGridUtils_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

std::vector<std::vector<double>> CppLaggedVar4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum
);

std::vector<std::vector<std::vector<double>>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat,
    int E
);

std::vector<std::vector<double>> GenGridEmbeddings2(
    const std::vector<std::vector<double>>& mat,
    int E
);

#endif // CppGridUtils_H
