#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Computes the directional symbolic causality strength between two spatial variables
 *        based on lattice neighborhood structure and symbolization.
 *
 * This function implements a symbolic causality test statistic for spatial data over lattice structures.
 * It follows a non-parametric entropy-based approach to evaluate whether variable x Granger-causes y,
 * as defined in the referenced paper using symbolized spatial neighbors.
 *
 * Specifically, the function calculates:
 *   - H(y): the entropy of the symbolized y values.
 *   - H(y | x): the conditional entropy of y given x.
 *   - H(x): the entropy of the symbolized x values.
 *   - H(x | y): the conditional entropy of x given y.
 * Then it computes:
 *   - sc_x_to_y = H(y) - H(y | x): symbolic causality strength from x to y.
 *   - sc_y_to_x = H(x) - H(x | y): symbolic causality strength from y to x.
 *
 * These correspond to the delta entropy test statistic described in Equation (24)
 * of the attached reference, where a significant positive value of sc_x_to_y indicates
 * that x contains additional information about y not present in y’s own spatial lags.
 *
 * @param x    Input vector of variable x, assumed to be on a lattice/grid.
 * @param y    Input vector of variable y, assumed to be on the same lattice/grid as x.
 * @param nb   Neighborhood list for each lattice point (e.g., rook or queen adjacency).
 * @param k    The number of neighbors considered in the lattice symbolization.
 * @param base The logarithm base used in entropy calculation (default is 2).
 *
 * @return A std::vector<double> containing two values:
 *         [0] sc_x_to_y — symbolic causality strength from x to y;
 *         [1] sc_y_to_x — symbolic causality strength from y to x.
 *
 */
std::vector<double> SCTSingle4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    double base = 2
){
  std::vector<double> sx = GenLatticeSymbolization(x,nb,k);
  std::vector<double> sy = GenLatticeSymbolization(y,nb,k);
  double Hx = CppEntropy_Disc(sx,base,false); // H(x)
  double Hy = CppEntropy_Disc(sy,base,false); // H(y)
  double Hxy = CppConditionalEntropy_Disc(sx, sy, base, false); // H(x | y)
  double Hyx = CppConditionalEntropy_Disc(sy, sx, base, false); // H(y | x)
  double sc_x_to_y = Hy - Hyx;
  double sc_y_to_x = Hx - Hxy;
  return {sc_x_to_y,sc_y_to_x};
}

// std::vector<double> SCT4Lattice(
//     const std::vector<double>& x,
//     const std::vector<double>& y,
//     const std::vector<std::vector<int>>& nb,
//     const std::vector<int>& block,
//     size_t k,
//     int boot = 399,
//     double base = 2
// ){
//
// }
