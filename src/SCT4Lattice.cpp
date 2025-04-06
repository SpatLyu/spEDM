#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

std::vector<double> SCT4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    size_t boot,
    size_t threads
){

}
