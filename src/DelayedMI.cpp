// #include <vector>
// #include <cmath>
// #include <limits>
// #include <numeric>
// #include <algorithm>
// #include <stdexcept>
// #include "NumericUtils.h"
// #include "CppStats.h"
// #include "CppDistances.h"
// #include <RcppThread.h>

// /***********************************************************
//  * Mutual Information (KSG estimator)
//  ***********************************************************/
// double KSGMI(
//     const std::vector<std::vector<double>>& d_xy,
//     const std::vector<double>& d_x,
//     const std::vector<double>& d_y,
//     size_t k = 3,
//     size_t alg = 0,
//     double base = 2.0,
//     bool normalize = false
// ){
//     const size_t n = d_xy.size();
//     const size_t p = d_xy[0].size();

//     double sum = 0.0;
//     double avg_log_eps = 0.0;

//     for (size_t i = 0; i < n; ++i) {   
//         auto& row = d_xy[i];

//         row.erase(
//             std::remove_if(
//                 row.begin(),
//                 row.end(),
//                 [](double v){ return std::isnan(v); }),
//             row.end());

//         if (row.size() < k)
//             throw std::runtime_error("k larger than valid neighbor count");

//         std::nth_element(row.begin(), row.begin()+k-1, row.end());
            
//         double eps = row[k-1];
//         // double eps = std::max(row[k-1], 1e-15);

//         avg_log_eps += (doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
//                         ? 0.0 : std::log(eps * 2.0);

//         size_t nx = 0, ny = 0;

//         for (size_t j = 0; j < p; ++j) {
//             if (alg == 0) {
//                 if (!std::isnan(d_x[i][j]) && d_x[i][j] < eps) nx++;
//                 if (!std::isnan(d_y[i][j]) && d_y[i][j] < eps) ny++;
//             } else {
//                 if (!std::isnan(d_x[i][j]) && d_x[i][j] <= eps) nx++;
//                 if (!std::isnan(d_y[i][j]) && d_y[i][j] <= eps) ny++;
//             } 
//         }

//         if (alg == 0)
//             sum += CppDigamma(nx+1)
//                 + CppDigamma(ny+1);
//         else
//             sum += CppDigamma(nx)
//                 + CppDigamma(ny);
//     }

//     avg_log_eps /= n;

//     double mival = CppDigamma(k)
//                  + CppDigamma(n)
//                  - sum / n;

//     if (alg == 1) mival -= 1.0 / k;

//     mival = std::max(0.0, mival);

//     if (!normalize) {
//         if (!doubleNearlyEqual(base,std::exp(1.0)))
//                 mival /= std::log(base);

//         return mival;
//     } 

//     double hxy = CppDigamma(n)
//                - CppDigamma(k)
//                + 2 * avg_log_eps;
//     if (alg == 1) hxy += 1.0 / k;

//     if (hxy <= 0) {
//         if (doubleNearlyEqual(base,std::exp(1.0)))
//             mival /= std::log(base);

//         return mival;
//     } 

//     return mival / hxy;
// }

// std::vector<double> CppDMI(const std::vector<std::vector<double>>& embedding,
//                            const std::vector<double>& vec,
//                            const std::vector<size_t>& lib,
//                            const std::vector<size_t>& pred,
//                            size_t k = 3,
//                            size_t alg = 0,
//                            double base = 2.0,
//                            bool normalize = false){
//     // Configure threads (cap at hardware concurrency)
//     size_t threads_sizet = static_cast<size_t>(std::abs(threads));
//     threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

//     const size_t n_row = pred.size();
//     const size_t n_col = lib.size();

//     // Compute pairwise distances for original vector
//     std::vector<std::vector<double>> Dx(
//         n_row, std::vector<double>(n_col, std::numeric_limits<double>::quiet_NaN()));
    
//     for (size_t i = 0; i < lib.size(); ++i) {
//         for (size_t j = 0; j < pred.size(); ++j) {
//             size_t li = lib[i];
//             size_t pi = pred[j];
//             if (pi == li) continue;
//             double dist = std::abs(vec[pi] - vec[li]);
//             if (!std::isnan(dist)) {
//                 Dx[j][i] = dist;  // assign distance; no mirroring required
//             }
//         }
//     }

//     const size_t tau_dim = embedding[0].size();
//     std::vector<double> results(tau_dim, std::numeric_limits<double>::quiet_NaN());

//     if (embedding.empty()) {
//         return results; 
//     }

//     if (parallel_level == 0){
//         for (size_t tau = 0; tau < tau_dim; ++tau) {
//             std::vector<std::vector<double>> Dy(
//                 n_row, 
//                 std::vector<double>(n_col, std::numeric_limits<double>::quiet_NaN()));
    
//             for (size_t i = 0; i < lib.size(); ++i) {
//                 for (size_t j = 0; j < pred.size(); ++j) {
//                     size_t li = lib[i];
//                     size_t pi = pred[j];
//                     if (pi == li) continue;
//                     double dist = std::abs(embedding[pi][tau] - embedding[li][tau]);
//                     if (!std::isnan(dist)) {
//                         Dx[j][i] = dist;  // assign distance; no mirroring required
//                     }
//                 }
//             }
//         }
//     } else {
//         // Parallel computation
//         RcppThread::parallelFor(1, max_E2, [&](size_t E1) {
//         size_t E2 = E1 + 1;
//         double fnn_ratio = CppSingleFNN(embedding, lib, pred, E1, E2, threads_sizet,
//                                         parallel_level, Rtol[E1 - 1], Atol[E1 - 1], L1norm);
//         results[E1 - 1] = fnn_ratio;
//         }, threads_sizet);
//     }

//     return results;
// }
