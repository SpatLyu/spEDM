#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "NumericUtils.h"
#include "CppStats.h"
#include <RcppThread.h>

/***********************************************************
 * Mutual Information (KSG estimator)
 ***********************************************************/
double KSGMI(
    const std::vector<std::vector<double>>& d_xy,
    const std::vector<double>& d_x,
    const std::vector<double>& d_y,
    size_t k = 3,
    size_t alg = 0,
    double base = 2.0)
{
    const size_t n = d_xy.size();

    double sum = 0.0;
    double avg_log_eps = 0.0;

    for (size_t i = 0; i < n; ++i) {   
        auto& row = d_xy[i];
        row[i] = std::numeric_limits<double>::quiet_NaN();

        row.erase(
            std::remove_if(
                row.begin(),
                row.end(),
                [](double v){ return std::isnan(v); }),
            row.end());

        if (row.size() < k)
            throw std::runtime_error("k larger than valid neighbor count");

        std::nth_element(row.begin(),row.begin()+k-1,row.end());
            
        double eps = row[k-1];
        // double eps = std::max(row[k-1], 1e-15);

        avg_log_eps += (doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                        ? 0.0 : std::log(eps * 2.0);

        size_t nx = 0, ny = 0;

        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;

            if (alg == 0) {
                if (!std::isnan(d_x[i][j]) && d_x[i][j] < eps) nx++;
                if (!std::isnan(d_y[i][j]) && d_y[i][j] < eps) ny++;
            } else {
                if (!std::isnan(d_x[i][j]) && d_x[i][j] <= eps) nx++;
                if (!std::isnan(d_y[i][j]) && d_y[i][j] <= eps) ny++;
            } 
        }

        if (alg == 0)
            sum += CppDigamma(nx+1)
                + CppDigamma(ny+1);
        else
            sum += CppDigamma(nx)
                + CppDigamma(ny);
    }

    avg_log_eps /= n;

    double mival = CppDigamma(k)
                 + CppDigamma(n)
                 - sum / n;

    if (alg == 1) mival -= 1.0 / k;

    mival = std::max(0.0, mival);

    double hxy = CppDigamma(n)
               - CppDigamma(k)
               + 2 * avg_log_eps;
    if (alg == 1) hxy += 1.0 / k;

    if (hxy <= 0) {
        if (doubleNearlyEqual(base,std::exp(1.0)))
            mival /= std::log(base);

        return mival;
    } 

    return mival / hxy;
}
