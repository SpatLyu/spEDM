#include <RcppThread.h>
// [[Rcpp::depends(RcppThread)]]

// [[Rcpp::export]]
unsigned int DetectMaxNumThreads(){
  unsigned int max_threads = std::thread::hardware_concurrency();
  return max_threads;
}
