#ifndef CppStats_H
#define CppStats_H

double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false);

double CppSignificance(double r, int n);

std::vector<double> CppConfidence(double r, int n,
                                  double level = 0.05);

std::vector<double> LinearTrendRM(const std::vector<double>& vec,
                                  const std::vector<double>& xcoord,
                                  const std::vector<double>& ycoord,
                                  bool NA_rm = false);

#endif // CppStats_H
