#ifndef spEDMDataStruct_H
#define spEDMDataStruct_H

#include <vector>

struct PartialCorRes {
  int first;
  double second;
  double third;

  // Default constructor
  PartialCorRes() : first(0), second(0.0), third(0.0) {}

  // Constructor to initialize all members
  PartialCorRes(int f, double s, double t) : first(f), second(s), third(t) {}
};

struct CMCRes {
  std::vector<double> causal_strength;
  std::vector<std::vector<double>> cross_mapping;

  CMCRes() = default;

  CMCRes(const std::vector<double>& causal_strength,
         const std::vector<std::vector<double>>& cross_mapping)
    : causal_strength(causal_strength), cross_mapping(cross_mapping) {}
};


struct IntersectionRes {
  size_t libsize;
  std::vector<double> Intersection;

  // Default constructor
  IntersectionRes() : libsize(0), Intersection() {}

  // Parameterized constructor
  IntersectionRes(size_t t, const std::vector<double>& x)
    : libsize(t), Intersection(x) {}
};

struct DeLongPlacementsRes {
  double theta;
  std::vector<double> X;
  std::vector<double> Y;

  // Default constructor
  DeLongPlacementsRes() : theta(0.0), X(), Y() {}

  // Parameterized constructor
  DeLongPlacementsRes(double t, const std::vector<double>& x, const std::vector<double>& y)
    : theta(t), X(x), Y(y) {}
};

#endif // spEDMDataStruct_H
