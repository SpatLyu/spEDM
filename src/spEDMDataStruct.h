#ifndef spEDMDataStruct_H
#define spEDMDataStruct_H

struct PartialCorRes {
  int first;
  double second;
  double third;

  // Default constructor
  PartialCorRes() : first(0), second(0.0), third(0.0) {}

  // Constructor to initialize all members
  PartialCorRes(int f, double s, double t) : first(f), second(s), third(t) {}
};

#endif // spEDMDataStruct_H
