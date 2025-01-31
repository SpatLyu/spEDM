#ifndef spEDMDataStruct_H
#define spEDMDataStruct_H

struct PartialCorRes {
  int first;
  double second;
  double third;

  // Constructor to initialize all members
  PartialCorRes(int f, double s, double t) : first(f), second(s), third(t) {}
};

#endif // spEDMDataStruct_H
