#ifndef __radphys_test_vectorops_hh__
#define __radphys_test_vectorops_hh__

#include <vector>
using std::vector;

vector<double> pow(const vector<double> &lhs, double n);
vector<double> sqrt(const vector<double> &lhs);
vector<double> log(const vector<double> &lhs);
vector<double> max(const vector<double> &lhs, const vector<double> &rhs);
vector<double> max(const vector<double> &lhs, double rhs);
vector<double> max(double lhs, const vector<double> &rhs);
vector<double> operator*(const vector<double> &lhs, const vector<double> &rhs);
vector<double> operator*(double lhs, const vector<double> &rhs);
vector<double> operator*(const vector<double> &lhs, double rhs);
vector<double> operator/(const vector<double> &lhs, const vector<double> &rhs);
vector<double> operator/(double lhs, const vector<double> &rhs);
vector<double> operator/(const vector<double> &lhs, double rhs);
vector<double> operator-(const vector<double> &lhs, const vector<double> &rhs);
vector<double> operator-(double lhs, const vector<double> &rhs);
vector<double> operator-(const vector<double> &lhs, double rhs);

#endif                          // __radphys_test_vectorops_hh__

