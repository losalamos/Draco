//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSortPermutation.cc
 * \author Randy M. Roberts
 * \date   Mon Feb 14 14:20:45 2000
 * \note   Copyright (c) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/SortPermutation.hh"
#include <iostream>
#include <list>

namespace {

using std::cout;
using std::endl;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void printStatus(const std::string &name, bool passed) {
  // Print the status of the test.

  std::string stars;
  for (size_t i = 0; i < name.length(); i++)
    stars += '*';

  cout << "\n********" << stars << "********************\n";
  if (passed)
    cout << "**** " << name << " Self Test: PASSED ****\n";
  else
    cout << "**** " << name << " Self Test: FAILED ****\n";
  cout << "********" << stars << "********************\n";
  cout << endl;
}

//---------------------------------------------------------------------------//
template <typename IT>
inline bool testit(const std::string & /*name*/, IT first, IT last) {
  rtt_dsxx::SortPermutation lfsp(first, last);

  typedef typename std::iterator_traits<IT>::value_type value_type;
  std::vector<value_type> vv1(first, last);
  std::vector<value_type> vv2;

  for (unsigned i = 0; i < vv1.size(); i++)
    vv2.push_back(vv1[lfsp[i]]);

  IT lfi = first;
  for (unsigned i = 0; lfi != last; i++, ++lfi)
    vv1[lfsp.inv(i)] = *lfi;

  std::copy(first, last, std::ostream_iterator<value_type>(std::cout, " "));
  std::cout << std::endl;

  std::copy(lfsp.begin(), lfsp.end(),
            std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;

  std::copy(vv2.begin(), vv2.end(),
            std::ostream_iterator<value_type>(std::cout, " "));
  std::cout << std::endl;

  std::copy(lfsp.inv_begin(), lfsp.inv_end(),
            std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;

  using rtt_dsxx::isSorted;
  bool passed =
      isSorted(vv2.begin(), vv2.end()) && isSorted(vv1.begin(), vv1.end());

  return passed;
}

template <typename IT, typename CMP>
inline bool testit(const std::string & /*name*/, IT first, IT last,
                   const CMP &comp) {
  rtt_dsxx::SortPermutation lfsp(first, last, comp);

  typedef typename std::iterator_traits<IT>::value_type value_type;
  std::vector<value_type> vv1(first, last);
  std::vector<value_type> vv2;

  for (unsigned i = 0; i < vv1.size(); i++)
    vv2.push_back(vv1[lfsp[i]]);

  IT lfi = first;
  for (unsigned i = 0; lfi != last; i++, ++lfi)
    vv1[lfsp.inv(i)] = *lfi;

  std::copy(first, last, std::ostream_iterator<value_type>(std::cout, " "));
  std::cout << std::endl;

  std::copy(lfsp.begin(), lfsp.end(),
            std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;

  std::copy(vv2.begin(), vv2.end(),
            std::ostream_iterator<value_type>(std::cout, " "));
  std::cout << std::endl;

  using rtt_dsxx::isSorted;
  bool passed = isSorted(vv2.begin(), vv2.end(), comp) &&
                isSorted(vv1.begin(), vv1.end(), comp);

  return passed;
}

} // end unnamed namespace

struct Foo {
  double d;
  Foo(double d_) : d(d_) { /* empty */
  }
  bool operator<(const Foo &rhs) const { return d < rhs.d; }
  friend std::ostream &operator<<(std::ostream &os, const Foo &f) {
    return os << f.d;
  }
};

struct FooGT {
  double d;
  FooGT(double d_) : d(d_) { /* empty */
  }
  bool operator>(const FooGT &rhs) const { return d > rhs.d; }
  friend std::ostream &operator<<(std::ostream &os, const FooGT &f) {
    return os << f.d;
  }
};

template <typename F> struct evenIsLess {
  bool operator()(const F &f1, const F &f2) const {
    int i1 = static_cast<int>(f1.d);
    int i2 = static_cast<int>(f2.d);

    return i1 % 2 == 0 ? ((i2 % 2 == 0 ? i1 < i2 : true))
                       : (i2 % 2 == 0 ? false : i1 < i2);
  }
};

//---------------------------------------------------------------------------//
int main(int /*argc*/, char * /*argv*/ []) {
  cout << "Initiating test of the SortPermutation.\n";

  std::string name = "SortPermutation";

  try {
    bool passed = false;

    name = "SortPermutation(empty vector<Foo>)";
    std::vector<Foo> evf;
    passed = testit("empty vector<Foo>", evf.begin(), evf.end());
    printStatus(name, passed);

    Foo caf[] = {64, 89, 64, 73, 14, 90, 63, 14};
    const int ncaf = sizeof(caf) / sizeof(*caf);

    name = "SortPermutation(const list<Foo>)";
    const std::list<Foo> lf(caf, caf + ncaf);
    passed = testit("const list<Foo>", lf.begin(), lf.end());
    printStatus(name, passed);

    name = "SortPermutation(vector<Foo>)";
    std::vector<Foo> vf(caf, caf + ncaf);
    passed = testit("vector<Foo>", vf.begin(), vf.end());
    printStatus(name, passed);

    name = "SortPermutation(C-Array<Foo>)";
    passed = testit("C-Array<Foo>", caf, caf + ncaf);
    printStatus(name, passed);

    name = "SortPermutation(const list<Foo>, evenIsLess<Foo>)";
    const std::list<Foo> lfeven(caf, caf + ncaf);
    passed = testit("const list<Foo>", lfeven.begin(), lfeven.end(),
                    evenIsLess<Foo>());
    printStatus(name, passed);

    const FooGT cafg[] = {64, 89, 64, 73, 14, 90, 63, 14};
    const int ncafg = sizeof(cafg) / sizeof(*cafg);

    name = "SortPermutation(list<FooGT>, greater<FooGT>)";
    std::list<FooGT> lfg(cafg, cafg + ncafg);
    passed =
        testit("list<FooGT>", lfg.begin(), lfg.end(), std::greater<FooGT>());
    printStatus(name, passed);

    name = "SortPermutation(const vector<FooGT>, greater<FooGT>)";
    const std::vector<FooGT> vfg(cafg, cafg + ncafg);
    passed = testit("const vector<FooGT>", vfg.begin(), vfg.end(),
                    std::greater<FooGT>());
    printStatus(name, passed);

    name = "SortPermutation(const C-Array<FooGT>, greater<FooGT>)";
    passed = testit("const C-Array<FooGT>", cafg, cafg + ncaf,
                    std::greater<FooGT>());
    printStatus(name, passed);

  } catch (rtt_dsxx::assertion &a) {
    cout << "Failed assertion: " << a.what() << endl;
    printStatus(name, false);
    return 1;
  } catch (...) {
    cout << "tstSortPermulation: Caught unknown exception." << endl;
    printStatus(name, false);
    return 1;
  }
  cout << "Done testing SortPermutation container.\n";
  return 0;
}

//---------------------------------------------------------------------------//
// end of tstSortPermutation.cc
//---------------------------------------------------------------------------//
