//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   ds++/test/cxx11example_move_sematics.cc
 * \author Tim M. Kelley <tkelley@lanl.gov>, Kelly G. Thompson <kgt@lanl.gov>
 * \date   Wednesday, May 24, 2017, 11:06 am
 * \brief  Demonstrate proper and improper use of std::move
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * \sa http://blog.smartbear.com/c-plus-plus/c11-tutorial-introducing-the-move-constructor-and-the-move-assignment-operator/
 */
//----------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

// forward declaration
void report_memory_locations(std::vector<double> const &v,
                             std::string const &name);

// Classes used by the example

//============================================================================//
/*!
 * \class A
 * \brief Improper use of move semantics in constructor.
 */
//============================================================================//
struct A {

public:
  /*!
   * \brief Constructor
   * \param[in] v_in v_in is lvalue; its type is rval ref
   */
  explicit A(std::vector<double> &&v_in)
      : v_(v_in) // regular vector copy ctor called b/c lvalue
  {
    /* empty */
  }

  // member data
  std::vector<double> v_;
};

//============================================================================//
/*!
 * \class B
 * \brief Proper use of move semantics in constructor.
 */
//============================================================================//
struct B {

  /*!
   * \brief Constructor
   * \param[in] v_in v_in is lvalue; its type is rval ref
   */
  explicit B(std::vector<double> &&v_in)
      : v_(std::move(v_in)) /* move casts to rval, move ctor called */
  {
    /* empty */
  }

  // member data
  std::vector<double> v_;
};

//----------------------------------------------------------------------------//
/*! \breif Demonstration of move semantics
 *
 * 1. Create a vector
 * 2. Attempt to construct a class (A), demonstrate that ownership is not
 *    transferred.
 * 3. Attempt to construct another class (B), demonstrate that ownership is
 *    transferred.
 */
//----------------------------------------------------------------------------//
void move_semantics_example(rtt_dsxx::UnitTest &ut) {
  using namespace std;

  cout << ">>> Begin demonstration...\n" << endl;

  // Create a complex object.
  cout << "Create a vector v1.\n";
  vector<double> v1 = {1, 2, 3};
  report_memory_locations(v1, "v1");

  auto const v1_loc = &v1;
  auto v1_data_loc = &v1[0];

  // Case 1:
  // Create an object, attempt to transfer ownership from v1 to a.
  // This will fail.
  cout << "\nCreate an instantiation of A that owns a copy of v1.";
  A a(move(v1));
  cout << "\nAfter call to A::ctor\n";
  report_memory_locations(v1, "v1");
  report_memory_locations(a.v_, "a.v_");

  // v1 remains unchanged! (not the behavior we want)
  if (v1_loc != &v1)
    ITFAILS;
  if (v1_data_loc != &v1[0])
    ITFAILS;
  else
    PASSMSG("v1 remains unchanged! (but we want v1 to be empty).");
  if (v1_loc == &a.v_)
    ITFAILS;
  if (v1_data_loc == &a.v_[0])
    ITFAILS;
  else
    PASSMSG("Object 'a' has made a copy of v1 (ownership not transfered).");
  if (rtt_dsxx::soft_equiv(v1.begin(), v1.end(), a.v_.begin(), a.v_.end()))
    PASSMSG("a.v_ matches v1.");
  else
    FAILMSG("A's constructor did not copy the vector's data correctly.");

  // change the data in the vector. Print the new state.
  cout << "\nExamine the behavior of 'swap'.\n";
  a.v_[0] = 4;
  a.v_[1] = 5;
  a.v_[2] = 6;
  a.v_.push_back(7); // may force the vector to resize (&a.v_[0] will change!)
  a.v_.swap(v1);
  cout << "After call to a.v_.swap(v1):\n";
  report_memory_locations(v1, "v1");
  report_memory_locations(a.v_, "a.v_");

  // a.v_ should be the memory location of v1's original data.
  if (v1_data_loc != &a.v_[0])
    ITFAILS;

  // The v1 vector's data store may have changed location when we used
  // 'push_back' above.
  v1_data_loc = &v1[0];

  // location of v1 remains unchanged (even though it's underlying data
  // container may be new).
  if (v1_loc != &v1)
    ITFAILS;

  // Case 2:
  // Create an object, attempt to transfer ownership from v1 to b.
  // This works.
  cout << "\nCreate an instantiation of B that takes ownership of v1's data.";
  B b(move(v1));
  cout << "\nAfter call to B::ctor\n";
  report_memory_locations(v1, "v1");
  report_memory_locations(b.v_, "b.v_");

  // v1 remains unchanged! (not the behavior we want)
  if (v1_loc != &v1)
    ITFAILS;
  if (v1.size() != 0)
    ITFAILS;
  else
    PASSMSG("v1 no longer has a data store (transfered to B)");
  if (v1_loc == &b.v_)
    ITFAILS;
  // v1_data_loc was set before B was constructed
  if (v1_data_loc == &b.v_[0] && b.v_.size() > 0)
    PASSMSG("Object 'b' has member data taken from v1.");
  else
    ITFAILS;
  if (rtt_dsxx::soft_equiv(v1.begin(), v1.end(), b.v_.begin(), b.v_.end()))
    FAILMSG("B's constructor did not invalidate the vector's data correctly.");
  else
    PASSMSG("'b' has taken full owndership of v1's data.");

  return;
}

//----------------------------------------------------------------------------//
void report_memory_locations(std::vector<double> const &v,
                             std::string const &name) {
  using namespace std;
  cout << name << " @ " << &v << ", " << name << " data @ ";
  if (v.size() > 0)
    cout << &v[0] << endl;
  else
    cout << "nullptr" << endl;
  cout << name << " = {";
  if (!v.empty()) {
    copy(v.begin(), prev(v.end()), ostream_iterator<double>(cout, ","));
    cout << v.back();
  }
  cout << "}" << endl;
  return;
}

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    move_semantics_example(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// end of cxx11example_move_semantics.cc
//----------------------------------------------------------------------------//
