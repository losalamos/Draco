//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   container/test/tstSlice.cc
 * \author Kent Budge
 * \date   Thu Jul  8 08:02:51 2004
 * \brief  Test the Slice subset container class.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Slice.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstSlice(UnitTest &ut) {
  vector<unsigned> v(5);
  for (unsigned i = 0; i < 5; i++)
    v[i] = i;
  Slice<vector<unsigned>::iterator> s = slice(v.begin() + 1, 2, 2);
  if (s.size() == 2)
    ut.passes("size of vector Slice is correct");
  else
    ut.failure("size of vector Slice is NOT correct");
  if (s[1] == 3)
    ut.passes("indexing of Slice is correct");
  else
    ut.failure("size of Slice is NOT correct");
  if (s.begin() < s.end())
    ut.passes("ordering of Slice is correct");
  else
    ut.failure("ordering of Slice is NOT correct");
  if (!(s.end() < s.begin()))
    ut.passes("ordering of Slice is correct");
  else
    ut.failure("ordering of Slice is NOT correct");

  Slice<vector<unsigned>::iterator>::iterator i = s.begin();
  if (*i == 1)
    ut.passes("deref of begin is correct");
  else
    ut.failure("deref of begin is NOT correct");
  if (i[0] == 1)
    ut.passes("operator[] of begin is correct");
  else
    ut.failure("operator[] of begin is NOT correct");
  if (*(i + 1) == 3)
    ut.passes("operator+ is correct");
  else
    ut.failure("operator+ is NOT correct");
  ++i;
  if (*i == 3)
    ut.passes("deref of ++begin is correct");
  else
    ut.failure("deref of ++begin is NOT correct");
  if (i - s.begin() == 1)
    ut.passes("operator- is correct");
  else
    ut.failure("operator- is NOT correct");
  i++;
  if (i != s.end())
    ut.failure("++ past end is NOT correct");
  else
    ut.passes("++ past end is correct");

  Slice<vector<unsigned>::iterator>::const_iterator ci = s.begin();
  if (ci.first() == v.begin() + 1)
    ut.passes("first is correct");
  else
    ut.failure("first is NOT correct");
  if (ci.offset() == 0)
    ut.passes("offset is correct");
  else
    ut.failure("offset is NOT correct");
  if (ci.stride() == 2)
    ut.passes("stride is correct");
  else
    ut.failure("stride is NOT correct");
  if (*ci == 1)
    ut.passes("deref of begin is correct");
  else
    ut.failure("deref of begin is NOT correct");
  if (ci[0] == 1)
    ut.passes("operator[] of begin is correct");
  else
    ut.failure("operator[] of begin is NOT correct");
  if (*(ci + 1) == 3)
    ut.passes("operator+ is correct");
  else
    ut.failure("operator+ is NOT correct");
  ++ci;
  if (*ci == 3)
    ut.passes("deref of ++begin is correct");
  else
    ut.failure("deref of ++begin is NOT correct");
  if (s.begin() < ci)
    ut.passes("correct const iterator ordering");
  else
    ut.failure("NOT correct const iterator ordering");
  if (ci - s.begin() == 1)
    ut.passes("operator- is correct");
  else
    ut.failure("operator- is NOT correct");
  ci++;
  if (ci != s.end())
    ut.failure("++ past end is NOT correct");
  else
    ut.passes("++ past end is correct");

  Slice<vector<unsigned>::iterator> s2(v.begin(), 3, 2);
  if (s2.size() == 3)
    ut.passes("size of Slice is correct");
  else
    ut.failure("size of Slice is NOT correct");
  if (s2[1] == 2)
    ut.passes("indexing of Slice is correct");
  else
    ut.failure("size of Slice is NOT correct");

  double da[6];
  Slice<double *> das(da, 2, 3);
  if (das.size() == 2)
    ut.passes("size of Slice is correct");
  else
    ut.failure("size of Slice is NOT correct");

  vector<double> db_vector(6);
  double *const db = &db_vector[0];
  db[0] = 0;
  Slice<vector<double>::iterator> dbs(db_vector.begin(), 2, 3);
  if (dbs.size() == 2)
    ut.passes("size of Slice is correct");
  else
    ut.failure("size of Slice is NOT correct");

  Slice<vector<unsigned>::iterator> const cs = s;
  if (cs[1] == 3)
    ut.passes("indexing of const Slice is correct");
  else
    ut.failure("size of const Slice is NOT correct");
  if (cs.front() == 1)
    ut.passes("front of const Slice is correct");
  else
    ut.failure("front of const Slice is NOT correct");
  if (cs.back() == 3)
    ut.passes("back of const Slice is correct");
  else
    ut.failure("back of const Slice is NOT correct");
  if (cs.begin() < cs.end())
    ut.passes("ordering of const Slice is correct");
  else
    ut.failure("ordering of const Slice is NOT correct");
  if (!(cs.end() < cs.begin()))
    ut.passes("ordering of const Slice is correct");
  else
    ut.failure("ordering of const Slice is NOT correct");
  if (!cs.empty())
    ut.passes("emptiness of const Slice is correct");
  else
    ut.failure("emptiness of const Slice is NOT correct");
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstSlice(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstSlice.cc
//---------------------------------------------------------------------------//
