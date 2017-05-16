//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstIndex_Counter.cc
 * \author Mike Buksas
 * \date   Wed Feb  1 08:58:48 2006
 * \brief  Unit test for Index_Counter
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Index_Converter.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
/*
 Issues with this test:

 tstIndex_Counter.cc will produce the following warning if inlining
 is enabled (gcc):
 
 \code 
 /.../source/src/ds++/test/../Index_Counter.hh: In member function 'int rtt_dsxx::Index_Converter<D, OFFSET>::get_next_index(const rtt_dsxx::Index_Counter<D, OFFSET>&, int) const [with unsigned int D = 3u, int OFFSET = 1]':
 /.../source/src/ds++/test/../Index_Counter.hh:69: warning: array subscript is above array bounds
 \endcode

 This appears to be an issue with speculative instructions issued by
 gcc due to the pattern of access indices found in the test.  This warning
 only appears when inlining is enabled (-02 and above) and DBC is turned on.
 Our normal test setup turns DBC off for release (-03) builds.
*/
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_index_counter(rtt_dsxx::UnitTest &ut) {
  unsigned dimensions[] = {3, 4, 5};
  Index_Converter<3, 1> box(dimensions);
  Index_Converter<3, 1>::Counter it = box.counter();

  if (it.get_index() != 1)
    ITFAILS;
  if (it.get_index(0) != 1)
    ITFAILS;
  if (it.get_index(1) != 1)
    ITFAILS;
  if (it.get_index(2) != 1)
    ITFAILS;
  if (!it.is_in_range())
    ITFAILS;
  if (it.get_indices()[0] != 1)
    ITFAILS;
  if (it.get_indices()[1] != 1)
    ITFAILS;
  if (it.get_indices()[2] != 1)
    ITFAILS;
  vector<unsigned> it_copy(3);
  it.get_indices(it_copy.begin());
  if (it_copy[0] != 1)
    ITFAILS;
  if (it_copy[1] != 1)
    ITFAILS;
  if (it_copy[2] != 1)
    ITFAILS;

  ++it;

  if (it.get_index() != 2)
    ITFAILS;
  if (it.get_index(0) != 2)
    ITFAILS;
  if (it.get_index(1) != 1)
    ITFAILS;
  if (it.get_index(2) != 1)
    ITFAILS;
  if (!it.is_in_range())
    ITFAILS;

  --it;
  --it;

  if (it.is_in_range())
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("done with test_index_counter().");
  return;
}

//---------------------------------------------------------------------------------------//
void test_looping(rtt_dsxx::UnitTest &ut) {
  unsigned dimensions[] = {3, 4, 5};
  Index_Converter<3, 1> box(dimensions);

  int index = 1;
  for (Index_Counter<3, 1> it(box); it.is_in_range(); ++it) {
    const int it_index = it.get_index();

    // Check the returned index against a manual count.
    if (it_index != index++)
      ITFAILS;

    // Check the first and last index directly.
    if ((it_index - 1) % 3 + 1 != it.get_index(0))
      ITFAILS;
    if ((it_index - 1) / 12 + 1 != it.get_index(2))
      ITFAILS;
  }
  return;
}

//---------------------------------------------------------------------------------------//
void test_next_index(rtt_dsxx::UnitTest &ut) {
  unsigned dimensions[] = {3, 4, 5};
  Index_Converter<3, 1> box(dimensions);
  Index_Counter<3, 1> it = box.counter();

  if (it.get_index() != 1)
    ITFAILS;
  if (box.get_next_index(it, 1) != -1)
    ITFAILS;
  if (box.get_next_index(it, 2) != 2)
    ITFAILS;
  if (box.get_next_index(it, 3) != -1)
    ITFAILS;
  if (box.get_next_index(it, 4) != 4)
    ITFAILS;
  if (box.get_next_index(it, 5) != -1)
    ITFAILS;
  if (box.get_next_index(it, 6) != 13)
    ITFAILS;
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_index_counter(ut);
    test_looping(ut);
    test_next_index(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstIndex_Counter.cc
//---------------------------------------------------------------------------//
