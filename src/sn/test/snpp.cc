//----------------------------------*-C++-*----------------------------------//
// snpp.cc
// Scott Turner
// 17 March 1998
//---------------------------------------------------------------------------//
// @> Main program for testing the structured sweeper.
//---------------------------------------------------------------------------//

// SN++ is a single node, ordered sweep, C++ version of:
//
// Method 3 Sweep - mesh is swept by a diagonal plane one angle at a time
//                  with block recursion, in data parallel.
//
// In other words, it is a major simplification of the F90 Method 3 Parallel
// Sweeper.
//

#include "sn/read_input.hh"
#include "sn/inner_iter.hh"

#include <iostream.h>

int main()
{

  // get input data from the input object

  read_input snpp_in;
  snpp_in.read_data();

  int  it, jt, kt, mm, isct;
  snpp_in.get_basic_data( it, jt, kt, mm, isct );

  // print output header to screen

  int isn;      // quadrature order, the N in SN

  if (mm==3)
    isn = 4;
  else
    isn = 6;

  cout << "Method 3 - Ordered Single Node C++ version" << endl;
  cout << " " << it << " x " << jt << " x " << kt << endl;
  cout << " S" << isn << "P" << isct << endl;

  // do the inner iterations

  inner_iter inner_iter_object;
  inner_iter_object.do_inner_iter();

  return (0);
}

//---------------------------------------------------------------------------//
//                              end of snpp.cc
//---------------------------------------------------------------------------//

