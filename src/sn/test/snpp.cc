// snpp.cc
// Scott Turner
// 19 February 1998
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

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

#include "sn/test/protos.hh"

int main()
{

#include "sn/test/snpp.hh"

  int it;       // total number of mesh cells in the x-direction
  int jt;       // total number of mesh cells in the y-direction
  int kt;       // total number of mesh cells in the z-direction
  int mm;       // number of   quadrature points (angles) per quadrant
  int nm;       // number of   flux and source angular moments
  int isct;     // legendre order of scattering
  int isctp;    // isct + 1
  int ibl;      // left boundary condition 0/1 = vacuum/reflective
  int ibb;      // bottom boundary condition 0/1 = vacuum/reflective
  int ibfr;     // front boundary condition 0/1 = vacuum/reflective
  int iprint;   // 0/1 = no/yes print fluxes
  int ifxg;     // 0/1 = no/yes set to zero flux fixup
  int isn;      // quadrature order, the N in SN

  REAL dx;      // mesh size in the x-direction
  REAL dy;      // mesh size in the y-direction
  REAL dz;      // mesh size in the z-direction
  REAL epsi;    // convergence precision or,
                //   if negative, then the number of iterations to do

// read input data

  ifstream input_file;
  input_file.open("snpp.in");

  if (input_file.bad())
  {
    cerr << "Error: Could not open snpp.in" << endl;
    exit (8);
  }

  input_file >> it   >> jt     >> kt   >> mm   >> isct;
  input_file >> dx   >> dy     >> dz   >> epsi;
  input_file >> ifxg >> iprint;
  input_file >> ibl  >> ibb    >> ibfr;

  input_file.close();

// test for input errors and set parameters based on input

  if (mm==3)
    isn = 4;                                  // Sn order
  else if (mm==6)
    isn = 6;
  else
  {
    cerr << "Error: mm must equal 3 or 6" << endl;
    exit (8);
  }

  if (!(isct==0 || isct==1))
  {
    cerr << "Error: isct must equal 0 or 1" << endl;
    exit (8);
  }
  isctp = isct + 1;

  nm = (isct+1) * (isct+1);                   // flux and source angular moments
      
// header for output

  cout << "Method 3 - Ordered Single Node C++ version" << endl;
  cout << " " << it << " x " << jt << " x " << kt << endl;
  cout << " S" << isn << "P" << isct << endl;

// call inner, the driver for the solution
      
  inner( it,     jt,    kt,     mm,    nm,
         isct,   isctp, ibl,    ibb,   ibfr,
         iprint, ifxg,  dx,    dy,     dz,
         epsi                               );

  return (0);
}

//---------------------------------------------------------------------------//
//                              end of snpp.cc
//---------------------------------------------------------------------------//

