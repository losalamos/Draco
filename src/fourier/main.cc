//----------------------------------*-C++-*----------------------------------//
// main.cc
// John Gulick
// Tue Aug  3 13:51:27 1999
//---------------------------------------------------------------------------//
// @> Temporary main for testing Fourier Analysis Package.
//---------------------------------------------------------------------------//

#ifndef main_H
#define main_H
#include "fourier.hh"

int main ()
{
Fourier fourier;
fourier.input();
fourier.solve();
fourier.plot();
//fourier.print();
}
 

#endif                          

//---------------------------------------------------------------------------//
//                              end of fourier/main.cc
//---------------------------------------------------------------------------//
