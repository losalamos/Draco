//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/PoomaAssign_pt.cc
 * \author Randy M. Roberts
 * \date   Thu Nov 18 14:30:16 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Pooma_pt.hh"
#include "PoomaAssign_pt.hh"

#define PETE_SCALAR_VEKTOR(TYPE, SIZE) PETE_Scalar<Vektor<TYPE,SIZE> >

#define VektorAssigner(TYPE, SIZE) \
Assigner(VEKTOR(TYPE,SIZE), PETE_SCALAR_VEKTOR(TYPE, SIZE));

VektorAssigner(int, 6);
VektorAssigner(double, 6);
VektorAssigner(long, 6);
VektorAssigner(int, 8);
VektorAssigner(double, 8);
VektorAssigner(long, 8);


//---------------------------------------------------------------------------//
//                              end of PoomaAssign_pt.cc
//---------------------------------------------------------------------------//
