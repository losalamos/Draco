//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/Pooma_pt.cc
 * \author Randy M. Roberts
 * \date   Wed Nov 17 15:48:09 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Pooma_pt.hh"

template class Cartesian<Dimension_s>;

#define FIELD(TYPE) \
template class Field<TYPE, Dimension_s, Cartesian<Dimension_s>, Vert>; \
template class Field<TYPE, Dimension_s, Cartesian<Dimension_s>, Cell>; \

FIELD(int);
FIELD(long);
FIELD(double);
FIELD(VEKTOR(int,6));
FIELD(VEKTOR(long,6));
FIELD(VEKTOR(double,6));
FIELD(VEKTOR(int,8));
FIELD(VEKTOR(long,8));
FIELD(VEKTOR(double,8));
    
template class CompressedBrickIterator<int, Dimension_s>;
template class CompressedBrickIterator<long, Dimension_s>;
template class CompressedBrickIterator<double, Dimension_s>;

#define VektorCBI(TYPE, SIZE) \
template class CompressedBrickIterator<Vektor<TYPE,SIZE>, Dimension_s>;

VektorCBI(int,8);
VektorCBI(long,8);
VektorCBI(double,8);
VektorCBI(int,6);
VektorCBI(long,6);
VektorCBI(double,6);


//---------------------------------------------------------------------------//
//                              end of Pooma_pt.cc
//---------------------------------------------------------------------------//
