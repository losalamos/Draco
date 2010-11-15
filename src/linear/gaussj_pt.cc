//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/gaussj.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of gaussj
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>

#include "gaussj.i.hh"

namespace rtt_linear
{
using std::vector;

template void gaussj(vector<double> &A,
		     unsigned n,
		     vector<double> &b,
		     unsigned m);

template void gaussj(vector<vector<double> > &A,
		     vector<double> &b);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
//                 end of gaussj.cc
//---------------------------------------------------------------------------//
