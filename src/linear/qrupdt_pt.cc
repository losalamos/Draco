//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/qrupdt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of qrupdt
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>

#include "qrupdt.i.hh"

namespace rtt_linear
{
using std::vector;

//---------------------------------------------------------------------------//
// T = RandomContainer = vector<double>
//---------------------------------------------------------------------------//

template void qrupdt(vector<double> &r, 
		     vector<double> &qt, 
		     const unsigned n,
		     vector<double> &u, 
		     vector<double> &v);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
//                 end of qrupdt.cc
//---------------------------------------------------------------------------//
