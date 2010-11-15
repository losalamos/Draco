//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/tqli.cc
 * \author Kent Budge
 * \date   Thu Sep  2 15:00:32 2004
 * \brief  Specializations of tqli
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>
#include "tqli.i.hh"

namespace rtt_linear
{
using std::vector;

//---------------------------------------------------------------------------//
// T=vector<double>
//---------------------------------------------------------------------------//

template void tqli(vector<double> &d,
                   vector<double> &e,
                   const unsigned n,
                   vector<double> &z);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
//                 end of tqli.cc
//---------------------------------------------------------------------------//
