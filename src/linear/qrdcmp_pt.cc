//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/qrdcmp_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of qrdcmp
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#include "qrdcmp.i.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

//---------------------------------------------------------------------------//
// RandomContainer=vector<double>
//---------------------------------------------------------------------------//

template DLL_PUBLIC_linear bool qrdcmp(vector<double> &a, unsigned n,
                                       vector<double> &c, vector<double> &d);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
// end of qrdcmp_pt.cc
//---------------------------------------------------------------------------//
