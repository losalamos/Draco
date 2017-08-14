//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/gaussj.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of gaussj
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "gaussj.i.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

template DLL_PUBLIC_linear void gaussj(vector<double> &A, unsigned n,
                                       vector<double> &b, unsigned m);

template DLL_PUBLIC_linear void gaussj(vector<vector<double>> &A,
                                       vector<double> &b);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
// end of gaussj.cc
//---------------------------------------------------------------------------//
