//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/rsolv_pt.cc
 * \author Kent Budge
 * \date   Tue Aug 10 13:01:02 2004
 * \brief  Specializations of rsolv
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "rsolv.hh"
#include <vector>

namespace rtt_linear {

using std::vector;

//---------------------------------------------------------------------------//
// T=vector<double> const &
//---------------------------------------------------------------------------//

template DLL_PUBLIC_linear void rsolv(const vector<double> &R, const unsigned n,
                                      vector<double> &b);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
// end of rsolv.cc
//---------------------------------------------------------------------------//
