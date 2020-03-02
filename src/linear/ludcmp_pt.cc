//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   linear/ludcmp_pt.cc
 * \author Kent Budge
 * \date   Thu Jul  1 10:54:20 2004
 * \brief  LU decomposition
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "ludcmp.i.hh"
#include "ds++/Slice.hh"

namespace rtt_linear {

template DLL_PUBLIC_linear void ludcmp(vector<double> &a,
                                       vector<unsigned> &indx, double &d);

template DLL_PUBLIC_linear void lubksb(vector<double> const &a,
                                       vector<unsigned> const &indx,
                                       vector<double> &b);

template DLL_PUBLIC_linear void
lubksb(vector<double> const &a, vector<unsigned> const &indx,
       rtt_dsxx::Slice<vector<double>::iterator> &b);

} // end namespace rtt_linear

//----------------------------------------------------------------------------//
// end of ludcmp_pt.cc
//----------------------------------------------------------------------------//
