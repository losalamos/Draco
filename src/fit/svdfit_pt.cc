//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fit/svdfit_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of svdfit
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "svdfit.i.hh"
#include <vector>

namespace rtt_fit
{
using std::vector;

template DLL_PUBLIC_fit
void svdfit(vector<double> const &x,
            vector<double> const &y,
            vector<double> const &sig,
            vector<double> &a,
            vector<double> &u,
            vector<double> &v,
            vector<double> &w,
            double &chisq,
            void (&funcs)(double, vector<double> const &),
            double TOL);

} // end namespace rtt_fit

//---------------------------------------------------------------------------//
// end of svdfit_pt.cc
//---------------------------------------------------------------------------//
