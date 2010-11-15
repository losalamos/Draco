//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fit/svdfit.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of svdfit
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>

#include "svdfit.i.hh"

namespace rtt_fit
{
using std::vector;

template
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
//                 end of svdfit.cc
//---------------------------------------------------------------------------//
