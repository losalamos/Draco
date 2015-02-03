//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fit/svdfit.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate a generalized least squares fit.
 * \note   Copyright (C) 2006-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef utils_svdfit_hh
#define utils_svdfit_hh

namespace rtt_utils
{

//! Compute a generalized least squares fit.
template<class RandomContainer,
         class Functor>
void svdfit(RandomContainer const &x,
            RandomContainer const &y,
            RandomContainer const &sig,
            RandomContainer &a,
            RandomContainer &u,
            RandomContainer &v,
            RandomContainer &w,
            double &chisq,
            Functor &funcs,
            double TOL=1.0e-13);

} // end namespace rtt_utils

#endif // utils_svdfit_hh

//---------------------------------------------------------------------------//
// end of utils/svdfit.hh
//---------------------------------------------------------------------------//
