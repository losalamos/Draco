//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/L2norm_pt.cc
 * \author Kent Budge
 * \date   Tue Sep 18 08:22:08 2007
 * \brief  Preinstantiate template function L2norm
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.  
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef norms_L2norm_hh
#define norms_L2norm_hh

#include "L2norm.i.hh"
#include "ds++/config.h"
#include <vector>

namespace rtt_norms {

//! Compute the L2-norm of a vector.
template DLL_PUBLIC_norms double L2norm(std::vector<double> const &x);

} // end namespace rtt_norms

#endif // norms_L2norm_hh

//---------------------------------------------------------------------------//
// end of norms/L2norm_pt.cc
//---------------------------------------------------------------------------//
