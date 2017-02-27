//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/norm.hh
 * \author Kent Budge
 * \date   Tue Sep 18 08:53:41 2007
 * \brief  Declare template class norm
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef norms_norm_hh
#define norms_norm_hh

namespace rtt_norms {

//---------------------------------------------------------------------------//
//! Compute the norm of a scalar.
template <class Field> double norm(Field const &x) { return x * x; }

} // end namespace rtt_norms

#endif // norms_norm_hh

//---------------------------------------------------------------------------//
// end of norms/norm.hh
//---------------------------------------------------------------------------//
