//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/FM12.hh
 * \author Kent Budge
 * \date   Tue Sep 21 12:06:09 2004
 * \brief  Compute Fermi-Dirac function of 1/2 order
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_FM12_hh
#define sf_FM12_hh

namespace rtt_sf
{
//! Calculate Fermi-Dirac integral of index -1/2.
template<class OrderedField>
OrderedField FM12(OrderedField const &x);

} // end namespace rtt_sf

#endif // sf_FM12_hh

//---------------------------------------------------------------------------//
//              end of sf/FM12.hh
//---------------------------------------------------------------------------//
