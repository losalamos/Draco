//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/F32.hh
 * \author Kent Budge
 * \date   Tue Sep 21 12:06:09 2004
 * \brief  Compute Fermi-Dirac function of 1/2 order
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_F32_hh
#define sf_F32_hh

namespace rtt_sf
{
//! Calculate Fermi-Dirac integral of index 3/2.
template<class OrderedField>
OrderedField F32(OrderedField const &x);

} // end namespace rtt_sf

#endif // sf_F32_hh

//---------------------------------------------------------------------------//
//              end of sf/F32.hh
//---------------------------------------------------------------------------//
