//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Release.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 09:48:39 2000
 * \brief  Release function for the quadrature library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __quadrature_Release_hh__
#define __quadrature_Release_hh__

//===========================================================================//
/*!
 * \page quadrature_overview Overview of the quadrature package
 * \version 1_0_0
 * \author  Kelly Thompson
 * 
 * The quadrature package provides services related to the angular
 * transport discretization.  This package provides many different
 * quadrature sets (Level Symmetric, Gauss Legendre, etc.) for
 * problems of varying dimension.
 *
 * This package also performs other related services including
 * integrating angularly dependent variables over the quadrature set
 * range (psi2phi).
 */
//===========================================================================//
/*!
 * \namespace rtt_quadrature
 *
 * \brief Namespace that contains the quadrature package classes and variables.
 *
 */
//===========================================================================//

#include <string>

namespace rtt_quadrature 
{
    //! Query package for the release number.
    const std::string release();
}

#endif                          // __quadrature_Release_hh__

//---------------------------------------------------------------------------//
//                           end of quadrature/Release.hh
//---------------------------------------------------------------------------//
