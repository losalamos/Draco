//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Fleck_Factors.hh
 * \author Thomas M. Evans
 * \date   Tue Feb 18 17:46:50 2003
 * \brief  Fleck_Factors class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Fleck_Factors_HH
#define RTT_imc_Fleck_Factors_HH

#include "ds++/SP.hh"

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \struct Fleck_Factors
 *
 * \brief Holds Fleck factors in frequency-independent container.
 *
 * This class allows a host of frequency-dependent opacity classes to access
 * Fleck factors, which are frequency-independent.  By placing the Fleck
 * factors in this container, we avoid multiple copies of the factors.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
struct Fleck_Factors 
{
    // Useful typedefs.
    typedef typename MT::template CCSF<double> ccsf_double;

    //! Fleck factors are stored in cell-centered-scalar-field.
    ccsf_double fleck;

    //! Constructor.
    Fleck_Factors(rtt_dsxx::SP<MT> mesh) : fleck(mesh) {}
};

} // end namespace rtt_imc

#endif                          // RTT_imc_Fleck_Factors_HH

//---------------------------------------------------------------------------//
//                              end of imc/Fleck_Factors.hh
//---------------------------------------------------------------------------//
