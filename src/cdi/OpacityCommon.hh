//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/OpacityCommon.hh
 * \author Kelly Thompson
 * \date   Mon Jan 19 13:41:01 2001
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_OpacityCommon_hh__
#define __cdi_OpacityCommon_hh__

namespace rtt_cdi
{
    // --------------------- //
    // Enumerated data types //
    // --------------------- //

    /*!
     * \brief Physics model used to compute the opacity values.
     */
    enum Model
    {
	Rosseland,
	Plank
    };

    /*!
     * \brief Opacity reaction type stored in this opacity object.
     */
    enum Reaction
    {
	Total,      /*!< Total opacity value (scattering plus absorption). */
	Absorption, /*!< Absorption cross sections only. */
	Scattering  /*!< Scattering cross sections only. */
    };

} // end namespace rtt_cdi

#endif // __cdi_OpacityCommon_hh__

//---------------------------------------------------------------------------//
// end of cdi/OpacityCommon.hh
//---------------------------------------------------------------------------//
