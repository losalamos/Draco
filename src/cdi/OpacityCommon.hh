//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/OpacityCommon.hh
 * \author Kelly Thompson
 * \date   Mon Jan 19 13:41:01 2001
 * \brief  Datatypes needed in GrayOpacity and MultigroupOpacity
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_OpacityCommon_hh__
#define __cdi_OpacityCommon_hh__

namespace rtt_cdi
{

//===========================================================================//
// NUMBER OF MODELS AND REACTIONS
//===========================================================================//

namespace constants
{
/*!
 * \brief Number of models contained in rtt_cdi::Model.
 */
const int num_Models = 6;

/*!
  \brief Number of reaction types contained in rtt_cdi::Reaction.
 */
const int num_Reactions = 3;
}

//===========================================================================//
// ENUMERATIONS USED BY OPACITY CLASSES IN CDI
//===========================================================================//
/*!
 * \brief Physics model used to compute the opacity values.  
 *
 * This enumeration \b must be unnumbered, ie it spans the set [0,N).  The
 * number of models is given by rtt_cdi::constants::num_Models.
 */
enum Model
{
    ROSSELAND, /*!< use Rosseland mean opacities. */
    PLANCK,    /*!< use Plank mean opacities. */
    ANALYTIC,  /*!< use Analytic model opacities. */
    ISOTROPIC, /*!< use Isotropic scattering opacities. */
    THOMSON,   /*!< use Thomson scattering opacities. */
    COMPTON    /*!< use Compton scattering opacities. */
};

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity reaction type stored in this opacity object.
 *
 * This enumeration \b must be unnumbered, ie it spans the set [0,N).  The
 * number of readtion types is given by rtt_cdi::constants::num_Reactions.
 */
enum Reaction
{
    TOTAL,      /*!< Total opacity value (scattering plus absorption). */
    ABSORPTION, /*!< Absorption cross sections only. */
    SCATTERING  /*!< Scattering cross sections only. */
};

} // end namespace rtt_cdi

#endif // __cdi_OpacityCommon_hh__

//---------------------------------------------------------------------------//
// end of cdi/OpacityCommon.hh
//---------------------------------------------------------------------------//
