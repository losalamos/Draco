//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Release.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 14 15:36:52 1999
 * \brief  Release function for imc library.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Release_hh
#define rtt_imc_Release_hh

//===========================================================================//
/*!
 * \page imc_overview Overview of the IMC package
 * \version 3_3_0
 * \author Tom Evans, Todd Urbatsch, Mike Buksas
 *
 * The imc package is a set of class that are used to perform cell-centered,
 * Fleck and Cummings Implicit Monte Carlo transport.  The imc package has
 * capabilities for two types of parallel decompositions: (1) full domain
 * replication; (2) full domain decomposition.  The primary components of the
 * imc package are:
 *
 * \arg Transporter the imc transport solver;
 * \arg Particle the imc particle class;
 * \arg Source the imc source class;
 * \arg Mat_State the imc material state class;
 * \arg Opacity the imc opacity class;
 * \arg Tally the imc tally calss;
 *
 * The frequency treatment is determined, in most cases, upon the Frequency
 * Type (FT) template parameter.  This parameter can be set to Gray_Frequency
 * or Multigroup_Frequency.  Classes that depend on FT use specialization to
 * perform the correct gray or multigroup operations.
 *
 * The imc package makes heavy use of the mc (rtt_mc) package.
 */
//===========================================================================//
/*!
 * \namespace rtt_imc

 * \brief Namespace that contains imc package classes and variables.

 * The rtt_imc namespace contains all classes and variables necessary to use
 * the imc package.  See \ref imc_overview for more information.

 */
//===========================================================================//

#include <string>

namespace rtt_imc 
{
    const std::string release();
}

#endif                          // rtt_imc_Release_hh

//---------------------------------------------------------------------------//
//                              end of Release.hh
//---------------------------------------------------------------------------//
