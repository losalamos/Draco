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

#ifndef __imc_Release_hh__
#define __imc_Release_hh__

//===========================================================================//
/*!
 * \page imc_overview Overview of the IMC package
 * \version 2_0_0
 * \author Tom Evans, Todd Urbatsch

 * The imc package is a set of class that are used to perform cell-centered,
 * Fleck and Cummings Implicit Monte Carlo transport.  The imc package has
 * capabilities for two types of parallel decompositions: (1) full
 * replication; (2) full domain decomposition.  The primary components of the
 * imc package are:

 * \arg Transporter the imc transport solver;
 * \arg Particle the imc particle class;
 * \arg Source the imc source class;
 * \arg Mat_State the imc material state class;
 * \arg Opacity the imc opacity class;
 * \arg Tally the imc tally calss;

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

#endif                          // __imc_Release_hh__

//---------------------------------------------------------------------------//
//                              end of Release.hh
//---------------------------------------------------------------------------//
