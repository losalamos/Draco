//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Random_Walk_Sub_Tally.cc
 * \author Thomas M. Evans
 * \date   Tue Jun 24 13:43:28 2003
 * \brief  Random_Walk_Sub_Tally implementation file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Random_Walk_Sub_Tally.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor.
 */
Random_Walk_Sub_Tally::Random_Walk_Sub_Tally()
    : n_random_walks(0),
      sphere_radius(0.0),
      n_spheres_generated(0),
      rw_step_length(0.0)
{
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                   end of Random_Walk_Sub_Tally.cc
//---------------------------------------------------------------------------//
