//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Rnd_Control.cc
 * \author Thomas M. Evans
 * \date   Wed Apr 29 16:08:59 1998
 * \brief  Rnd_Control class implementation file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Rnd_Control_Inline.hh"

namespace rtt_rng 
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!  
 * \brief Rnd_Control class constructor.
 *
 * Constructor for the Rnd_Control class only requires one argument with 3
 * additional options.
 *
 * \param s seed value 
 * \param n total number of independent streams; default = 1.0e9 
 * \param sn stream index; default = 0 
 */
Rnd_Control::Rnd_Control(int s, int n, int sn, int p)
    : d_seed(s), d_number(n), d_streamnum(sn)
{
    Require (d_streamnum <= d_number);
}



} // end namespace rtt_rng




//---------------------------------------------------------------------------//
//                              end of Rnd_Control.cc
//---------------------------------------------------------------------------//
