//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/GrayOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:05:22 2001
 * \brief  GrayOpacity class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GrayOpacity.hh"

namespace rtt_cdi
{

/*!
 * \brief Default GrayOpacity() destructor.
 *
 * This is required to correctly release memory when any
 * object derived from GrayOpacity is destroyed.
 */
GrayOpacity::~GrayOpacity() 
{
    // empty
}

} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
// end of GrayOpacity.cc
//---------------------------------------------------------------------------//
