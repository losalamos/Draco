//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/MultigroupOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:06:13 2001
 * \brief  MultigroupOpacity class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MultigroupOpacity.hh"

namespace rtt_cdi
{

    /*!
     * \brief Default Opacity() destructor.
     *
     * This is required to correctly release memory when any
     * object derived from MultigroupOpacity is destroyed.
     */
     MultigroupOpacity::~MultigroupOpacity()
	 {
	     // empty
	 }

} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
// end of MultigroupOpacity.cc
//---------------------------------------------------------------------------//
