//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/EoS.cc
 * \author Kelly Thompson
 * \date   Mon April 16 10:08:42 2001
 * \brief  EoS class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "EoS.hh"

namespace rtt_cdi
{
    /*!
     * \brief Default EoS() destructor.
     *
     * This is required to correctly release memory when any
     * object derived from EoS is destroyed.
     */
    EoS::~EoS() 
	{
	    // empty
	}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
// end of EoS.cc
//---------------------------------------------------------------------------//
