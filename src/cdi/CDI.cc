//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/SP.hh"

#include "GrayOpacity.hh"
#include "MultigroupOpacity.hh"

#include "CDI.hh"

namespace rtt_cdi
{
    
    // ---------------- //
    // CDI constructors //
    // ---------------- //
    
    CDI::CDI( const rtt_dsxx::SP< const GrayOpacity >& in_spGrayOpacity )
	: spGrayOpacity( in_spGrayOpacity )
	{
	    // empty
	}
    
    CDI::CDI( const rtt_dsxx::SP< const MultigroupOpacity >& in_spMultigroupOpacity )
	: spMultigroupOpacity( in_spMultigroupOpacity ) 
	{
	    // empty
	}
    
    CDI::CDI( const rtt_dsxx::SP< const GrayOpacity >& in_spGrayOpacity,
	      const rtt_dsxx::SP< const MultigroupOpacity >& in_spMultigroupOpacity )
	: spGrayOpacity( in_spGrayOpacity ), 
	spMultigroupOpacity( in_spMultigroupOpacity ) 
	{
	    // empty
	}
    
    CDI::~CDI() 
	{
	    // empty
	}
    
    // --------- //  
    // Accessors //
    // --------- //
    
    // Provide CDI with access to the full interfaces defined by
    // GrayOpacity.hh and MultigroupOpacity.hh
    
    const rtt_dsxx::SP< const GrayOpacity >& CDI::gray() const { 
	return spGrayOpacity; }
    
    const rtt_dsxx::SP< const MultigroupOpacity >& CDI::mg() const { 
	return spMultigroupOpacity; }
    
} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
