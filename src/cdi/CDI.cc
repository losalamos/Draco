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
#include "EoS.hh"

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

    CDI::CDI( const rtt_dsxx::SP< const EoS >& in_spEoS )
	: spEoS( in_spEoS ) 
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
    
    // "set" functions

    void CDI::setGrayOpacity( 
	const rtt_dsxx::SP< const GrayOpacity >& in_spGrayOpacity )
	{
	    spGrayOpacity = in_spGrayOpacity;
	}

    void CDI::setMultigroupOpacity( 
	const rtt_dsxx::SP< const MultigroupOpacity >&
	in_spMultigroupOpacity )
	{
	    spMultigroupOpacity = in_spMultigroupOpacity;
	}

    void CDI::setEoS( const rtt_dsxx::SP< const EoS >& in_spEoS )
	{
	    spEoS = in_spEoS;
	}

    // "get" functions

    // Provide CDI with access to the full interfaces defined by
    // GrayOpacity.hh and MultigroupOpacity.hh
    
    const rtt_dsxx::SP< const GrayOpacity >& CDI::gray() const { 
	return spGrayOpacity; }
    
    const rtt_dsxx::SP< const MultigroupOpacity >& CDI::mg() const { 
	return spMultigroupOpacity; }

    // Provide CDI with access to the full interfaces defined by
    // EoS.hh
    
    const rtt_dsxx::SP< const EoS >& CDI::eos() const { 
	return spEoS; }

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
