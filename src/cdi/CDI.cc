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
    
CDI::CDI( const rtt_dsxx::SP< GrayOpacity > _spGrayOpacity )
    : spGrayOpacity( _spGrayOpacity )
    {
	// empty
    }

CDI::CDI( const rtt_dsxx::SP< MultigroupOpacity > _spMultigroupOpacity )
    : spMultigroupOpacity( _spMultigroupOpacity ) 
    {
	// empty
    }

CDI::CDI(  const rtt_dsxx::SP< GrayOpacity > _spGrayOpacity,
	   const rtt_dsxx::SP< MultigroupOpacity > _spMultigroupOpacity )
    : spGrayOpacity( _spGrayOpacity ), 
      spMultigroupOpacity( _spMultigroupOpacity ) 
    {
	// empty
    }

    // --------- //  
    // Accessors //
    // --------- //

// Provide CDI with access to the full interfaces defined by
// GrayOpacity.hh and MultigroupOpacity.hh

rtt_dsxx::SP< GrayOpacity > CDI::gray() { 
    return spGrayOpacity; }

rtt_dsxx::SP< MultigroupOpacity > CDI::mg() { 
    return spMultigroupOpacity; }
 
} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
