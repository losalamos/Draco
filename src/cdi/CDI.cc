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

#include "CDI.hh"

#include "ds++/Assert.hh"
#include "ds++/SP.hh"

namespace rtt_cdi
{

using rtt_dsxx::SP;
    
// CDI constructor:
// Links the opacity object parameter to the CDI class.
    
CDI::CDI( SP<Opacity> _spOpacity ) :
    spOpacity( _spOpacity ) { }
 
 
// Return a list of material ids found in the current file.
vector<int> const CDI::getMatIDs()
    {
	return spOpacity->getMatIDs();
    }
 
 
// Return the interpolated Rosseland Gray Opacity for the specified 
// temperature and density.
double CDI::getGrayOpacity( const double targetTemperature, 
			    const double targetDensity )
    {
	return spOpacity->getGray( targetTemperature, targetDensity );
    }
 
// Return the interpolated Rosseland Gray Opacity for the specified 
// temperature and density.
vector<double> CDI::getMGOpacity( const double targetTemperature, 
				  const double targetDensity )
    {
	return spOpacity->getMG( targetTemperature, targetDensity );
    }
 
} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
