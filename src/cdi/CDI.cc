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

//#include "GrayOpacity.hh"

namespace rtt_cdi
{

// CDI constructor:
// Links the opacity object parameter to the CDI class.
    
// CDI::CDI( const rtt_dsxx::SP< double > _spOpacity ) :
//     spOpacity( _spOpacity ) { }
 
// CDI::~CDI()
//     {
// 	// empty.
//     }
 
// // Return the name of the Opacity data file.
// std::string CDI::getOpacityDataFilename() const
//     {
// 	return spOpacity->getDataFilename();
//     }
 
// // Return the interpolated Rosseland Gray Opacity for the specified 
// // temperature and density.
//  double CDI::getGrayRosselandOpacity( 
//     const double targetTemperature, 
//     const double targetDensity ) const
//     {
//  	return spOpacity->getGrayRosseland( targetTemperature, 
//  					    targetDensity );
//     }
 
// // Return the interpolated Rosseland Gray Opacity for the specified 
// // temperature and density.
// std::vector<double> CDI::getMGRosselandOpacity( 
//     const double targetTemperature, 
//     const double targetDensity,
//     const std::string skey ) const
//     {
// 	return spOpacity->getMGRosseland( targetTemperature, 
// 					  targetDensity, skey );
//     }
  
// // Return the interpolated Plank Gray Opacity for the specified 
// // temperature and density.
//  double CDI::getGrayPlankOpacity( 
//     const double targetTemperature, 
//     const double targetDensity ) const
//     {
// 	return spOpacity->getGrayPlank( targetTemperature, 
// 					targetDensity );
//     }
 
// // Return the interpolated Plank Gray Opacity for the specified 
// // temperature and density.
// std::vector<double> CDI::getMGPlankOpacity( 
//     const double targetTemperature, 
//     const double targetDensity ) const
//     {
// 	return spOpacity->getMGPlank( targetTemperature, 
// 				      targetDensity );
//     }

} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
