//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "CDI.hh"

#include "ds++/Assert.hh"

#include <fstream>

namespace rtt_cdi
{
using std::string;
using std::cout;
using std::endl;

// CDI constructor
//----------------------------------------------------------------------
CDI::CDI( OpType _opacity_type, string _opacity_data_filename ) :
    opacityType ( _opacity_type ), 
    opacityDataFilename ( _opacity_data_filename )
    {
	cout << "In CDI::CDI() constructor." << endl;

	// Create the appropriate opacity type.
	// ------------------------------------------------------------
	switch ( opacityType ) 
	    {
	    case Gandolf: // Use Gandolf to obtain opacity data.
		spOpacity = new GandolfOpacity( opacityDataFilename );
		break;
	    default:
		Insist( 0, "CDI: Invalid entry for opacityType." );
		break;
	    }
    }


    // Return the inerpolated Rosseland Gray Opacity for the specified 
    // temperature and density.
    //----------------------------------------------------------------------
    double CDI::getGrayOpacity( const double temp, const double density )
	{
	    double grayOpacity = spOpacity->getGray( temp, density );

	    return grayOpacity;
	}
	


} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
