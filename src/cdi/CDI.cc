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
#include "ds++/SP.hh"

//#include <fstream>

namespace rtt_cdi
{

using std::string;
using std::cout;
using std::endl;
using rtt_dsxx::SP;

// CDI constructor
//----------------------------------------------------------------------
CDI::CDI( SP<Opacity> _spOpacity ) :
    spOpacity( _spOpacity )
    {
	cout << "In CDI::CDI() constructor." << endl;
    }


// Return a list of material ids found in the current file.
 vector<int> CDI::getMatIDs()
     {
	 cout << "In CDI::getMatIDs()" << endl;
	 return spOpacity->getMatIDs();
     }
 

    // Return the interpolated Rosseland Gray Opacity for the specified 
    // temperature and density.
    //----------------------------------------------------------------------
    double CDI::getGrayOpacity( const double temp, 
				const double density )
	{
	    double grayOpacity = spOpacity->getGray( temp, density );
	    return grayOpacity;
	}

    // Return the interpolated Rosseland Gray Opacity for the specified 
    // temperature and density.
    //----------------------------------------------------------------------
    vector<double> CDI::getMGOpacity( const double temp, 
				      const double density )
	{
	    vector<double> MGOpacity = spOpacity->getMG( temp, density );
	    return MGOpacity;
	}

} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
