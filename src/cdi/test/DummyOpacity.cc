//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOpacity.cc
 * \author Kelly Thompson
 * \date   Wed Jul 13 16:11:55 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DummyOpacity.hh"

#include <string>

namespace rtt_dummy_opacity
{
/*!
 * \brief blah
 */
DummyOpacity::DummyOpacity( )
    //    : dataFilename ( _data_filename )
    {
	//	matIDs.resize(1);

	cout << "In DummyOpacity::DummyOpacity()" << endl;
	cout << endl;
    } // end DummyOpacity::DummyOpacity()


/*!
 * \brief blah
 */
 double DummyOpacity::getGray( const double temp, const double density )
     {
	 double grayOpacity = 0.0;
	 
	 return grayOpacity;
     }

} // end namespace rtt_dummy_opacity


//---------------------------------------------------------------------------//
//                          end of DummyOpacity.cc
//---------------------------------------------------------------------------//
