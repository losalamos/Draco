//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOpacity.hh
 * \author Kelly Thompson
 * \date   Wed Jul 13 16:11:55 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyOpacity_hh__
#define __cdi_DummyOpacity_hh__

#include <iostream>
#include <vector>

#include "../Opacity.hh"

namespace rtt_dummy_opacity
{

using std::string;
using std::cout;
using std::endl;
using std::vector;
 
//===========================================================================//
/*!
 * \class DummyOpacity
 *
 * \brief
 *
 * \sa Short Abstract goes here
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class DummyOpacity : public rtt_cdi::Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    vector<int> matIDs; // empty for class DummyOpacity.
    
  public:

    // CREATORS
    
    DummyOpacity( );
    //defaulted DummyOpacity(const DummyOpacity &rhs);
    ~DummyOpacity()
    {
	// empty
	cout << "Destroying DummyOpacity Object." << endl;
    }

    // MANIPULATORS
    
    //defaulted DummyOpacity& operator=(const DummyOpacity &rhs);

    // ACCESSORS

    string getDataFilename() 
    { 
	return "no data file associated with this class"; 
    }

    vector<int> getMatIDs()
    {
	cout << "In DummyOpacity::getMatIDs()" << endl;
	return matIDs;
    }

    double getGray( const double temp, const double density );

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_dummy_opacity

#endif // __cdi_DummyOpacity_hh__

//---------------------------------------------------------------------------//
//             end of cdi/test/DummyOpacity.hh
//---------------------------------------------------------------------------//
