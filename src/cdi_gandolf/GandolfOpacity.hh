//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.hh
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  CDI component for accessing IPCRESS files through Gandolf
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfOpacity_hh__
#define __cdi_gandolf_GandolfOpacity_hh__

#include <iostream>
#include <vector>

#include "cdi/Opacity.hh"

namespace rtt_cdi_gandolf
{

using std::string;
using std::cout;
using std::endl;
using std::vector;
 
//===========================================================================//
/*!
 * \class GandolfOpacity
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

class GandolfOpacity : public rtt_cdi::Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
    const string dataFilename;
    vector<int> matIDs;
    
  public:

    // CREATORS
    
    GandolfOpacity( string _data_filename );
    //defaulted GandolfOpacity(const GandolfOpacity &rhs);
    ~GandolfOpacity()
    {
	// empty
	cout << "Destroying GandolfOpacity Object." << endl;
    }

    // MANIPULATORS
    
    //defaulted GandolfOpacity& operator=(const GandolfOpacity &rhs);

    // ACCESSORS

    string getDataFilename() { return dataFilename; };
    vector<int> getMatIDs() 
    { 
	cout << "In GandolfOpacity::getMatIDs()" << endl;
	return matIDs; 
    };

     double getGray( const double temp, const double density );

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi_gandolf

#endif                          // __cdi_gandolf_GandolfOpacity_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_gandolf/GandolfOpacity.hh
//---------------------------------------------------------------------------//
