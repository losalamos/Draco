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


// Vars defining max size of data elements;
 const int maxMaterials = 10; // max number of materials in data file.
 const int maxKeys = 25; // max number of keys per material
 const int key_length = 24; // length of each key descriptor.

 const int maxTemps = 10;
 const int maxDensities = 10;
 const int maxGroupBoundaries = 35;
 const int maxGrayOpacities = maxTemps * maxDensities;
 const int maxMGOpacities = 
     maxGrayOpacities * ( maxGroupBoundaries - 1 );

class GandolfOpacity : public rtt_cdi::Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
    // IPCRESS data filename
    const string dataFilename;
    // list of all material ID's in the data file
    vector<int> matIDs;
    // The material ID that this instance of GandolfOpacity
    // represents.
    int matID;

    // keywords for current material
    char keys[maxKeys][key_length];
    int numKeys; // number of keys for current material
    // size of data grid in Gandolf file
    int numTemps;
    int numDensities;
    int numGroupBoundaries;             // = numGroups+1
    int numGrayOpacities;               // = numTemps * numDensities
    int numMGOpacities;                // = numTemps * numDensities * numGroups    
    
    // Opacity Grid Information
    vector<double> temperatures;
    vector<double> densities;
    vector<double> groupBoundaries;
    vector<double> grayOpacities;
    vector<double> MGOpacities;

  public:

    // CREATORS
    
    // Open the data file and retrieve the matIDs
    GandolfOpacity( string _data_filename );
    
    // Create an Opacity object for a specific material ID.
    GandolfOpacity( string _data_filename, int _matid );

    //defaulted GandolfOpacity(const GandolfOpacity &rhs);
    ~GandolfOpacity()
    {
	// empty
	// cout << "Destroying GandolfOpacity Object." << endl;
    }

    // MANIPULATORS
    
    //defaulted GandolfOpacity& operator=(const GandolfOpacity &rhs);

    // ACCESSORS

    string getDataFilename() { return dataFilename; }
    vector<int> getMatIDs() 
    { 
	// cout << "In GandolfOpacity::getMatIDs()" << endl;
	return matIDs; 
    }

  
    // Assumes Rosseland Opacities
    double getGray( const double temp,
		    const double density );
    vector<double> getMG( const double temp, 
			  const double density );

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfOpacity_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_gandolf/GandolfOpacity.hh
//---------------------------------------------------------------------------//
