//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.cc
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfOpacity.hh"
#include "GandolfWrapper.hh"

#include "ds++/Assert.hh"

#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>

namespace rtt_cdi_gandolf
{
/*!
 * \brief Constructor for Gandolf Opacity reader (as a part of CDI).
 */

    // An abrieviated Opacity class is created.  The user must still
    // specify a matID before opacities may be extracted.

GandolfOpacity::GandolfOpacity( string _data_filename )
    : dataFilename ( _data_filename )

    {
	// cout << "In GandolfOpacity::GandolfOpacity(string)" << endl;
	
	// Assert that the opacity data file exists and that it can be read.
	// Let Gandolf do this for now.
	//
	// ofstream infile( gandolfFilename );
	// Insist(!infile,"Could not open Gandalf data file for reading.");
	
	int arrayMatIDs[maxMaterials], numMaterials=0, errorCode=0;

// 	for ( int i=0; i<maxMaterials; ++i ) {
// 	    cout << "GandolfOpacity::GandolfOpacity() arrayMatIDs[" 
// 		 << i << "] = " << arrayMatIDs[i] << endl;
// 	}

	gmatids( dataFilename, arrayMatIDs, maxMaterials, 
		 numMaterials, errorCode );

	// copy arrayMatIDs into the vector matIDs
	matIDs.resize(numMaterials);

	// use sdnolen's example to copy arrayMatIDs into matIDs (STL
	// mechanism).
	for ( int i=0; i<numMaterials; ++i ) {
	    matIDs[i] = arrayMatIDs[i];
// 	    cout << "GandolfOpacity::GandolfOpacity() matIDs[" 
// 		 << i << "] = " << matIDs[i] << endl;
	}
//	cout << endl;
    } // end GandolfOpacity::GandolfOpacity(string)


/*!
 * \brief Constructor for Gandolf Opacity reader (as a part of CDI).
 */
GandolfOpacity::GandolfOpacity( string _data_filename, int _matid )
    : dataFilename ( _data_filename ), matID( _matid )
	  
    {
// 	cout << "In GandolfOpacity::GandolfOpacity(string,int)" << endl;
	
	// Assert that the opacity data file exists and that it can be read.
	// Let Gandolf do this for now.
	//
	// ofstream infile( gandolfFilename );
	// Insist(!infile,"Could not open Gandalf data file for reading.");
	
	int arrayMatIDs[maxMaterials], numMaterials=0, errorCode=0;

	gmatids( dataFilename, arrayMatIDs, maxMaterials,
		 numMaterials, errorCode );

	// copy arrayMatIDs into the vector matIDs
	matIDs.resize(numMaterials);

	// copy the integer array back into our matIDs vector.
	std::copy( arrayMatIDs, arrayMatIDs+numMaterials,
		   matIDs.begin() );
// 	for ( int i=0; i<numMaterials; ++i ) {
// 	    //matIDs[i] = arrayMatIDs[i];
// 	    cout << "GandolfOpacity::GandolfOpacity() matIDs[" 
// 		 << i << "] = " << matIDs[i] << endl;
// 	    cout << "                                arrayMatIDs["
// 		 << i << "] = " << arrayMatIDs[i] << endl;
// 	}
// 	cout << endl;
	
	// search for matID in the vector matIDs.
	bool matID_found = false;
	for (int i=0; i<numMaterials; ++i ) {
	    //	    cout << matID << " == " << matIDs[i] << " ?" << endl;
	    if ( matID == matIDs[i] ) {
		matID_found = true;
		break;
	    }
	}
	// cout << "matID_found = " << matID_found << endl << endl;
	Ensure( matID_found );  // die if we can't find the requeted matID.


	// retrieve keys available for this material from data file.
	//-----------------------------------------------------------
	gkeys( dataFilename, matID, keys, maxKeys, numKeys,
	       errorCode );

// 	cout << "Keys available for material " << matID << "." << endl;
// 	for ( int i=0; i<numKeys; ++i ) {
// 	    cout << "   Keys[" << i << "] = ";
// 	    for( int j=0; j<key_length; ++j ) {
// 		cout << keys[i][j];
// 	    } 
// 	    cout << endl;
// 	}
// 	cout << endl;

	// retrieve size of the data set.
	// -------------------------------
	gchgrids( dataFilename, matID, numTemps, numDensities,
		  numGroupBoundaries, numGrayOpacities,
		  numMGOpacities, errorCode );


    } // end GandolfOpacity::GandolfOpacity(string,int)

/*!
 * \brief Return a Rosseland Mean Gray Opacity value for the user
 *        specified temperature and density.
 */
 double GandolfOpacity::getGray( const double targetTemp, 
				 const double targetDensity )
     {
	 // cout << "In GandolfOpacity::getGray(double,double)." << endl;

	 double grayOpacity = 0.0;

	 // request rosseland gray opacity information.
	 char key[key_length] = "rgray";

	 // ensure that key is available in keys[].
	 bool key_available = false;
	 for ( int i=0; i<numKeys; ++i )
	     // memcmp() compares key[] to keys[i][].  we only compare
	     // the first strlen(key) characters of these strings.  If 
	     // the two match then a zero is returned by memcmp().
	     // This logic may not work if key[] = "tfree" since this
	     // keyword matches both "tfree" and "tfree2".
	     if ( ! memcmp( key, keys[i], strlen(key) ) )
		  {
		      key_available = true;
		      break;
		  }
	 Require( key_available );

	 // if needed, obtain the opacity grid information from file.
	 // do this every time?
	 // if ( temperatures.size() < numTemps ) {
		 int errorCode = 0;
		 ggetgray( dataFilename, matID, key, 
			   temperatures, maxTemps, numTemps,
			   densities, maxDensities, numDensities,
			   grayOpacities, maxGrayOpacities, numGrayOpacities,
			   errorCode );
	 // }

	 // DEBUG PRINTS

// 	 cout << endl
// 	      << "Tabulated values from the IPCRESS file:" << endl;
// 	 cout << endl;
// 	 for( int i=0; i<numTemps; ++i)
// 		 cout << "temperatures[" << i << "] = " 
// 		      << std::scientific << std::setprecision(2)
// 		      << temperatures[i] << " keV" << endl;
// 	 cout << endl;
// 	 for( int i=0; i<numDensities; ++i)
// 		 cout << "densities[" << i << "] = " 
// 		      << std::setw(4)
// 		      << densities[i] << " g/cm^3" << endl;
// 	 cout << endl;
// 	 for( int i=0; i<numGrayOpacities; ++i)
// 		 cout << "GrayOpacities[" << i << "] = " 
// 		      << std::scientific << std::setprecision(8)
// 		      << grayOpacities[i] << " cm^2/g" << endl;
// 	 cout << endl;

	 // gintgrlog expects everything to be in log form.
 	 vector<double> logTemperatures(numTemps);
 	 vector<double> logDensities(numDensities);
 	 vector<double> logGrayOpacities(numGrayOpacities);
	 // replace loop with foreach? or transform? 
 	 for( int i=0; i<numTemps; ++i )
 	     logTemperatures[i] = log( temperatures[i] );
	 for( int i=0; i<numDensities; ++i)
	     logDensities[i] = log( densities[i] );
	 for( int i=0; i<numGrayOpacities; ++i )
	     logGrayOpacities[i] = log( grayOpacities[i] );

	 // send the opacity grid information to the interpolation routine.
	 gintgrlog( logTemperatures, numTemps, logDensities, numDensities,
		    logGrayOpacities, numGrayOpacities, 
		    log(targetTemp), log(targetDensity), grayOpacity );
	 
	 return grayOpacity;
     }


/*!
 * \brief Return Rosseland Multi-group Opacity values for the user
 *        specified temperature and density.
 */
 vector<double> GandolfOpacity::getMG( const double targetTemp, 
				       const double targetDensity )
     {
	 //	 cout << "In GandolfOpacity::getMG(double,double)." << endl;
	 
	 // request rosseland multigroup opacity information.
	 char key[key_length] = "ramg";

	 // ensure that key is available in keys[].
//	 Require ( key_available(key) );
	 bool key_available = false;
	 for ( int i=0; i<numKeys; ++i )
	     // memcmp() compares key[] to keys[i][].  we only compare
	     // the first strlen(key) characters of these strings.  If 
	     // the two match then a zero is returned by memcmp().
	     // This logic may not work if key[] = "tfree" since this
	     // keyword matches both "tfree" and "tfree2".
	     if ( ! memcmp( key, keys[i], strlen(key) ) )
		 // add 1 to strlen(key) to check for space at end of keyword?
		  {
		      key_available = true;
		      break;
		  }
	 Require( key_available );

	 // Test data size:
	 // create a functions to do this?
// 	 bool ReReadData = false;
// 	 if ( numGroupBoundaries > groupBoundaries.size() ) 
// 	     {
// 		 cout << "   Resizing vector<double> groupBoundaries"
// 		      << endl;
// 		 groupBoundaries.resize(numGroupBoundaries);
// 		 ReReadData = true;
// 	     }

// 	 if ( numGroupBoundaries-1 > MGOpacities.size() )
// 	     {
// 		 cout << "   Resizing vector<double> MGOpacities"
// 		      << endl;
// 		 MGOpacities.resize(numGroupBoundaries-1);
// 		 ReReadData = true;
// 	     }

	 // if needed, obtain the opacity grid information from file.
	 // do this every time???
// 	 if ( ReReadData ) {
	 //	 cout << "   Reading MG Opacity Data from file." << endl;
	 int errorCode = 0;
	 ggetmg( dataFilename, matID, key, 
		 temperatures, maxTemps, numTemps,
		 densities, maxDensities, numDensities,
		 groupBoundaries, maxGroupBoundaries, numGroupBoundaries,
		 MGOpacities, maxMGOpacities, numMGOpacities,
		 errorCode );
	 //	 cout << "   Back from reading opacity data file." << endl;
// 	 }

	 // DEBUG PRINTS

// 	 cout << endl
// 	      << "Tabulated values from the IPCRESS file:" << endl;
// 	 cout << endl;
// 	 for( int i=0; i<numTemps; ++i )
// 	     cout << "temperatures[" << i << "] = " 
// 		  << std::scientific << std::setprecision(2)
// 		  << temperatures[i] << " keV" << endl;
// 	 cout << endl;
// 	 for( int i=0; i<numDensities; ++i )
// 	     cout << "densities[" << i << "] = " 
// 		  << std::scientific << std::setprecision(2)
// 		  << densities[i] << " g/cm^3" << endl;
// 	 cout << endl;
// 	 for( int i=0; i<numGroupBoundaries; ++i )
// 	     cout << "groupBoundaries[" << i << "] = " 
// 		  << std::scientific << std::setprecision(2)
// 		  << groupBoundaries[i] << " keV" << endl;
// 	 cout << endl;
//  	 for( int i=0; i<numMGOpacities; ++i)
// 	     cout << "MGOpacities[" << i << "] = " 
// 		  << std::scientific << std::setprecision(8)
// 		  << MGOpacities[i] << " cm^2/g" << endl;
// 	 cout << endl;
// 	 int temp_index=1;
// 	 int dens_index=3;
// 	 int mgop_offset = 
// 	     temp_index * ( numDensities * ( numGroupBoundaries - 1 ) )
// 	     + dens_index * ( numGroupBoundaries - 1 );	     

// 	 cout << endl 
// 	      << "{ ti, di, mi } = { " 
// 	      << temp_index << ", " 
// 	      << dens_index << ", " 
// 	      << mgop_offset << " }" << endl << endl;

// 	 cout << "MGOpacities for T = " << temperatures[temp_index] 
// 	      << " keV and rho = "
// 	      << densities[dens_index] << " g/cm^3 are:" << endl << endl
// 	      << "     Group     E_low           E_high          Opacity" << endl
// 	      << "     -----     -----------     -----------     -----------" << endl;

// 	 for( int i=0; i<numGroupBoundaries-1; ++i)
// 	     cout << "     " << std::setw(5) << i 
// 		  << "     " << std::setprecision(5) << groupBoundaries[i]
// 		  << "     " << std::setprecision(5) << groupBoundaries[i+1]
// 		  << "     " << std::setprecision(5) << MGOpacities[mgop_offset+i]
// 		  << endl;
// 	 cout << endl;

	 // gintmglog expects everything to be in log form.
 	 vector<double> logTemperatures( numTemps );
 	 vector<double> logDensities( numDensities );
	 vector<double> logGroupBoundaries( numGroupBoundaries );
 	 vector<double> logMGOpacities( numMGOpacities );

	 // replace loop with foreach? or transform? 
 	 for( int i=0; i<numTemps; ++i )
 	     logTemperatures[i] = log( temperatures[i] );
	 for( int i=0; i<numDensities; ++i)
	     logDensities[i] = log( densities[i] );
	 for( int i=0; i<numGroupBoundaries; ++i)
	     logGroupBoundaries[i] = log( groupBoundaries[i] );
	 for( int i=0; i<numMGOpacities; ++i )
	     logMGOpacities[i] = log( MGOpacities[i] );
	 
	 // This is the vector we are looking for:
	 vector<double> MGOpacity( numGroupBoundaries-1 );

	 // send the opacity grid information to the interpolation routine.
	 gintmglog( logTemperatures, numTemps, logDensities, numDensities,
		    numGroupBoundaries,
		    logMGOpacities, numMGOpacities, 
		    log(targetTemp), log(targetDensity), MGOpacity );

//  	 cout << endl
// 	      << "MGOpacities for T = " << targetTemp << " keV and rho = "
//  	      << targetDensity << " g/cm^3 are:" << endl << endl;
// 	 for ( int i=0; i<numGroupBoundaries-1; ++i )
// 	     cout << "   Group " << std::setw(2) << i << "  :"
// 		  << "   Opacity = " << std::scientific << std::setprecision(10) 
// 		  << MGOpacity[i] << " kev"
// 		  << endl;
// 	 cout << endl;
	 
	 return MGOpacity;
     }

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
//                              end of GandolfOpacity.cc
//---------------------------------------------------------------------------//
