//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.cc
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  GandolfOpacity class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfOpacity.hh"

#include "ds++/Assert.hh"

#include <string>
#include <cstring>
#include <cmath>

namespace rtt_cdi_gandolf
{
/*!
 * \brief Constructor for Gandolf Opacity reader (as a part of CDI).
 *
 * This is an abrieviated Opacity class is created.  The user must
 * still specify a matID before opacities may be extracted.  Currently,
 * we do not provide a mechanism to specify the material identifier
 * after the cosntruction of a GandolfOpacity object.  This needs to
 * be implemented soon.  For now use the constructor for GandolfOpacity 
 * that requires both a filename and a material identifier.
 *
 */
GandolfOpacity::GandolfOpacity( string _data_filename )
    : dataFilename( _data_filename )
    {
	int arrayMatIDs[maxMaterials];
	int numMaterials = 0;
	int errorCode    = 0;
	
	// Call the wrapper routine.
	gmatids( dataFilename, arrayMatIDs, maxMaterials, 
		 numMaterials, errorCode );

	// copy arrayMatIDs into the vector matIDs
	matIDs.resize(numMaterials);
	std::copy( arrayMatIDs, arrayMatIDs+numMaterials, 
		   matIDs.begin() );

    } // end GandolfOpacity::GandolfOpacity(string)


/*!
 * \brief Constructor for Gandolf Opacity reader (as a part of CDI).
 */
GandolfOpacity::GandolfOpacity( string _data_filename, int _matid )
    : dataFilename ( _data_filename ), matID( _matid )
	  
    {
	int arrayMatIDs[maxMaterials];
	int numMaterials = 0;
	int errorCode    = 0;

	// Call the wrapper routine.
	gmatids( dataFilename, arrayMatIDs, maxMaterials,
		 numMaterials, errorCode );

	// copy arrayMatIDs into the vector matIDs
	matIDs.resize(numMaterials);
	std::copy( arrayMatIDs, arrayMatIDs+numMaterials,
		   matIDs.begin() );
	
	// search for the requested matID in the vector of available matIDs.
	bool matID_found = false;
	for (int i=0; i<numMaterials; ++i ) {
	    if ( matID == matIDs[i] ) {
		matID_found = true;
		break;
	    }
	}
	Require( matID_found );  // die if we can't find the requeted matID.

	// Retrieve keys available for this material from data file.
	gkeys( dataFilename, matID, keys, maxKeys, numKeys,
	       errorCode );

	// Retrieve size of the data set.
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
	 // request rosseland gray opacity information.
	 char key[key_length] = "rgray";

	 // Require that key is available in keys[].
	 Require( key_available( key, keys, numKeys ) );

	 // Obtain the opacity grid information from file.
	 // This only needs to be done the first time getGray is called.
	 int errorCode = 0;
	 if (    temperatures.size()  != numTemps 
	      || densities.size()     != numDensities
	      || grayOpacities.size() != numGrayOpacities )
	     ggetgray( dataFilename, matID, key, 
		       temperatures, maxTemps, numTemps,
		       densities, maxDensities, numDensities,
		       grayOpacities, maxGrayOpacities, numGrayOpacities,
		       errorCode );

	 // The interpolation routine (gintgrlog) expects everything
	 // to be in log form so we create some temporary vectors to
	 // hold the logorithmic temperature, density and opacity data.
 	 vector<double> logTemperatures(numTemps);
 	 vector<double> logDensities(numDensities);
 	 vector<double> logGrayOpacities(numGrayOpacities);
	 // ( should we replace these loops with a foreach? or a transform? call).
 	 for( int i=0; i<numTemps; ++i )
 	     logTemperatures[i] = log( temperatures[i] );
	 for( int i=0; i<numDensities; ++i)
	     logDensities[i] = log( densities[i] );
	 for( int i=0; i<numGrayOpacities; ++i )
	     logGrayOpacities[i] = log( grayOpacities[i] );

	 // Send the opacity grid information to the interpolation routine.
	 double grayOpacity;
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
	 // request rosseland multigroup opacity information.
	 char key[key_length] = "ramg";

	 // Require that key is available in keys[].
	 Require( key_available( key, keys, numKeys ) );

	 // if needed, obtain the opacity grid information from file.
	 int errorCode = 0;
	 if (    temperatures.size()    != numTemps 
	      || densities.size()       != numDensities
	      || groupBoundaries.size() != numGroupBoundaries
	      || MGOpacities.size()     != numMGOpacities )
	 ggetmg( dataFilename, matID, key, 
		 temperatures, maxTemps, numTemps,
		 densities, maxDensities, numDensities,
		 groupBoundaries, maxGroupBoundaries, numGroupBoundaries,
		 MGOpacities, maxMGOpacities, numMGOpacities,
		 errorCode );

	 // The interpolation routine (gintmglog) expects everything
	 // to be in log form so we create some temporary vectors to
	 // hold the logorithmic temperature, density, energy group
	 // boundary and multigroup opacity data.
 	 vector<double> logTemperatures( numTemps );
 	 vector<double> logDensities( numDensities );
	 vector<double> logGroupBoundaries( numGroupBoundaries );
 	 vector<double> logMGOpacities( numMGOpacities );
	 // ( should we replace these loops with a foreach? or a transform? call).
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

	 // Send the opacity grid information to the interpolation routine.
	 gintmglog( logTemperatures, numTemps, logDensities, numDensities,
		    numGroupBoundaries,
		    logMGOpacities, numMGOpacities, 
		    log(targetTemp), log(targetDensity), MGOpacity );

	 return MGOpacity;
     }


 /*! 
  * \brief This function returns "true" if "key" is found in the list
  *        of "keys".
  */
 bool key_available( char key[], char keys[][key_length], int numKeys )
     {
	 // Assume that there is no match and then look for one.
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

	 return key_available;
     } // end of GandolfOpacity::key_available( char, char )

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
//                              end of GandolfOpacity.cc
//---------------------------------------------------------------------------//
