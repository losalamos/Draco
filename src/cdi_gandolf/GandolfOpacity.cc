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
#include "GandolfWrapper.hh"

#include "ds++/Assert.hh"

namespace rtt_cdi_gandolf
{

/*!
 * \brief Constructor for Gandolf Opacity reader (as a part of CDI).
 */
GandolfOpacity::GandolfOpacity( const string& _data_filename, 
				const int _matid )
    : dataFilename ( _data_filename ), matID( _matid ),
    numMaterials( 0 ), numKeys( 0 ), numTemps( 0 ), numDensities( 0 ), 
    numGroupBoundaries( 0 ), numGrayOpacities( 0 ), numMGOpacities( 0 ),
    grayRosselandTableLoaded( false ), mgRosselandTableLoaded( false )
	  
    {
	// Gandolf will only look at the first "maxDataFilenameLength"
	// characters of the data filename.  We must require that the
	// given name is less than "maxDataFilenameLength" characters.
	Require( dataFilename.length() < wrapper::maxDataFilenameLength );

	// local variables
	vector<int> matIDs;
	int errorCode = 0;

	// This call to Gandolf validates the datafile and if
	// successful returns a list of material identifiers for which 
	// the materials that exist in the data file.
	wrapper::gmatids( dataFilename, matIDs, wrapper::maxMaterials,
			  numMaterials, errorCode ); 

	// Abort if Gandolf returned an error.
	switch ( errorCode ) {
	case 0: // no errors
	    break;
	case 1: // IPCRESS file not found.
	    Insist( false, "The IPCRESS file was not found.");
	    break;
	case 2: // File is not IPCRESS.
	    Insist( false, "The file does not appear to be in IPCRESS format");
	    break;
	case 3: // Problem reading file
	    Insist( false, "Having trouble reading the IPCRESS file.");
	    break;
	case 4: // No material ID's found in file.
	    Insist( false, "No material ID's were found in the IPCRESS data file.");
	    break;
	case 5: // too many matids found ( nmat > kmat )
	    Insist( false, "Too many materials were found in the data file ( nmat > kmat ).");
	    break;
	default: // unknown error.
	    Insist( false, "Unknown error returned from Gandolf::gmatids().");
	    break;
	}

	// search for the requested matID in the vector of available matIDs.
	Require( key_available( matID, matIDs ) );

	// Retrieve keys available for this material from data file.
	wrapper::gkeys( dataFilename, matID, vkeys, wrapper::maxKeys,
			numKeys, errorCode );

	// Abort if Gandolf returns an error.
	switch ( errorCode ) {
	case 0: // no errors
	    break;
	case 1: // IPCRESS file not found.
	    Insist( false, "The IPCRESS file was not found.");
	    break;
	case 2: // File is not IPCRESS.
	    Insist( false, "The file does not appear to be in IPCRESS format");
	    break;
	case 3: // Problem reading file
	    Insist( false, "Having trouble reading the IPCRESS file.");
	    break;
	case 4: // No keys found for this material.
	    Insist( false, "No keys were found for this material");
	    break;
	case 5: // Too many keys found.
	    Insist( false, "Too many keys for array ( nkeys > kkeys )." );
	    break;
	default: // unknown error.
	    Insist( false, "Unknown error returned from Gandolf::gkeys().");
	    break;
	}

	// Retrieve size of the data set.
	wrapper::gchgrids( dataFilename, matID, numTemps, numDensities,
			   numGroupBoundaries, numGrayOpacities,
			   numMGOpacities, errorCode );

	// Abort if Gandolf returns an error.
	switch ( errorCode ) {
	case 0: // no errors
	    break;
	case -1: // return with etas, not densities
	    Insist( false, "IPCRESS file returned ETAs not Densities.");
	    break;		
	case 1: // IPCRESS file not found.
	    Insist( false, "The IPCRESS file was not found.");
	    break;
	case 2: // File is not IPCRESS.
	    Insist( false, "The file does not appear to be in IPCRESS format");
	    break;
	case 3: // Problem reading file
	    Insist( false, "Having trouble reading the IPCRESS file.");
	    break;
	case 4: // Inconsistent gray grids, mg not checked
	    Insist( false, "Gray grid inconsistent with the temp/density grid.");
	    break;
	case 5: // ngray != nt*nrho, mg not checked
	    Insist( false, "Wrong number of gray opacities found (ngray != nt*nrho)." );
	    break;
	case 6: // inconsistent mg grid.
	    Insist( false, "MG grid inconsistent with the temp/density/hnu grid.");
	    break;
	case 7: //  nmg != nt*nrho*(nhnu-1).
	    Insist( false, "Wrong number of MG opacities found (nmg != nt*nrho*(nhnu-1)).");
	    break;
	default: // unknown error.
	    Insist( false, "Unknown error returned from Gandolf::gchgrids().");
	    break;
	}
	
	// Resize the temperature and density member data.  We only
	// resize the opacity grid when it is first requested.
	logTemperatures.resize(numTemps);
	logDensities.resize(numDensities);

    } // end GandolfOpacity::GandolfOpacity(string,int)

/*!
 * \brief Return a Rosseland Mean Gray Opacity value for the user
 *        specified temperature and density.
 */
 double GandolfOpacity::getGrayRosseland( 
     const double targetTemp, const double targetDensity )
     {
	 // request rosseland gray opacity information.
	 string skey = "rgray";

	 // Require that key is available in keys[].
	 Require( key_available( skey, vkeys ) );

	 // Resize member containers and fill them with opacity grid
	 // data from the data file.  We only need to load this table
	 // once. 
	 if ( ! grayRosselandTableLoaded ) 
	     {
		 // Resize opacity member data as needed
		 logGrayOpacities.resize(numGrayOpacities);

		 // I'm not sure if the temperature/density grid is identical
		 // for the MG and the Gray set.  To be safe I will load the
		 // Gray temp/density grid into different arrays and compare.
		 vector<double> logGrayTemperatures(numTemps);
		 vector<double> logGrayDensities(numDensities);
		 
		 // Retrieve the gray data 
		 int errorCode;
		 wrapper::ggetgray( dataFilename, matID, skey, 
				    logGrayTemperatures, numTemps, numTemps,
				    logGrayDensities, numDensities, numDensities,
				    logGrayOpacities, numGrayOpacities, numGrayOpacities,
				    errorCode );
// 		 wrapper::ggetgray( dataFilename, matID, skey, 
// 				    logGrayTemperatures, wrapper::maxTemps, numTemps,
// 				    logGrayTemperatures, wrapper::maxTemps, numTemps,
// 				    logGrayDensities, wrapper::maxDensities, numDensities,
// 				    logGrayOpacities, wrapper::maxGrayOpacities, numGrayOpacities,
// 				    errorCode );

		 // abort if Gandolf returned an error.
		 switch ( errorCode ) {
		 case 0: // no errors
		     break;
		 case -1: // return with etas, not densities
		     Insist( false, "IPCRESS file returned ETAs not Densities.");
		     break;		
		 case 1: // IPCRESS file not found.
		     Insist( false, "The IPCRESS file was not found.");
		     break;
		 case 2: // File is not IPCRESS.
		     Insist( false, "The file does not appear to be in IPCRESS format");
		     break;
		 case 3: // Problem reading file
		     Insist( false, "Having trouble reading the IPCRESS file.");
		     break;
		 case 4: // Data not found
		     Insist( false, "Requested data not found.  Check nt, nrho, ngray.");
		     break;
		 case 5: // Data larger than allocated arrays.
		     Insist( false, "Data found is larger than allocated array size." );
		     break;
		 case 6: // Data size not equal to nt*nrho
		     Insist( false, "Data size not equal to expected size (ndata != nt*nrho)");
		     break;
		 default: // unknown error.
		     Insist( false, "Unknown error returned from Gandolf::ggetgray().");
		     break;
		 }

		 // The interpolation routine (gintgrlog) expects everything
		 // to be in log form so we only store the logorithmic
		 // temperature, density and opacity data. 
		 
		 for ( int i=0; i<numTemps; ++i )
		     logGrayTemperatures[i] = log( logGrayTemperatures[i] );
		 for ( int i=0; i<numDensities; ++i )
		     logGrayDensities[i] = log( logGrayDensities[i] );
		 for ( int i=0; i<numGrayOpacities; ++i )
		     logGrayOpacities[i] = log( logGrayOpacities[i] );
		 
		 // if we have previously loaded the multigroup
		 // Rosseland opacity table then compare the
		 // temperature and density grids.
		 if ( mgRosselandTableLoaded )
		     {
			 Ensure( isSame( logGrayTemperatures, 
					 logTemperatures ) );
			 Ensure( isSame( logGrayDensities,
					 logDensities ) );
		     }
		 else
		     {
			 logTemperatures = logGrayTemperatures;
			 logDensities = logGrayDensities;
		     }

		 // if the table was loaded sucessfully then set the
		 // appropriate flag to prevent reloading the table.
		 grayRosselandTableLoaded = true;
	     }

	 // Send the opacity grid information to the interpolation routine.
	 double grayOpacity;
	 wrapper::gintgrlog( logTemperatures, numTemps, logDensities, numDensities,
		    logGrayOpacities, numGrayOpacities, 
		    log(targetTemp), log(targetDensity), grayOpacity );
	 
	 return grayOpacity;
     }


/*!
 * \brief Return Rosseland Multi-group Opacity values for the user
 *        specified temperature and density.
 */
 vector<double> GandolfOpacity::getMGRosseland( 
     const double targetTemp, const double targetDensity )
     {
	 // request rosseland multigroup opacity information.
	 string skey = "ramg";

	 // Require that key is available in key.
	 Require( key_available( skey, vkeys ) );

	 // Resize member containers and fill them with opacity grid
	 // data from the data file.  We only need to load this table
	 // once. 
	 if ( ! mgRosselandTableLoaded ) 
	     {
		 // Resize opacity member data as needed
		 groupBoundaries.resize(numGroupBoundaries);
		 logMGOpacities.resize(numMGOpacities);

		 // I'm not sure if the temperature/density grid is identical
		 // for the MG and the Gray set.  To be safe I will load the
		 // MG temp/density grid into different arrays and compare.
		 vector<double> logMGtemperatures(numTemps);
		 vector<double> logMGdensities(numDensities);
		 
		 // Retrieve the multi-group data
		 int errorCode;
// 		 wrapper::ggetmg( dataFilename, matID, skey, 
// 			 logMGtemperatures, wrapper::maxTemps, numTemps,
// 			 logMGdensities, wrapper::maxDensities, numDensities,
// 			 groupBoundaries, wrapper::maxGroupBoundaries, numGroupBoundaries,
// 			 logMGOpacities, wrapper::maxMGOpacities, numMGOpacities,
// 			 errorCode );
		 wrapper::ggetmg( dataFilename, matID, skey, 
			 logMGtemperatures, numTemps, numTemps,
			 logMGdensities, numDensities, numDensities,
			 groupBoundaries, numGroupBoundaries, numGroupBoundaries,
			 logMGOpacities, numMGOpacities, numMGOpacities,
			 errorCode );
		 
		 // abort if Gandolf returned an error.
		 switch ( errorCode ) {
		 case 0: // no errors
		     break;
		 case -1: // return with etas, not densities
		     Insist( false, "IPCRESS file returned ETAs not Densities.");
		     break;		
		 case 1: // IPCRESS file not found.
		     Insist( false, "The IPCRESS file was not found.");
		     break;
		 case 2: // File is not IPCRESS.
		     Insist( false, "The file does not appear to be in IPCRESS format");
		     break;
		 case 3: // Problem reading file
		     Insist( false, "Having trouble reading the IPCRESS file.");
		     break;
		 case 4: // Data not found
		     Insist( false, "Requested data not found.  Check nt, nrho, nhnu and ndata.");
		     break;
		 case 5: // Data larger than allocated arrays.
		     Insist( false, "Data found is larger than allocated array size." );
		     break;
		 case 6: // Data size not equal to nt*nrho
		     Insist( false, "Data size not equal to expected size (ndata != nt*nrho*(nhnu-1))");
		     break;
		 default: // unknown error.
		     Insist( false, "Unknown error returned from Gandolf::ggetmg().");
		     break;
		 }

		 // The interpolation routine (gintmglog) expects everything
		 // to be in log form so we only store the logorithmic
		 // temperature, density and opacity data. 

		 for ( int i=0; i<numTemps; ++i)
		     logMGtemperatures[i] = log( logMGtemperatures[i] );
		 for ( int i=0; i<numDensities; ++i)
		     logMGdensities[i] = log( logMGdensities[i] );
		 for ( int i=0; i<numMGOpacities; ++i)
		     logMGOpacities[i] = log( logMGOpacities[i] );
		 
		 // if we have previously loaded the gray Rosseland
		 // opacity table then compare the temperature and
		 // density grids.
		 if ( grayRosselandTableLoaded ) 
		     {
			 Ensure( isSame( logMGtemperatures,
					 logTemperatures ) );
			 Ensure( isSame( logMGdensities, 
					 logDensities    ) ) ;
		     }
		 else
		     {
			 logTemperatures = logMGtemperatures;
			 logDensities    = logMGdensities;
		     }

		 // if the table was loaded sucessfully then set the
		 // appropriate flat to prevent reloading the table.
		 mgRosselandTableLoaded = true;
	     }

	 // This is the vector we are looking for:
	 vector<double> MGOpacity( numGroupBoundaries-1 );

	 // Send the opacity grid information to the interpolation routine.
	 wrapper::gintmglog( logTemperatures, numTemps, logDensities,
			     numDensities, numGroupBoundaries,
		    logMGOpacities, numMGOpacities, 
		    log(targetTemp), log(targetDensity), MGOpacity );

	 return MGOpacity;
     }


 /*! 
  * \brief This function returns "true" if "key" is found in the list
  *        of "keys".  This is a static member function.
  */
 template < typename T >
 bool GandolfOpacity::key_available( T key, vector<T> keys )
     {
	 // Loop over all available keys.  If the requested key
	 // matches one in the list return true.  If we reach the end
	 // of the list without a match return false.
	 for ( int i=0; i<keys.size(); ++i )
	     if ( key == keys[i] ) return true;
	 return false;

     } // end of GandolfOpacity::key_available( string, vector<string> )

 /*!
  * \brief This function returns "true" if the two vectors are the
  *        same.
  */
 bool GandolfOpacity::isSame( const vector<double> &v1, 
			      const vector<double> &v2 )
     {
	 const double TOL = 1.0e-10;

	 Require( v1.size() == v2.size() );

	 for ( int i=0; i < v1.size(); ++i )
	     if ( fabs(v1[i]-v2[i])/v1[i] > TOL )
		 return false;

	 return true;

     } // end of GandolfOpacity::isSame( vector<double>, vector<double> )
  

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
//                              end of GandolfOpacity.cc
//---------------------------------------------------------------------------//
