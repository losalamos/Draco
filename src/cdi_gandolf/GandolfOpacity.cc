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
#include "GandolfFile.hh"
#include "GandolfException.hh"

#include "ds++/Assert.hh"

#include <iostream>

namespace rtt_cdi_gandolf
{
/*!
 * \brief Constructor for GandolfOpacity object.
 */
GandolfOpacity::GandolfOpacity( 
    const rtt_dsxx::SP<GandolfFile> _spGandolfFile, 
    const int _matid )
    : spGandolfFile( _spGandolfFile ), matID( _matid ), 
      numMaterials( 0 ), numKeys( 0 ), numTemps( 0 ), 
      numDensities( 0 ), numGroupBoundaries( 0 ), 
      numGrayOpacities( 0 ), numMGOpacities( 0 ),
      anyTableLoaded( false),
      grayRosselandTableLoaded( false ), mgRosselandTableLoaded( false ),
      grayPlankTableLoaded( false ), mgPlankTableLoaded( false )
    {
	const std::vector<int> matids = spGandolfFile->getMatIDs();
	if ( ! key_available( matID, matids ) )
	    throw gkeysException( -1 );

	// Retrieve keys available for this material from data file.
	int errorCode;
	wrapper::wgkeys( spGandolfFile->getDataFilename(),
			matID, vkeys, wrapper::maxKeys,
			numKeys, errorCode ); 

	if ( errorCode != 0 ) throw gkeysException ( errorCode );
	
	// Retrieve size of the data set.
	wrapper::wgchgrids( spGandolfFile->getDataFilename(), matID,
			   numTemps, numDensities, numGroupBoundaries,
			   numGrayOpacities, numMGOpacities, errorCode );

	if ( errorCode != 0 ) throw gchgridsException( errorCode );

	// Resize the temperature and density member data.  We only
	// resize the opacity grid when it is first requested.
	logTemperatures.resize(numTemps);
	logDensities.resize(numDensities);
	
// Option: Load complete opacity table here so that the accessor
// functions below can be const.

    } // end of GandolfOpacity::GandolfOpacity(SP<GandolfFile>,int)


/*!
 * \brief Returns the IPCRESS data filename.
 */
 const std::string& GandolfOpacity::getDataFilename() const 
     {
	 return spGandolfFile->getDataFilename();
     }


/*!
 * \brief Return a Rosseland Mean Gray Opacity value for the user
 *        specified temperature and density.
 */
 double GandolfOpacity::getGrayRosseland( 
     const double targetTemperature, const double targetDensity)
     {
	 // request rosseland gray opacity information.
	 const std::string skey = "rgray";
	 
	 // if any table has been loaded and it is no the Gray
	 // Rosseland table then otherTableLoaded will be true.
	 // Otherwise it is false.
	 bool otherTableLoaded = anyTableLoaded &&
	     ( ! grayRosselandTableLoaded );

	 double grayOpacity = getGray( skey,
				       grayRosselandTableLoaded,
				       otherTableLoaded, 
				       targetTemperature,
				       targetDensity ); 
	 
	 // set the appropriate flags to prevent reloading the table.
	 grayRosselandTableLoaded = true;
	 mgRosselandTableLoaded = false;
	 grayPlankTableLoaded = false;
	 mgPlankTableLoaded = false;
	 anyTableLoaded = true;

	 return grayOpacity;

     } // end of getGrayRosseland(double,double);


/*!
 * \brief Return Rosseland Multi-group Opacity values for the user
 *        specified temperature and density.
 */
 std::vector<double> GandolfOpacity::getMGRosseland( 
     const double targetTemperature, 
     const double targetDensity,
     const std::string skey )
     {
	 // request rosseland multigroup opacity information.
	 // const std::string skey = "rtmg";

	 // if any table has been loaded and it is not the multigroup
	 // Rosseland table then otherTableLoaded will be true.
	 // Otherwise it is false.
	 bool otherTableLoaded = anyTableLoaded &&
	     ( ! mgRosselandTableLoaded );

	 // This is the vector we are looking for:
	 std::vector<double> MGOpacity = 
	     getMG( skey, mgRosselandTableLoaded, otherTableLoaded,
		    targetTemperature, targetDensity );
	 
	 // set the appropriate flags to prevent reloading the table.
	 grayRosselandTableLoaded = false;
	 mgRosselandTableLoaded = true;
	 grayPlankTableLoaded = false;
	 mgPlankTableLoaded = false;	
	 anyTableLoaded = true;

	 return MGOpacity;
     }

double GandolfOpacity::getGrayPlank( const double targetTemperature,
				     const double targetDensity )
    {
	// request rosseland gray opacity information.
	const std::string skey = "pgray";
	
	// if any table has been loaded and it is no the Gray
	// Plank table then otherTableLoaded will be true.
	// Otherwise it is false.
	bool otherTableLoaded = anyTableLoaded &&
	    ( ! grayPlankTableLoaded );
	
	double grayOpacity = getGray( skey, 
				      grayPlankTableLoaded,
				      otherTableLoaded,  
				      targetTemperature,
				      targetDensity );  
	
	// set the appropriate flags to prevent reloading the table.
	grayRosselandTableLoaded = false;
	mgRosselandTableLoaded = false;
	
	grayPlankTableLoaded = true;
	mgPlankTableLoaded = false;
	
	anyTableLoaded = true;
	
	return grayOpacity;
    }
 
std::vector<double> GandolfOpacity::getMGPlank( 
    const double targetTemperature, 
    const double targetDensity )
    {
	 // request rosseland multigroup opacity information.
	 const std::string skey = "pmg";

	 // if any table has been loaded and it is not the multigroup
	 // Plank table then otherTableLoaded will be true.
	 // Otherwise it is false.
	 bool otherTableLoaded = anyTableLoaded &&
	     ( ! mgRosselandTableLoaded );

	 std::vector<double> MGOpacity =
	     getMG( skey, mgPlankTableLoaded, otherTableLoaded,
		    targetTemperature, targetDensity );

	 // set the appropriate flags to prevent reloading the table.
	 grayRosselandTableLoaded = false;
	 mgRosselandTableLoaded = false;
	 grayPlankTableLoaded = false;
	 mgPlankTableLoaded = true;	
	 anyTableLoaded = true;

	 return MGOpacity;
    }

std::vector<double> GandolfOpacity::getTemperatureGrid() const
     {
	 std::vector<double> temperatureGrid(numTemps);
	 for ( int i=0; i<numTemps; ++i )
	     temperatureGrid[i] = exp(logTemperatures[i]);
	 return temperatureGrid;
     } 

std::vector<double> GandolfOpacity::getDensityGrid() const
     {
	 std::vector<double> densityGrid(numDensities);
	 for ( int i=0; i<numDensities; ++i )
	     densityGrid[i] = exp(logDensities[i]);
	 return densityGrid;
     } 

/*!
 * \brief Generic routine used to retrieve gray opacity data.  
 */
 double GandolfOpacity::getGray( const std::string &skey,
				 bool grayTableLoaded,
				 const bool otherTableLoaded,
				 const double targetTemperature,
				 const double targetDensity )
    {
	// Require the requested key to be available in the keys list.
	if ( ! key_available( skey, vkeys ) )
	    throw gkeysException( -2 );

	// Resize member containers and fill them with opacity grid
	// data from the data file.  We only need to load this table
	// once. 
	if ( ! grayTableLoaded )
	    {
		// Resize opacity member data as needed
		logGrayOpacities.resize(numGrayOpacities);
		
		// I'm not sure if the temperature/density grid is identical
		// for the MG and the Gray set.  To be safe I will load the
		// Gray temp/density grid into different arrays and compare.
		std::vector<double> logGrayTemperatures(numTemps);
		std::vector<double> logGrayDensities(numDensities);
		
		// Retrieve the gray data 
		int errorCode;
		wrapper::wggetgray( spGandolfFile->getDataFilename(), matID, skey, 
				    logGrayTemperatures, numTemps, numTemps,
				    logGrayDensities, numDensities, numDensities,
				    logGrayOpacities, numGrayOpacities, numGrayOpacities,
				    errorCode );
		
		if ( errorCode != 0 ) 
		    throw ggetgrayException( errorCode );
		
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
		if ( otherTableLoaded )
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
		
	    } // end if ( ! tableLoaded )

	// Send the opacity grid information to the interpolation routine.
	double grayOpacity;
	wrapper::wgintgrlog( logTemperatures,  numTemps, logDensities, numDensities,
			     logGrayOpacities, numGrayOpacities, 
			     log(targetTemperature), log(targetDensity), grayOpacity );
	
	return grayOpacity;

    } // end of getGray()

/*!
 * \breif Generic routine used to retrieve multigroup opacity data.
 */
std::vector<double> GandolfOpacity::getMG( const std::string &skey,
					   bool mgTableLoaded,
					   const bool otherTableLoaded,
					   const double targetTemperature,
					   const double targetDensity )
    {
	// Require that key is available in key.
	if( ! key_available( skey, vkeys ) )
	    throw gkeysException( -2 );

	// Resize member containers and fill them with opacity grid
	// data from the data file.  We only need to load this table
	// once. 
	if ( ! mgTableLoaded ) 
	    {
		// Resize opacity member data as needed
		groupBoundaries.resize(numGroupBoundaries);
		logMGOpacities.resize(numMGOpacities);

		// I'm not sure if the temperature/density grid is identical
		// for the MG and the Gray set.  To be safe I will load the
		// MG temp/density grid into different arrays and compare.
		std::vector<double> logMGtemperatures(numTemps);
		std::vector<double> logMGdensities(numDensities);
		 
		// Retrieve the multi-group data
		int errorCode;
		wrapper::wggetmg( spGandolfFile->getDataFilename(),
				  matID, skey, logMGtemperatures,
				  numTemps, numTemps, logMGdensities,
				  numDensities, numDensities,
				  groupBoundaries, numGroupBoundaries,
				  numGroupBoundaries, logMGOpacities,
				  numMGOpacities, numMGOpacities,
				  errorCode ); 
		
		if ( errorCode != 0 ) 
		    throw ggetmgException( errorCode );

		// The interpolation routine (gintmglog) expects everything
		// to be in log form so we only store the logorithmic
		// temperature, density and opacity data. 
		
		for ( int i=0; i<numTemps; ++i)
		    logMGtemperatures[i] = log( logMGtemperatures[i] );
		for ( int i=0; i<numDensities; ++i)
		    logMGdensities[i] = log( logMGdensities[i] );
		for ( int i=0; i<numMGOpacities; ++i)
		    logMGOpacities[i] = log( logMGOpacities[i] );
		
		// if we have previously loaded the gray
		// opacity table then compare the temperature and
		// density grids.
		if ( otherTableLoaded ) 
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
	     }

	// Send the opacity grid information to the interpolation routine.
	std::vector<double> MGOpacity( numGroupBoundaries-1 );
	wrapper::wgintmglog( logTemperatures, numTemps, logDensities, 
			     numDensities, numGroupBoundaries, 
			     logMGOpacities, numMGOpacities,  
			     log(targetTemperature),
			     log(targetDensity), MGOpacity ); 
	return MGOpacity;

    } // end of vector<double> getMG();

 /*! 
  * \brief This function returns "true" if "key" is found in the list
  *        of "keys".  This is a static member function.
  */
 template < typename T >
 bool GandolfOpacity::key_available( T key, std::vector<T> keys )
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
 bool GandolfOpacity::isSame( const std::vector<double> &v1, 
			      const std::vector<double> &v2 )
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
