//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/EospacWrapper.cc
 * \author Kelly Thompson
 * \date   Fri Mar 30 15:07:48 2001
 * \brief  Implementation File for EosopacWrapper
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

// cdi_eospac package
#include "EospacWrapper.hh"

// DRACO Packages
#include "ds++/Assert.hh"

namespace rtt_cdi_eospac
{
    namespace wrapper 
	{
	    //====================================================================
	    /*! 
	     * \brief C++ EOSPAC wrapper routines.
	     *
	     * The EOSPAC routines are written in FORTRAN.  The following
	     * are C++ prototypes that mimic the F77 Gandolf functions.
	     * Each of these routines flattens the data types and then
	     * calls the Gandolf library's F77 functions.  
	     */
	    //====================================================================

	    /*!
	     * \brief Based on the requested returnTypes array, allocate
	     * 	  space and cache the required EoS Tables in the array
	     *    eosTable. 
	     */
	    int es1tabs( int numRegions,  int numReturnTypes, 
			 int *returnTypes, int *matIDs, 
			 int &eosTableLength, double **eosTable )
		{
		    // Make some assumptions
		    
		    // don't use log data.
		    int *llog1 = new int [ numRegions*numReturnTypes ];
		    for ( int i=0; i<numRegions*numReturnTypes; ++i )
			llog1[i] = 0; 
		    
		    // input integer of smoothing and splitting
		    int *iopt   =  new int [ numRegions*numReturnTypes ];
		    for ( int i=0; i<numRegions*numReturnTypes; ++i )
			iopt[i] = 0;   

		    // options for cat-1 data (see http://laurel.lanl.gov/XCI
		    // /PROJECTS/DATA/eos/UsersDocument/HTML/Overloaded_IOPT-Java1.1.html
		    //  for values.)
		    
		    int lprnt = 1;   // print summary table
		    int iprnt = 0;   // print to screen (not used)
		    int idtab = 0;   // (not used)
		    
		    // Don't use unit convertions
		    double *unitConversion = new double [ numReturnTypes * 3 ];
		    for ( int i=0; i<numReturnTypes*3; ++i )
			unitConversion[i] = 1.0;
		    
		    // Error reporting.
		    int *errorCodes = new int [ numReturnTypes * numRegions ];
		    for ( int i=0; i<numReturnTypes*numRegions; ++i )
			errorCodes[i] = 0;
		    
		    // Call the fortran routine
		    es1tabs_( llog1, iopt, lprnt, iprnt,
			      numReturnTypes, numRegions,
			      returnTypes, unitConversion, matIDs, idtab, 
			      eosTableLength, eosTable, errorCodes );
		    
		    // Check error code and return
		    int errorCode = 0;
		    for ( int i=0; i<numReturnTypes*numRegions; ++i)
			if ( errorCodes[i] != 0 ) 
			    {
				errorCode = errorCodes[i];
				break;
			    }
		    
		    // Release memory
		    delete [] unitConversion;
		    delete [] errorCodes;
		    delete [] llog1;
		    delete [] iopt;
		    
		    // Return the first error code encountered (or 0 if no errors found).
		    return errorCode;
		}
	    
	    /*!
	     * \brief Return a text error message associated that is
	     *        associated with an EOSPAC error code.
	     */	  
	    std::string es1errmsg( int errorCode )
		{
		    int len = 80;

		    // use this string to init the errormessage (to avoid problems
		    // with f90 interface).
		    // offset by 1 so we dont' kill the trailing \0.
		    std::string errorMessage(len,'_');  // init string with 80 spaces.
		    
		    char *cErrorMessage = new char [len];
		    std::copy( errorMessage.begin(), errorMessage.end(), cErrorMessage );
 		    const char *ccem = cErrorMessage;
		 
		    // Retrieve the text description of errorCode
		    // I would like to call es1errmsg_() directly but
		    // the C doesn't talk to the fortran correctly.
		    // The fortran code can't figure out how long ccem 
		    // is.  Instead I have to call an intermediary F90 
		    // routine (that I compiled into the eospac
		    // library).  This routine takes an additinal
		    // argument, len, so that the fortran knows that
		    // ccem is a character*(len) value and then calls
		    // the EOSPAC routine es1errmsg.
		    kt1errmsg_( errorCode, ccem, len );

		    // Copy to a C++ string container.
		    std::copy( cErrorMessage, cErrorMessage+len,
			       errorMessage.begin() );
		    
		    // Trim trailing whitespace from string.
 		    errorMessage.erase( errorMessage.find_last_of(
 			"abcdefghijklmnopqrstuvwxyz" )+1 );

		    delete [] cErrorMessage;

		    return errorMessage;
		}

	    /*
	     * \brief Retrieve information about the cached EoS data.
	     */
	    int es1info( int &tableIndex, double **eosTable, int &llogs, 
			 int &matID, double &atomicNumber, double &atomicMass,
			 double &density0 )
		{
		    // This is always uniquely unity for my
		    // implementation of Eospac.
		    int regionIndex=1;

		    // Init the error code.
		    int errorCode = 0;

		    // throw this stuff away.
		    int iname, ifile; 

		    // we don't use unit conversions
		    double xcnvt, ycnvt, fcnvt; 
		    
		    // Call the fortran info routine
		    es1info_( tableIndex, regionIndex,
			      eosTable, iname, llogs,
			      xcnvt, ycnvt, fcnvt, matID, atomicNumber, 
			      atomicMass, density0, ifile, errorCode );
		    
		    return errorCode;
		}
	    
	    /*!
	     * \brief Retrieve the table name associated with given
	     *        tableIndex. 
	     */
	    std::string es1name( int &tableID )
		{
		    const int len = 80;
		    // use this string to init the tableName (to avoid problems with
		    // the f90 interface).
		    // offset by 1 so we don't kill the trailing \0.
		    std::string tableName(len-1,' ');
		    char cTableName[len];
		    std::copy( tableName.begin(), tableName.end(), cTableName );
		    const char *cctn = cTableName;
		    
		    // Retrieve the table name from EOSPAC.
		    es1name_( tableID, cctn );
		    
		    // Copy from the const char* container to the string container.
		    std::copy( cTableName, cTableName+len-1, tableName.begin() );
		    
		    // Trim trailing whitespace from string.
		    tableName.erase( tableName.find_last_of(
			"abcdefghijklmnopqrstuvwxyz" )+1 );
		    
		    return tableName;
		}
	    
	    /*!
	     * \brief Retrive EoS values for this material (using the
	     *        specified Sesame tables) corresponding to the
	     *        specified density and temperature values.
	     */
	    int es1vals( const int c_returnType, const int c_derivatives, 
			 const int c_interpolation, double *eosTable, 
			 const int c_eosTableLength,
			 const std::vector< double > &xVals, 
			 const std::vector< double > &yVals, 
			 double *returnVals, const int c_returnSize )
		{
		    // init some values
		    int errorCode = 0;

		    // For our implementation of EOSPAC regionIndex is 
		    // always unity.
		    int regionIndex = 1;
		    
		    // xVals and yVals are a tuple so they must have
		    // the same length.
		    Assert ( xVals.size() == yVals.size() );
		    int numZones = xVals.size();

		    // Also returnVals should be (numZones,:) where : is determined by 
		    // the value of derivatives and interpolation.
		    int nvalsi = c_returnSize / numZones;
		    
		    // convert xVals and yVals into arrays.
		    double *a_xVals = new double [ numZones ];
		    double *a_yVals = new double [ numZones ];
		    
		    std::copy( xVals.begin(), xVals.end(), a_xVals );
		    std::copy( yVals.begin(), yVals.end(), a_yVals );

		    // Remove "const-ness" from const int variables
		    // Neither Fortran library nor the extern "C"
		    // block know about "const" so the C++ compiler
		    // cannot guarantee const-ness once this call is made.

		    int returnType = c_returnType;
		    int derivatives = c_derivatives;
		    int interpolation = c_interpolation;
		    int eosTableLength = c_eosTableLength;

		    // Call the fortran interpolation routine.
		    es1vals_( returnType, derivatives, interpolation, eosTable,
			      eosTableLength, numZones, regionIndex, a_xVals, a_yVals,
			      returnVals, nvalsi, errorCode );
		    
		    // clean up
		    delete [] a_xVals;
		    delete [] a_yVals;
		    
		    return errorCode;
		}
	    
	} // end namespace wrapper

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
// end of EospacWrapper.cc
//---------------------------------------------------------------------------//
