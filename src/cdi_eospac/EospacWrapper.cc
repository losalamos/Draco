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

#include "EospacWrapper.hh"

// for DEBUG only
#include <iostream>

namespace rtt_cdi_eospac
{
    namespace wrapper 
	{

	    int es1tabs( int numRegions,  int numReturnTypes, 
			 int *returnTypes, int *matIDs, 
			 int &eosTableLength, double **eosTable )
		{
		    // Make some assumptions
		    int *llog1 = new int [ numRegions*numReturnTypes ];
		    for ( int i=0; i<numRegions*numReturnTypes; ++i )
			llog1[i] = 0;  // don't use log data.
		    
		    int *iopt   =  new int [ numRegions*numReturnTypes ];
		    for ( int i=0; i<numRegions*numReturnTypes; ++i )
			iopt[i] = 0;      // input integer of smoothing and splitting

		    // options for cat-1 data (see http://laurel.lanl.gov/XCI
		    // /PROJECTS/DATA/eos/UsersDocument/HTML/Overloaded_IOPT-Java1.1.html
		    //  for values.)
		    
		    int lprnt = 1;   // print summary table
		    int iprnt = 0;   // print to screen (not used)
		    int idtab = 0;   // (not used)
		    
		    // Unit convertions:
		    double *unitConversion = new double [ numReturnTypes * 3 ];
		    for ( int i=0; i<numReturnTypes*3; ++i )
			unitConversion[i] = 1.0;
		    
		    // Error reporting:
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
	    
	    int es1info( int &tableIndex, int &regionIndex, double **eosTable, int &llogs, 
			 int &matID, double &atomicNumber, double &atomicMass,
			 double &density0 )
		{
		    // throw this stuff away.
		    int iname,ifile,errorCode;
		    double xcnvt, ycnvt, fcnvt;
		    
		    // Call the info routine
		    es1info_( tableIndex, regionIndex,
			      eosTable, iname, llogs,
			      xcnvt, ycnvt, fcnvt, matID, atomicNumber, 
			      atomicMass, density0, ifile, errorCode );
		    
		    //    if ( xcnvt != 1.0 || ycnvt != 1.0 || fcnvt != 1.0 )
		    //{
		    //	    std::cout << "Unit Conversions are not unity." << std::endl;
		    //     std::cout << "   iname = " << iname << std::endl
		    // 	 << "   ifile = " << ifile << std::endl
		    // 	 << "   errorCode = " << errorCode << std::endl;
		    
		    // 	std::cout << "   xcnvt = " << xcnvt << std::endl
		    // 	 << "   ycnvt = " << ycnvt << std::endl
		    // 	 << "   fcnvt = " << fcnvt << std::endl << std::endl;
		    //	}
		    
		    return errorCode;
		}
	    
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
	    
	    int es1vals( int returnType, int &derivatives, int &interpolation, 
			 double *eosTable, int eosTableLength, int regionIndex,
			 double &xVals, double &yVals, double *returnVals, int &returnSize )
		{
		    // init some values
		    int errorCode = 0;
		    
		    int numZones = 1; // = xVals.size() =?= yVals.size()
		    // Also returnVals should be (numZones,:) where : is determined by 
		    // the value of derivatives and interpolation.
		    
		    int nvalsi = returnSize; // divided by xVals.size()
		    
		    // convert xVals and yVals into an array.
		    double *a_xVals = new double [1];
		    double *a_yVals = new double [1];
		    
		    a_xVals[0] = xVals;
		    a_yVals[0] = yVals;
		    
		    // Call the interpolation routine.
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
//                              end of EospacWrapper.cc
//---------------------------------------------------------------------------//
