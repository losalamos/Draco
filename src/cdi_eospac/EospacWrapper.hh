//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/EospacWrapper.t.hh
 * \author Kelly Thompson
 * \date   Fri Mar 30 15:07:48 2001
 * \brief  Header file for EospacWrapper
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_EospacWrapper_t_hh__
#define __cdi_eospac_EospacWrapper_t_hh__

#include "EospacWrapper.hh"

#include <string>

// For DEBUG only
#include <iostream>

// EOSPAC uses the following definitions for IRIX64 64-bit:
//
// typedef int BOOLEAN;
// typedef int INTEGER;
// typedef double REAL;
// 
// If another architecture is used we may need to use these
// typedefs:
//
// SunOS:
//
// typedef long BOOLEAN;
// typedef long INTEGER;
// typedef double REAL;

namespace rtt_cdi_eospac
{
    namespace wrapper {

	int es1tabs( int numRegions,  int numReturnTypes, 
		     int *returnTypes, int *matIDs, 
		     int &eosTableLength, double **eosTable );

	std::string es1errmsg( int errorCode );

	int es1info( int &tableIndex, int &regionIndex, 
		     double **eosTable, int &llogs, int &matID,
		     double &atomicNumber, double &atomicMass, 
		     double &density0 );

	std::string es1name( int &tableID );

	int es1vals( int returnType, int &derivatives, 
		     int &interpolation, double *eosTable, 
		     int eosTableLength, int regionIndex, 
		     double &xVals, double &yVals, double *returnVals,
		     int &returnSize );

    } // end namespace wrapper

} // end namespace rtt_cdi_eospac

// Function prototypes for EOSPAC F77 subroutines
//----------------------------------------------------------------------------//

extern "C" {

    void es1tabs_( int *llog1, int *iopt, int &lprnt, int &iprnt, 
		   int &numReturnTypes, int &numRegions,    
		   int *returnTypes, double *ucons, int *matIDs, 
		   int &idtab, int &eosTableLength, double **eosTable,
		   int *errorCodes ); 

    void es1errmsg_( int &errorCode, const char *errorMessage );

    void es1info_( int &tableIndex, int &regionIndex, 
		   double **eosTable, int &iname, int &llogs, 
		   double &xcnvt, double &ycnvt, double &fcnvt, 
		   int &matid, double &znbar, double &anbar, 
		   double &dens0, int &ifile, int &errorCode );

    void es1name_( int &tableID, const char *tableName );

    void es1vals_( int &returnType, int &derivatives, 
		   int &interpolation, double *eosTable, 
		   int &eosTableLength, int &numZones, 
		   int &regionIndex, double *xVals, double *yVals,
		   double *returnVals, int &nvalsi, int &errorCode ); 

} // end of extern "C" block

#endif                          // __cdi_eospac_EospacWrapper_t_hh__

//---------------------------------------------------------------------------//
//                        end of cdi_eospac/EospacWrapper.t.hh
//---------------------------------------------------------------------------//
