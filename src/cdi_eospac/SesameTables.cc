//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/SesameTables.cc
 * \author Kelly Thompson
 * \date   Fri Apr  6 08:57:48 2001
 * \brief  Implementation file for SesameTables (mapping material IDs
 *         to Sesame table indexes).
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "SesameTables.hh"

#include "ds++/Assert.hh"

// Need for DEBUG only
#include <iostream>

namespace rtt_cdi_eospac
{

    // Constructor.

    SesameTables::SesameTables()
	: numReturnTypes( 37 )
	{  
	    // Initialize the material map;
	    matMap.resize( numReturnTypes );
	    for ( int i=0; i<numReturnTypes; ++i )
		matMap[i] = ES4null;
	    
	    // Init a list of return types
	    rtMap.resize( numReturnTypes );
	    for ( int i=0; i<numReturnTypes; ++i )
		rtMap[i] = ES4null;

 	    // EOSPAC has numReturnTypes=37 datatypes.  See
	    // http://laurel.lanl.gov/XCI/PROJECTS/DATA/eos/
	    // UsersDocument/HTML/EOSPAC.html#5.4 for details.
	}

    // Set functions

    SesameTables& SesameTables::prtot( int matID )
	{
	    matMap[ ES4prtot ] = matID;
	    rtMap[ ES4prtot ] = ES4prtot;
	    return *this;
	}

	SesameTables& SesameTables::entot( int matID )
	{
	    matMap[ ES4entot ] = matID;
	    rtMap[ ES4entot ] = ES4entot;
	    return *this;
	}
	SesameTables& SesameTables::tptot( int matID ) 
	{
	    matMap[ ES4tptot ] = matID;
	    rtMap[ ES4tptot ] = ES4tptot;
	    return *this;
	}
	SesameTables& SesameTables::tntot( int matID ) 
	{
	    matMap[ ES4tntot ] = matID;
	    rtMap[ ES4tntot ] = ES4tntot;
	    return *this;
	}
	SesameTables& SesameTables::pntot( int matID ) 
	{
	    matMap[ ES4pntot ] = matID;
	    rtMap[ ES4pntot ] = ES4pntot;
	    return *this;
	}
	SesameTables& SesameTables::eptot( int matID ) 
	{
	    matMap[ ES4eptot ] = matID;
	    rtMap[ ES4eptot ] = ES4eptot;
	    return *this;
	}
	SesameTables& SesameTables::prion( int matID ) 
	{
	    matMap[ ES4prion ] = matID;
	    rtMap[ ES4prion ] = ES4prion;
	    return *this;
	}
	SesameTables& SesameTables::enion( int matID ) 
	{
	    matMap[ ES4enion ] = matID;
	    rtMap[ ES4enion ] = ES4enion;
	    return *this;
	}
	SesameTables& SesameTables::tpion( int matID ) 
	{
	    matMap[ ES4tpion ] = matID;
	    rtMap[ ES4tpion ] = ES4tpion;
	    return *this;
	}
	SesameTables& SesameTables::tnion( int matID ) 
	{
	    matMap[ ES4tnion ] = matID;
	    rtMap[ ES4tnion ] = ES4tnion;
	    return *this;
	}
	SesameTables& SesameTables::pnion( int matID ) 
	{
	    matMap[ ES4pnion ] = matID;
	    rtMap[ ES4pnion ] = ES4pnion;
	    return *this;
	}
	SesameTables& SesameTables::epion( int matID ) 
	{
	    matMap[ ES4enion ] = matID;
	    rtMap[ ES4enion ] = ES4enion;
	    return *this;
	}
	SesameTables& SesameTables::prelc( int matID ) 
	{
	    matMap[ ES4prelc ] = matID;
	    rtMap[ ES4prelc ] = ES4prelc;
	    return *this;
	}
	SesameTables& SesameTables::enelc( int matID ) 
	{
	    matMap[ ES4enelc ] = matID;
	    rtMap[ ES4enelc ] = ES4enelc;
	    return *this;
	}
	SesameTables& SesameTables::tpelc( int matID ) 
	{
	    matMap[ ES4tpelc ] = matID;
	    rtMap[ ES4tpelc ] = ES4tpelc;
	    return *this;
	}
	SesameTables& SesameTables::tnelc( int matID ) 
	{
	    matMap[ ES4tnelc ] = matID;
	    rtMap[ ES4tnelc ] = ES4tnelc;
	    return *this;
	}
	SesameTables& SesameTables::pnelc( int matID ) 
	{
	    matMap[ ES4pnelc ] = matID;
	    rtMap[ ES4pnelc ] = ES4pnelc;
	    return *this;
	}
	SesameTables& SesameTables::epelc( int matID ) 
	{
	    matMap[ ES4epelc ] = matID;
	    rtMap[ ES4epelc ] = ES4epelc;
	    return *this;
	}
	SesameTables& SesameTables::prcld( int matID ) 
	{
	    matMap[ ES4prcld ] = matID;
	    rtMap[ ES4prcld ] = ES4prcld;
	    return *this;
	}
	SesameTables& SesameTables::encld( int matID ) 
	{
	    matMap[ ES4encld ] = matID;
	    rtMap[ ES4encld ] = ES4encld;
	    return *this;
	}
	SesameTables& SesameTables::opacr( int matID ) 
	{
	    matMap[ ES4opacr ] = matID;
	    rtMap[ ES4opacr ] = ES4opacr;
	    return *this;
	}
	SesameTables& SesameTables::opacc2( int matID )
	{
	    matMap[ ES4opacc2 ] = matID;
	    rtMap[ ES4opacc2 ] = ES4opacc2;
	    return *this;
	}
	SesameTables& SesameTables::zfree2( int matID )
	{
	    matMap[ ES4zfree2 ] = matID;
	    rtMap[ ES4zfree2 ] = ES4zfree2;
	    return *this;
	}
	SesameTables& SesameTables::opacp(  int matID )
	{
	    matMap[ ES4opacp ] = matID;
	    rtMap[ ES4opacp ] = ES4opacp;
	    return *this;
	}
	SesameTables& SesameTables::zfree3( int matID )
	{
	    matMap[ ES4zfree3 ] = matID;
	    rtMap[ ES4zfree3 ] = ES4zfree3;
	    return *this;
	}
	SesameTables& SesameTables::econde( int matID )
	{
	    matMap[ ES4econde ] = matID;
	    rtMap[ ES4econde ] = ES4econde;
	    return *this;
	}
	SesameTables& SesameTables::tconde( int matID )
	{
	    matMap[ ES4tconde ] = matID;
	    rtMap[ ES4tconde ] = ES4tconde;
	    return *this;
	}
	SesameTables& SesameTables::therme( int matID )
	{
	    matMap[ ES4therme ] = matID;
	    rtMap[ ES4therme ] = ES4therme;
	    return *this;
	}
	SesameTables& SesameTables::opacc3( int matID )
	{
	    matMap[ ES4opacc3 ] = matID;
	    rtMap[ ES4opacc3 ] = ES4opacc3;
	    return *this;
	}
	SesameTables& SesameTables::tmelt(  int matID )
	{
	    matMap[ ES4tmelt ] = matID;
	    rtMap[ ES4tmelt ] = ES4tmelt;
	    return *this;
	}
	SesameTables& SesameTables::pmelt(  int matID )
	{
	    matMap[ ES4pmelt ] = matID;
	    rtMap[ ES4pmelt ] =ES4pmelt ;
	    return *this;
	}
	SesameTables& SesameTables::emelt(  int matID )
	{
	    matMap[ ES4emelt ] = matID;
	    rtMap[ ES4emelt ] = ES4emelt;
	    return *this;
	}
	SesameTables& SesameTables::tfreez( int matID )
	{
	    matMap[ ES4tfreez ] = matID;
	    rtMap[ ES4tfreez ] = ES4tfreez;
	    return *this;
	}
	SesameTables& SesameTables::pfreez( int matID )
	{
	    matMap[ ES4pfreez ] = matID;
	    rtMap[ ES4pfreez ] = ES4pfreez;
	    return *this;
	}
	SesameTables& SesameTables::efreez( int matID )
	{
	    matMap[ ES4efreez ] = matID;
	    rtMap[ ES4efreez ] = ES4efreez;
	    return *this;
	}
	SesameTables& SesameTables::shearm( int matID )
	{
	    matMap[ ES4shearm ] = matID;
	    rtMap[ ES4shearm ] = ES4shearm;
	    return *this;
	}

    // Get Functions

	// Return the enumerated data type associated with the
	// provided integer index
	ES4DataType SesameTables::returnTypes( int index ) const
	    {
		Assert( index >= 0 && index < numReturnTypes );
 		return rtMap[ index ];
	    }

    int SesameTables::matID( ES4DataType returnType ) const
	{
	    return matMap[ returnType ];
	}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
// end of SesameTables.cc
//---------------------------------------------------------------------------//
