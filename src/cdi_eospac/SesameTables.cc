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

// Need for DEBUG only
#include <iostream>

namespace rtt_cdi_eospac
{

    // Constructor.

    SesameTables::SesameTables()
	: t301( 0 ), t303( 0 ), t304( 0 ), t306( 0 ), t502( 0 ), 
	t503( 0 ), t504( 0 ), t505( 0 ), t601( 0 ), t602( 0 ), 
	t603( 0 ), t604( 0 ), t605( 0 ), t411( 0 ), t412( 0 ),
	t431( 0 ), numTables( 16 ), numReturnTypes( 37 )
	{  
 	    // EOSPAC has 37 datatypes.  See
	    // http://laurel.lanl.gov/XCI/PROJECTS/DATA/eos/
	    // UsersDocument/HTML/EOSPAC.html#5.4 for details.
	}

    // Set functions

    SesameTables& SesameTables::table301( int matID )
	{
	    t301 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table303( int matID )
	{
	    t303 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table304( int matID )
	{
	    t304 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table306( int matID )
	{
	    t306 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table502( int matID )
	{
	    t502 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table503( int matID )
	{
	    t503 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table504( int matID )
	{
	    t504 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table505( int matID )
	{
	    t505 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table601( int matID )
	{
	    t601 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table602( int matID )
	{
	    t602 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table603( int matID )
	{
	    t603 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table604( int matID )
	{
	    t604 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table605( int matID )
	{
	    t605 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table411( int matID )
	{
	    t411 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table412( int matID )
	{
	    t412 = matID;
	    return *this;
	}

    SesameTables& SesameTables::table431( int matID )
	{
	    t431 = matID;
	    return *this;
	}

    // Get Functions

    int SesameTables::matID( int returnType ) const
	{
	    // EOSPAC has 37 datatypes.  See
	    // http://laurel.lanl.gov/XCI/PROJECTS/DATA/eos/
	    // UsersDocument/HTML/EOSPAC.html#5.4 for details.
	    switch( returnType )
		{
		case 1: case 2: case 3: case 4: case 5: case 6:
		    return t301;
		case 7: case 8: case 9: case 10: case 11: case 12:
		    return t303;
		case 13: case 14: case 15: case 16: case 17: case 18:
		    return t304;
		case 19: case 20:
		    return t306;
		case 21:
		    return t502;
		case 22:
		    return t503;
		case 23:
		    return t504;
		case 24:
		    return t505;
		case 25:
		    return t601;
		case 26:
		    return t602;
		case 27:
		    return t603;
		case 28:
		    return t604;
		case 29:
		    return t605;
		case 30: case 31: case 32:
		    return t411;
		case 33: case 34: case 35:
		    return t412;
		case 36:
		    return t431;
		}
	    // No match -- something is wrong.
	    return 0;
	}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
// end of SesameTables.cc
//---------------------------------------------------------------------------//
