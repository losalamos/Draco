//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/SesameTables.hh
 * \author Kelly Thompson
 * \date   Fri Apr  6 08:57:48 2001
 * \brief  Header file for SesameTables (mapping material IDs
 *         to Sesame table indexes).
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_SesameTables_hh__
#define __cdi_eospac_SesameTables_hh__

#include <vector>

namespace rtt_cdi_eospac
{
    
    //===========================================================================//
    /*!
     * \class SesameTables
     *
     */
    // revision history:
    // -----------------
    // 0) original
    // 
    //===========================================================================//
    
    class SesameTables 
    {
	
	// NESTED CLASSES AND TYPEDEFS
	
	// DATA
	
	const int numTabs; // should be 37!

	int t301;  // Pressure, Temperature, Energy (total)
	int t303;  // Pressure, Temperature, Energy (ion)
	int t304;  // Pressure, Temperature, Energy (electron)
	int t306;  // cold curve P and E.
	int t502;  // Rosseland opacity
	int t503;  // electron conductive opacity.
	int t504;  // free electrons (cat 2)
	int t505;  // Plank opacity
	int t601;  // free electrons (cat 3)
	int t602;  // electrical conductivity (electron)
	int t603;  // thremal conductivity (eletron)
	int t604;  // thrmo-electric coef
	int t605;  // electron conductive opacity
	int t411;  // melting P, T, E
	int t412;  // freezing P, T, E
	int t431;  // shear modulus

	// Keep a mapping between EOSPAC returnTypes and the matIDs
	// associated with each table.
	std::vector< int > RT2MatID;
	
      public:
	
	// CREATORS
	
	SesameTables();
	//SesameTables(const SesameTables &rhs);
	//~SesameTables();
	
	// MANIPULATORS
	
	//SesameTables& operator=(const SesameTables &rhs);
	
	// ACCESSORS

	// set functions by table

	SesameTables& table301( int matID );
	SesameTables& table303( int matID );
	SesameTables& table304( int matID );
	SesameTables& table306( int matID );
	SesameTables& table502( int matID );
	SesameTables& table503( int matID );
	SesameTables& table504( int matID );
	SesameTables& table505( int matID );
	SesameTables& table601( int matID );
	SesameTables& table602( int matID );
	SesameTables& table603( int matID );
	SesameTables& table604( int matID );
	SesameTables& table605( int matID );
	SesameTables& table411( int matID );
	SesameTables& table412( int matID );
	SesameTables& table431( int matID );

	// set functions by data type

	SesameTables& Cve( int matID )
	{
	    return table304( matID );
	}

	SesameTables& Cvi( int matID )
	{
	    return table303( matID );
	}

	SesameTables& zfree( int matID )
	{
	    return table601( matID );
	}

	SesameTables& chie( int matID )
	{ 
	    return table603( matID );
	}

	int matID( int returnType ) const;

	int getNumTables() const 
	{
	    return numTabs;
	}

      private:
	
	// IMPLEMENTATION

	void constructRT2MatID();
    };
    
} // end namespace rtt_cdi_eospac

#endif  // __cdi_eospac_SesameTables_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_eospac/SesameTables.hh
//---------------------------------------------------------------------------//
