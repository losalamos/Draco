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
     * \brief This is a helper class for Eospac.  It tells Eospac what 
     *        Sesame data is being requested and what lookup tables to 
     *        use.
     *
     * \sa The web page for <a 
     *     href="http://laurel.lanl.gov/XCI/PROJECTS/DATA/eos/eos.html">EOSPAC</a>.
     *
     * \sa The web page for <a 
     *     href="http://int.lanl.gov/projects/sdm/win/materials/">Eos
     *     Material Identifiers</a>.  This web site also does dynamic
     *     plotting of EoS values.
     *
     * Each sesame material definition has 16 data tables (actually
     * material identifiers) that define its state.  At least one
     * table must be defined for this to be a valid object.  This list 
     * of tables is used by the Eospac constructor to determine what
     * Sesame table data to cache.  There are 37 return types defined
     * by EOSPAC.  Some data tables provide information for more than
     * one return type.
     *
     * \example cdi_eospac/test/tEospac.cc
     */

    // revision history:
    // -----------------
    // 0) original
    // 
    //===========================================================================//
    
    class SesameTables 
    {
	// DATA
	
	/*!
	 * \brief There are 37 return types defined by EOSPAC.
	 */
	const int numReturnTypes; // should be 37

	/*!
	 * \brief There are 16 unique data tables defined by EOSPAC.
	 */
	const int numTables;      // should be 16

	/*!
	 * \brief These table numbers hold the associated material
	 *        identifier associated with this material and
	 *        requested data type. 
	 */
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

      public:
	
	// CREATORS
	
	explicit SesameTables();
	
	// ACCESSORS

	// set functions

	/*!
	 * \brief The set functions assign a material identifier to
	 *        the sesame table specified by the set function
	 *        name. 
	 *
	 * These maybe set using the following conventsion
	 * <pre>
	 *    SesameTables SesTab();
	 *    SesTab.table301( 3717 ).Cve( 23717)...
	 *</pre>
	 */
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

	// Get functions
	
	/*!
	 * \brief Return the material identifier associated with a
	 *        Sesame return type.
	 */
	int matID( int returnType ) const;

	// Return the number of tables (always 16).
	int getNumTables() const 
	{
	    return numTables;
	}

	// Return the number of return types (always 37)
	int getNumReturnTypes() const 
	{
	    return numReturnTypes;
	}

    };
    
} // end namespace rtt_cdi_eospac

#endif  // __cdi_eospac_SesameTables_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_eospac/SesameTables.hh
//---------------------------------------------------------------------------//
