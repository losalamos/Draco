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
 
    // ---------------------------- //
    // Enumerated EOSPAC data types //
    // ---------------------------- //

	/*!
	 * \brief These are the data types known by EOSPAC
	 */
	enum ES4DataType
	{
	    ES4null,      /*!< null table entry. */
	    
	    // Sesame Catalog 301 entries:
	    
	    ES4prtot,     /*!< temperature-based total pressure. */
	    ES4entot,     /*!< temperature-based total internal energy/mass. */
	    ES4tptot,     /*!< pressure-based total temperature. */
	    ES4tntot,     /*!< energy-based total temperature. */
	    ES4pntot,     /*!< energy-based total pressure. */
	    ES4eptot,     /*!< pressure-based total internal energy/mass. */
	    
	    // Sesame Catalog 303 entries:
	    
	    ES4prion,    /*!< temperature-based ion pressure. */
	    ES4enion,    /*!< temperature-based ion internal energy/mass. */
	    ES4tpion,    /*!< pressure-based ion temperature. */
	    ES4tnion,    /*!< energy-based ion temperature. */
	    ES4pnion,    /*!< energy-based ion pressure. */
	    ES4epion,    /*!< pressure-based ion internal energy/mass. */
	    
	    // Sesame Catalog 304 entries:
	    
	    ES4prelc,   /*!< temperature-based electron pressure. */
	    ES4enelc,   /*!< temperature-based electron internal energy/mass. */
	    ES4tpelc,   /*!< pressure-based electron temperature. */
	    ES4tnelc,   /*!< energy-based electron temperature. */
	    ES4pnelc,   /*!< energy-based electron pressure. */
	    ES4epelc,   /*!< pressure-based electron internal energy/mass. */
	    
	    // Sesame Catalog 306 entries:
	    
	    ES4prcld,  /*!< temperature-based cold curve pressure. */
	    ES4encld,  /*!< temperature-based cold curve internal energy/mass. */
	    
	    // Sesame Catalogs 502-505
	    
	    ES4opacr,  /*!< temperature-based Rosseland mean opacity. */
	    ES4opacc2, /*!< temperature-based electron conductive opacity. */
	    ES4zfree2, /*!< temperature-based number of free electrons per ion. */
	    ES4opacp,  /*!< temperature-based Planck mean opacity. */
	    
	    // Sesame Catalogs 601-605
	    
	    ES4zfree3, /*!< temperature-based number of free electrons per ion. */
	    ES4econde, /*!< temperature-based electron electrical conductivity. */
	    ES4tconde, /*!< temperature-based electron thermal conductivity. */
	    ES4therme, /*!< temperature-based electron thermo-electric coef. */
	    ES4opacc3, /*!< temperature-based electron conductive opacity. */
	    
	    // Sesame Catalog 411
	    
	    ES4tmelt,  /*!< temperature-based melting temperature. */
	    ES4pmelt,  /*!< temperature-based melting pressure. */
	    ES4emelt,  /*!< temperature-based melting internal energy/mass. */
	    
	    // Sesame Catalog 412
	    
	    ES4tfreez, /*!< temperature-based freezing temperature. */
	    ES4pfreez, /*!< temperature-based freezing pressure. */
	    ES4efreez, /*!< temperature-based freezing internal energy/mass. */
	    
	    // Sesame Catalog 431
	    
	    ES4shearm /*!< temperature-based shear modulus. */
	};

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
	 * \brief Map from EOSPAC data type to material identifier.
	 *
	 * Each of the enumerated EOSPAC data types can have a
	 * different SesameTable material identifier.  This vector
	 * contains a list of these material IDs.
	 */
	std::vector< int > matMap;

	/*!
	 * \brief Toggle list to identify which data types have been
	 *        requested. 
	 */
	std::vector< ES4DataType > rtMap;

      public:
	
	// CREATORS
	
	explicit SesameTables();
	
	// ACCESSORS

	// set functions

	SesameTables& prtot( int matID ); 
	SesameTables& entot( int matID ); 
	SesameTables& tptot( int matID ); 
	SesameTables& tntot( int matID ); 
	SesameTables& pntot( int matID ); 
	SesameTables& eptot( int matID ); 
	SesameTables& prion( int matID ); 
	SesameTables& enion( int matID ); 
	SesameTables& tpion( int matID ); 
	SesameTables& tnion( int matID ); 
	SesameTables& pnion( int matID ); 
	SesameTables& epion( int matID ); 
	SesameTables& prelc( int matID ); 
	SesameTables& enelc( int matID ); 
	SesameTables& tpelc( int matID ); 
	SesameTables& tnelc( int matID ); 
	SesameTables& pnelc( int matID ); 
	SesameTables& epelc( int matID ); 
	SesameTables& prcld( int matID ); 
	SesameTables& encld( int matID ); 
	SesameTables& opacr( int matID ); 
	SesameTables& opacc2( int matID );
	SesameTables& zfree2( int matID );
	SesameTables& opacp(  int matID );
	SesameTables& zfree3( int matID );
	SesameTables& econde( int matID );
	SesameTables& tconde( int matID );
	SesameTables& therme( int matID );
	SesameTables& opacc3( int matID );
	SesameTables& tmelt(  int matID );
	SesameTables& pmelt(  int matID );
	SesameTables& emelt(  int matID );
	SesameTables& tfreez( int matID );
	SesameTables& pfreez( int matID );
	SesameTables& efreez( int matID );
	SesameTables& shearm( int matID );

	// More Aliases

	SesameTables& Cve(   int matID ) { return enelc( matID );}
	SesameTables& Cvi(   int matID ) { return enion( matID );}
// 	SesameTables& zfree( int matID ) { return zfree3( matID );}
// 	SesameTables& chie(  int matID ) { return tconde( matID );}

	// Get functions
	
	/*!
	 * \brief Return the material identifier associated with a
	 *        Sesame return type.
	 */
	int matID( ES4DataType returnType ) const;

	// Return the enumerated data type associated with the
	// provided integer index
	ES4DataType returnTypes( int index ) const;

	// Return the number of return types (always 37)
	int getNumReturnTypes() const 
	{
	    return numReturnTypes;
	}

    }; // end class SesameTables
    
} // end namespace rtt_cdi_eospac

#endif  // __cdi_eospac_SesameTables_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_eospac/SesameTables.hh
//---------------------------------------------------------------------------//
