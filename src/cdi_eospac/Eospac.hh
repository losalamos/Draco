//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/Eospac.hh
 * \author Kelly Thompson
 * \date   Mon Apr  2 14:14:29 2001
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_Eospac_hh__
#define __cdi_eospac_Eospac_hh__

#include <vector>

namespace rtt_cdi_eospac
{
    
    //===========================================================================//
    /*!
     * \class Eospac
     *
     */
    // revision history:
    // -----------------
    // 0) original
    // 
    //===========================================================================//
    
    class Eospac 
    {
	
	// NESTED CLASSES AND TYPEDEFS
	
	// DATA
	
	// The number of material regions.  For a single Eospac object 
	// there is always exactly one region.      
	const int numRegions;
	// regionIndex is uniquely "1" because there will only be one
	// material region per eospac object.
	const int regionIndex;

	const int numReturnTypes;

	// List of numeric identifiers that specify what tables are
	// loaded (P(T,rho), internal energy, etc.)
	int *returnTypes;

	// List of materierial IDs
	int *matIDs;

	// The eos table is cached with this pointer.
	int eosTableLength;
	double *eosTable;
	
      public:
	
	// CREATORS
	
	Eospac( int matID );
	Eospac(const Eospac &rhs);
	~Eospac();
	
	// MANIPULATORS
	
	Eospac& operator=(const Eospac &rhs);
	
	// ACCESSORS
	
	double getSpecificElectronInternalEnergy(
	    double density, double temperature ) const;

      private:
	
	// IMPLEMENTATION

	bool typeFound( int returnType ) const;
    };
    
} // end namespace rtt_cdi_eospac

#endif // __cdi_eospac_Eospac_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/Eospac.hh
//---------------------------------------------------------------------------//
