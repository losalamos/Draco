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

#include "SesameTables.hh"

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
	
	const SesameTables SesTabs;

	// The number of material regions.  For a single Eospac object 
	// there is always exactly one region.      
	const int numRegions;

	// regionIndex is uniquely "1" because there will only be one
	// material region per eospac object.
	const int regionIndex;

	// List of materierial IDs
	mutable std::vector< int > matIDs;

	// These next three data members are mutalbe because they
	// specify what data is cached by the Eospac object.  The
	// cached data set may be changed when the user calls a
	// get... function.
	
	// List of numeric identifiers that specify what tables are
	// loaded (P(T,rho), internal energy, etc.)
	mutable std::vector< int > returnTypes;

	// The eos table is cached with this pointer.  The length and
	// contents are set by es1tabs().
	mutable int eosTableLength;
	mutable double *eosTable;
	
      public:
	
	// CREATORS
	
	// The default return type (14) is set so that the data
	// required for obtaining Cv(e-) is already cached in the
	// EOSPAC object.  The matID values can be obtained from the
	// web: http://int.lanl.gov/projects/sdm/win/materials/
	Eospac( const SesameTables& SesTabs );
	//Eospac( int matID, int returnType = 14 );
	//Eospac(const Eospac &rhs);
	~Eospac();
	
	// MANIPULATORS
	
	//Eospac& operator=(const Eospac &rhs);
	
	// ACCESSORS
	
	double getSpecificElectronInternalEnergy(
	    double density, double temperature ) const;

	double getElectronHeatCapacity(
	    double density, double temperature ) const;

	double getSpecificIonInternalEnergy(
	    double density, double temperature ) const;

	double getIonHeatCapacity(
	    double density, double temperature ) const;

	double getNumFreeElectronsPerIon(
	    double density, double temperature ) const;

	double getElectronBasedThermalConductivity(
	    double density, double temperature ) const;

      private:
	
	// IMPLEMENTATION

	double getF( double density, double temperature, 
		     int returnType ) const; 

	double getdFdT( double density, double temperature, 
			int returnType ) const;

 	void expandEosTable () const;

	bool typeFound( int returnType ) const;
    };
    
} // end namespace rtt_cdi_eospac

#endif // __cdi_eospac_Eospac_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/Eospac.hh
//---------------------------------------------------------------------------//
