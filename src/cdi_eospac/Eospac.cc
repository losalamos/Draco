//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/Eospac.cc
 * \author Kelly Thompson
 * \date   Mon Apr  2 14:14:29 2001
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Eospac.hh"
#include "EospacWrapper.hh"
#include "ds++/Assert.hh"

namespace rtt_cdi_eospac
{

    Eospac::Eospac( int matID )
	: numRegions( 1 ),
	numReturnTypes( 1 ),
	regionIndex( 1 )
	{

	    // one default return type
	    returnTypes = new int [ numReturnTypes ];

	    // The default table value is the temperature-based
	    // electron internal energy per mass E(rho,T) with units
	    // of GJ/Mg.
	    returnTypes[0] = 14; // (enelc)
	    
	    // for now only allow 1 matID but we need an array to
	    // store it in.
	    const int numMatIDs = 1;
	    matIDs = new int [ numMatIDs ];
	    matIDs[0] = matID;

	    // Initialize the eosTable to length 1.
	    eosTableLength = 1;
	    eosTable = new double [ eosTableLength ];

	    // Initialize eosTable and find it's required length
	    int errorCode 
		= wrapper::es1tabs( numRegions, numReturnTypes,
				    returnTypes, matIDs,
				    eosTableLength, &eosTable ); 

	    // Check for errors
	    // Report Errors (if any).
	    if ( errorCode == 0 )
		std::cout << "   Success!" << std::endl;
	    else
		{
		    std::cout << "   Error in es1tabs.  Error Code = " 
			      << errorCode << std::endl;
		    std::cout << wrapper::es1errmsg( errorCode ) 
			      << std::endl;
		    Assert( false );
		}
	    
	    // ---------------------------------------------------- //
	    // May want to use "es1info()" to get info about table? //
	    // May want to use "es1name()" to get table name?       //
	    // ---------------------------------------------------- //

	} // end Eospac::Eospac()

    Eospac::~Eospac()
	{
	    delete [] eosTable;
	    delete [] matIDs;
	    delete [] returnTypes;
	}

    double Eospac::getSpecificElectronInternalEnergy(
	double density, double temperature ) const
	{
	    std::cout << std::endl
		      << "Given a density and temperature return the specific electron internal energy"
		      << std::endl 
		      << "and Heat Capacity (dE/dT)." 
		      << std::endl << std::endl;

	    // Internal Energy has returnTypeIndex == 14
	    const int returnType = 14;

	    // Make sure that Internal energy is available for this
	    // EoS data set.
	    if ( ! typeFound( returnType ) ) 
		Assert( false );

	    // return EOS value plus first derivatives.
	    int derivatives = 2;
	    
	    // use bi-linear interpolation.
	    int interpolation = 1;

// 	    double density = 1.0; // (Mg/m^3)
// 	    double temperature = 5800; // K

	    // int returnSize = 1; // if derivatives==1
	    int returnSize = 3; // if derivatives==2

	    double *returnVals = new double [ returnSize ];

 	    int errorCode 
		= wrapper::es1vals( returnType, derivatives,
				    interpolation, eosTable,
				    eosTableLength, regionIndex,
				    density, temperature, 
				    returnVals, returnSize );
	    
	    Assert ( errorCode == 0 );

	    std::cout << "Density (given)     = " << density << " Mg/m^3" << std::endl
		      << "Temperature (given) = " << temperature << " K" << std::endl
		      << std::endl;
	    
	    double specificElectronInternalEnergy = returnVals[0];

	    std::cout << "Electron Internal Energy (result) = " 
		      <<  specificElectronInternalEnergy // EOS value
		      << " GJ/Mg" << std::endl << std::endl;
	    
	    std::cout << "Heat Capacity (result) = " 
		      << 1000.0 * returnVals[2] / temperature
		      << " J/g/K" << std::endl << std::endl;

	    delete [] returnVals;

	    return specificElectronInternalEnergy;
	}



    bool Eospac::typeFound( int returnType ) const
	{
	    // Loop over all available types.  If the requested
	    // type id matches on in the list then return true.
	    // If we reach the end of the list without a match return
	    // false. 
	    for ( int i=0; i<numReturnTypes; ++i )
		if ( returnType == returnTypes[i] ) return true;
	    return false;
	}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
//                              end of Eospac.cc
//---------------------------------------------------------------------------//
