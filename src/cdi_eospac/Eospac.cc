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
#include "EospacException.hh"
//#include "ds++/Assert.hh"
#include <sstream>

// For debug only
#include <iostream>

namespace rtt_cdi_eospac
{

    Eospac::Eospac( const SesameTables& in_SesTabs )
	: numRegions( 1 ), regionIndex( 1 ), eosTableLength( 1 ),
	SesTabs( in_SesTabs )
	{
	    // The initial size for the EOS table is 1.  The call to
	    // expandEosTable() below will reallocate the table to the 
	    // appropriate size.
	    // eosTable = new double [ eosTableLength ];

	    // PreCache the default data type
	    expandEosTable();
		
	    // May want to use "es1info()" to get info about table? //
	    // May want to use "es1name()" to get table name?       //
	    
	} // end Eospac::Eospac()

    Eospac::~Eospac()
	{
	    // allocated by expandEosTable()
	    delete [] eosTable;
	}

    double Eospac::getSpecificElectronInternalEnergy(
	double density, double temperature ) const
	{
	    // Internal Energy has returnType == 14
	    const int returnType = 14; 
	    return getF( density, temperature, returnType );
	}

    double Eospac::getElectronHeatCapacity(
	double density, double temperature ) const
	{
	    // specific Heat capacity is dE/dT at constant pressure.
	    // To obtain the specific electron heat capacity we load
	    // the specific electron internal energy (E) and it's
	    // first derivative w.r.t temperature.

	    // Internal Energy has returnTypeIndex == 14
	    const int returnType = 14;
	    return getdFdT( density, temperature, returnType );
	}

    double Eospac::getSpecificIonInternalEnergy(
	double density, double temperature ) const
	{
	    const int returnType = 8; // (enion)
	    return getF( density, temperature, returnType );
	}

    double Eospac::getIonHeatCapacity(
	double density, double temperature ) const
	{
	    // specific Heat capacity is dE/dT at constant pressure.
	    // To obtain the specific electron heat capacity we load
	    // the specific electron internal energy (E) and it's
	    // first derivative w.r.t temperature.

	    const int returnType = 8; // (enion)
	    return getdFdT( density, temperature, returnType );
	}

    double Eospac::getNumFreeElectronsPerIon(
	double density, double temperature ) const
	{
	    const int returnType = 25; // (zfree3)
	    return getF( density, temperature, returnType );
	}

    double Eospac::getElectronBasedThermalConductivity(
	double density, double temperature ) const
	{
	    const int returnType = 27; // (tconde)
	    return getF( density, temperature, returnType );
	}

    double Eospac::getF( double density, double temperature, 
			 int returnType ) const
	{
	    // Throw an exception if the required return type has not
	    // been loaded by Eospac.
	    if ( ! typeFound( returnType ) ) 
		{
		    std::ostringstream outputString;
		    outputString << "\n\tA request was made for data by getF() "
				 << "for which EOSPAC does not have an\n"
				 << "\tassociated material identifier.\n"
				 << "\tRequested returnType = \"" 
				 << returnType << "\"\n";
		    throw EospacException( outputString.str() );
		}

	    
	    // we don't need derivative values.
	    int derivatives = 1;
	    
	    // use bi-linear interpolation
	    int interpolation = 1;
	    
	    int returnSize = 1; // must be "1" for deriviatives == 1

	    double *returnVals = new double [ returnSize ];

 	    int errorCode 
		= wrapper::es1vals( returnType, derivatives,
				    interpolation, eosTable,
				    eosTableLength, regionIndex,
				    density, temperature, 
				    returnVals, returnSize );
	    
	    if ( errorCode != 0 )
		{
		    std::ostringstream outputString;
		    outputString << "\n\tAn unsuccessful request for EOSPAC data "
				 << "was made by es1vals() from within getF().\n"
				 << "\tThe requested returnType was \"" 
				 << returnType << "\"\n"
				 << "\tThe error code returned was \""
				 << errorCode << "\".\n"
				 << "\tThe associated error message is:\n\t\""
				 << wrapper::es1errmsg( errorCode ) << "\"\n";
		    throw EospacException( outputString.str() );
		}
	    
	    double F = returnVals[0];

	    delete [] returnVals;

	    return F;
	}

    double Eospac::getdFdT( double density, double temperature, 
			    int returnType ) const
	{
	    // dF/dT - assume T is the second independent variable.

	    // EOSPAC actually returns derivative values w.r.t log(T)
	    // so that:

	    // dF/dT = dF/d(log(T)) / T
	    //         ^^^^^^^^^^^^
	    //         This is what EOSPAC returns (we must still
	    //         divide by T)

	    // Throw an exception if the required return type has not
	    // been loaded by Eospac.
	    if ( ! typeFound( returnType ) ) 
		{
		    std::ostringstream outputString;
		    outputString << "\n\tA request was made for data by getdFdT() "
				 << "for which EOSPAC does not have an\n"
				 << "\tassociated material identifier.\n"
				 << "\tRequested returnType = \"" 
				 << returnType << "\"\n";
		    throw EospacException( outputString.str() );
		}
	    
	    // return EOS value plus first derivatives.
	    int derivatives = 2;
	    
	    // use bi-linear interpolation
	    int interpolation = 1;
	    
	    int returnSize = 3; // must be "3" for deriviatives == 2
	    // (EOS value plus 2 derivatives values)

	    double *returnVals = new double [ returnSize ];

 	    int errorCode 
		= wrapper::es1vals( returnType, derivatives,
				    interpolation, eosTable,
				    eosTableLength, regionIndex,
				    density, temperature, 
				    returnVals, returnSize );
	    
	    if ( errorCode != 0 )
		{
		    std::ostringstream outputString;
		    outputString << "\n\tAn unsuccessful request for EOSPAC data "
				 << "was made by es1vals() from within getdFdT().\n"
				 << "\tThe error code returned was \""
				 << errorCode << "\".\n"
				 << "\tThe associated error message is:\n\t\""
				 << wrapper::es1errmsg( errorCode ) << "\"\n";
		    throw EospacException( outputString.str() );
		}
	    
	    double dFdT = returnVals[2]/temperature;

	    delete [] returnVals;

	    return dFdT;
	}

    void Eospac::expandEosTable() const
	{
	    // loop over all possible table definitions.  If a matid
	    // has been assigned to a table then add this information
	    // to the vectors returnTypes[] and matIDs[] which are
	    // used by EOSPAC.	  
	    
	    // MatIDs[] and returnTypes[] are a tuple.

	    for ( int i=0; i<SesTabs.getNumTables(); ++i )
		if ( SesTabs.matID( i ) != 0 )
		    {
			returnTypes.insert( returnTypes.begin(), i ); 
			matIDs.insert( matIDs.begin(),
				       SesTabs.matID( i ) );
		    }

	    // the EOSPAC routines need arrays not vectors so we make
	    // temporary copies of these to vectors.
	    int *returnTypes_a = new int [ returnTypes.size() ];
	    int *matIDs_a = new int [ matIDs.size() ];
	    std::copy( returnTypes.begin(), returnTypes.end(),
		       returnTypes_a );
	    std::copy( matIDs.begin(), matIDs.end(), matIDs_a );

	    eosTable = new double [ eosTableLength ];

	    // Initialize eosTable and find it's required length
	    int errorCode 
		= wrapper::es1tabs( numRegions, returnTypes.size(),
				    returnTypes_a, matIDs_a,
				    eosTableLength, &eosTable ); 

	    // Check for errors
	    if ( errorCode != 0 )
		{
		    std::cout << "   ErrorCode = " << errorCode << std::endl;
 		    std::cout << "   ErrorMessage = " 
			      << wrapper::es1errmsg( errorCode ) << std::endl;

		    std::ostringstream outputString;
		    outputString << "\n\tAn unsuccessful request was made to"
				 << "initialize the EOSPAC table area by "
				 << "expandEosTable().\n"
				 << "\tThe error code returned by es1tabs() was \""
				 << errorCode << "\".\n"
				 << "\tThe associated error message is:\n\t\""
				 << wrapper::es1errmsg( errorCode ) << "\"\n";
		    throw EospacException( outputString.str() );
		}
	    
	    // Clean up temporaries
	    delete [] returnTypes_a;
	    delete [] matIDs_a;
	}

    bool Eospac::typeFound( int returnType ) const
	{
	    // Loop over all available types.  If the requested
	    // type id matches on in the list then return true.
	    // If we reach the end of the list without a match return
	    // false. 
	    for ( int i=0; i<returnTypes.size(); ++i )
		if ( returnType == returnTypes[i] ) return true;
	    return false;
	}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
//                              end of Eospac.cc
//---------------------------------------------------------------------------//
