//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/Eospac.cc
 * \author Kelly Thompson
 * \date   Mon Apr  2 14:14:29 2001
 * \brief  
 * \note   Copyright (C) 2001-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

// cdi_eospac dependencies
#include "Eospac.hh"
#include "EospacWrapper.hh"
#include "EospacException.hh"

// Other Draco dependencies
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"

// C++ standard library dependencies
#include <sstream>

// For debug only
#include <iostream>

namespace rtt_cdi_eospac
{
    // ------------ //
    // Constructors //
    // ------------ //
    
    /*!
     * \brief The constructor for Eospac.
     *
     * \sa The definition of rtt_cdi_eospac::SesameTables.
     *
     * \param SesTabs A rtt_cdi_eospac::SesameTables object that
     * defines what data tables will be available for queries from
     * the Eospac object. 
     *
     */

    // In this implementation of EOSPAC there is never more than one
    // region per Eospac object.  The associated values are set here
    // and never changed.

    // The eosTableLength is initalized to 1 and then changed by
    // es1tabs_() to be the required length.
    Eospac::Eospac( SesameTables const & in_SesTabs )
	: SesTabs( in_SesTabs ),
          numRegions( 1 ),
          regionIndex( 1 ),
          eosTableLength( 1 )
	{
	    // Eospac can only be instantiated if SesameTables is
	    // provided.  If SesameTables is invalid this will be
	    // caught in expandEosTable();

	    // PreCache the default data type
	    expandEosTable();
		
	    // May want to use "es1info()" to get info about table? //
	    // May want to use "es1name()" to get table name?       //
	    
	} // end Eospac::Eospac()

    /*!
     * \brief Default Eospac() destructor.
     *
     * This is required to correctly release memeroyt when an
     * Eospac object is destroyed.  We define the destructor in
     * the implementation file to avoid including the unnecessary
     * header files.
     */
    Eospac::~Eospac()
	{
	    // allocated by expandEosTable()
	    delete [] eosTable;
	}
    
    // --------- //
    // Accessors //
    // --------- //

    double Eospac::getSpecificElectronInternalEnergy(
	double temperature, double density ) const
	{
	    const ES4DataType returnType = ES4enelc; 
	    return getF( dbl_v1(temperature), dbl_v1(density),
			 returnType )[0];
	}

    std::vector< double > Eospac::getSpecificElectronInternalEnergy(
	const std::vector< double >& vtemperature,
	const std::vector< double >& vdensity ) const
	{
	    const ES4DataType returnType = ES4enelc; 
	    return getF( vtemperature, vdensity, returnType );
	}

    double Eospac::getElectronHeatCapacity(
	double temperature, double density ) const
	{
	    // specific Heat capacity is dE/dT at constant pressure.
	    // To obtain the specific electron heat capacity we load
	    // the specific electron internal energy (E) and it's
	    // first derivative w.r.t temperature.

	    const ES4DataType returnType = ES4enelc;
 	    return getdFdT( dbl_v1(temperature), dbl_v1(density),
			    returnType )[0];
	}

    std::vector< double > Eospac::getElectronHeatCapacity(
	const std::vector< double >& vtemperature, 
	const std::vector< double >& vdensity ) const
	{
	    // specific Heat capacity is dE/dT at constant pressure.
	    // To obtain the specific electron heat capacity we load
	    // the specific electron internal energy (E) and it's
	    // first derivative w.r.t temperature.

	    const ES4DataType returnType = ES4enelc;
 	    return getdFdT( vtemperature, vdensity, returnType );
	}

    double Eospac::getSpecificIonInternalEnergy(
	double temperature, double density ) const
	{
	    const ES4DataType returnType = ES4enion;
	    return getF( dbl_v1(temperature), dbl_v1(density),
			 returnType )[0];
	}

    std::vector< double > Eospac::getSpecificIonInternalEnergy(
	const std::vector< double >& vtemperature, 
	const std::vector< double >& vdensity ) const
	{
	    const ES4DataType returnType = ES4enion;
	    return getF( vtemperature, vdensity, returnType );
	}

    double Eospac::getIonHeatCapacity(
	double temperature, double density ) const
	{
	    // specific Heat capacity is dE/dT at constant pressure.
	    // To obtain the specific electron heat capacity we load
	    // the specific electron internal energy (E) and it's
	    // first derivative w.r.t temperature.

	    const ES4DataType returnType = ES4enion; 
	    return getdFdT( dbl_v1(temperature), dbl_v1(density), 
			    returnType )[0];
	}

    std::vector< double > Eospac::getIonHeatCapacity(
	const std::vector< double >& vtemperature,
	const std::vector< double >& vdensity ) const
	{
	    // specific Heat capacity is dE/dT at constant pressure.
	    // To obtain the specific electron heat capacity we load
	    // the specific electron internal energy (E) and it's
	    // first derivative w.r.t temperature.

	    const ES4DataType returnType = ES4enion;
	    return getdFdT( vtemperature, vdensity, returnType );
	}
    
    double Eospac::getNumFreeElectronsPerIon(
	double temperature, double density ) const
	{
	    const ES4DataType returnType = ES4zfree3; // (zfree3)
	    return getF( dbl_v1(temperature), dbl_v1(density), 
			 returnType )[0];
	}
    
    std::vector< double > Eospac::getNumFreeElectronsPerIon(
	const std::vector< double >& vtemperature,
	const std::vector< double >& vdensity ) const
	{
	    const ES4DataType returnType = ES4zfree3; // (zfree3)
	    return getF( vtemperature, vdensity, returnType );
	}
    
    double Eospac::getElectronThermalConductivity(
	double temperature, double density ) const
	{
	    const ES4DataType returnType = ES4tconde; // (tconde)
	    return getF( dbl_v1(temperature), dbl_v1(density), 
			 returnType )[0];
	}
    
    std::vector< double > Eospac::getElectronThermalConductivity(
	const std::vector< double >& vtemperature, 
	const std::vector< double >& vdensity ) const
	{
	    const ES4DataType returnType = ES4tconde; // (tconde)
	    return getF( vtemperature, vdensity, returnType );
	}

// ------- //
// Packing //
// ------- //

/*!
 * Pack the Eospac state into a char string represented by a
 * vector<char>. This can be used for persistence, communication, etc. by
 * accessing the char * under the vector (required by implication by the
 * standard) with the syntax &char_string[0]. Note, it is unsafe to use
 * iterators because they are \b not required to be char *.
 */
std::vector<char> Eospac::pack() const
{
    using std::vector;
    using std::string;

    string msg = "eospac::pack not fully implemented";
    vector<char> packed_descriptor;
    rtt_dsxx::pack_data(msg, packed_descriptor);

    // ************************************************************
    // See GandolfGrayOpacity.cc for how to finish this function.
    // ************************************************************

    // determine the total size: 3 ints (reaction, model, material id) + 2
    // ints for packed_filename size and packed_descriptor size + char in
    // packed_filename and packed_descriptor
    int size = packed_descriptor.size() + 1;

    // make a container to hold packed data
    vector<char> packed(size);

    // make a packer and set it
    rtt_dsxx::Packer packer;
    packer.set_buffer(size, &packed[0]);

    // pack the descriptor
    packer << static_cast<int>(packed_descriptor.size());
    for( size_t i = 0; i < packed_descriptor.size(); i++ )
	packer << packed_descriptor[i];

    Ensure (packer.get_ptr() == &packed[0] + size);

    return packed;
}


    // -------------- //
    // Implementation //
    // -------------- //

    /*!
     * \brief Retrieves the EoS data associated with the returnType 
     *        specified and the given (density, temperature) tuples.
     *
     * Each of the public access functions calls either getF() or
     * getdFdT() after assigning the correct value to
     * "returnType".
     *
     * \param vdensity A vector of density values (g/cm^3).
     * \param vtemperature A vector of temperature values (K).
     * \param returnType The integer index that corresponds to the 
     *        type of data being retrieved from the EoS tables.
     */
    std::vector< double > Eospac::getF( 
	const std::vector< double >& vtemperature,
	const std::vector< double >& vdensity,
	ES4DataType returnType ) const
	{
	    // The density and vector parameters must be a tuple.
	    Require( vdensity.size() == vtemperature.size() );

	    // The returnType must be in the range [1,36]
	    // Require( returnType > 0 && returnType <= 36 );

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
		    throw EospacUnknownDataType( outputString.str() );
		}

	    // Convert temperatures from keV to degrees Kelvin.
	    std::vector< double > vtempsKelvin = vtemperature;
	    for( size_t i=0; i<vtemperature.size(); ++i )
		vtempsKelvin[i] = keV2K( vtemperature[i] );
	    
	    // we don't need derivative values.
	    const int derivatives = 1;
	    
	    // use bi-linear interpolation
	    const int interpolation = 1;
	    
	    // There is one piece of returned information for each
	    // (density, temperature) tuple.
	    const int returnSize = vdensity.size();

	    std::vector< double > returnVals( returnSize );
	    
 	    int errorCode 
		= wrapper::es1vals( returnType, derivatives,
				    interpolation, eosTable,
				    eosTableLength,
				    vdensity, vtempsKelvin, 
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
				 << wrapper::es1errmsg( errorCode )
				 << "\"\n";

		    // This is a fatal exception right now.  It might
		    // be useful to throw a specific exception that is 
		    // derived from EospacException.  The host code
		    // could theoretically catch such an exception,
		    // fix the problem and then continue.
		    throw EospacException( outputString.str() );
		}
	    
	    return returnVals;
	}

    /*!
     * \brief Retrieves the EoS data associated with the returnType 
     *        specified and the given (density, temperature) tuples.
     *
     * Each of the public access functions calls either getF() or
     * getdFdT() after assigning the correct value to
     * "returnType".
     *
     * \param vdensity A vector of density values (g/cm^3).
     * \param vtemperature A vector of temperature values (K).
     * \param returnType The integer index that corresponds to the 
     *        type of data being retrieved from the EoS tables.
     */
    std::vector< double > Eospac::getdFdT( 
	const std::vector< double >& vtemperature, 
	const std::vector< double >& vdensity,
	ES4DataType returnType ) const
	{
	    // The density and vector parameters must be a tuple.
	    Require( vdensity.size() == vtemperature.size() );

	    // The returnType must be in the range [1,36]
	    // Require( returnType > 0 && returnType <= 36 );

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
		    throw EospacUnknownDataType( outputString.str() );
		}

	    // Convert temperatures from keV to degrees Kelvin.
	    std::vector< double > vtempsKelvin = vtemperature;
	    for ( size_t i=0; i<vtemperature.size(); ++i )
		vtempsKelvin[i] = keV2K( vtemperature[i] );
	    
	    // return EOS value plus first derivatives.
	    const int derivatives = 2;
	    
	    // use bi-linear interpolation
	    const int interpolation = 1;
	    
	    // EOS value plus 2 derivatives values for each
	    // (density,temperature) tuple.
	    const int returnSize = 3 * vdensity.size(); 

	    std::vector< double > returnVals( returnSize );

 	    int errorCode 
		= wrapper::es1vals( returnType, derivatives,
				    interpolation, eosTable,
				    eosTableLength,
				    vdensity, vtempsKelvin,
				    returnVals, returnSize );
	    
	    if ( errorCode != 0 )
		{
		    std::ostringstream outputString;
		    outputString << "\n\tAn unsuccessful request for EOSPAC data "
				 << "was made by es1vals() from within getdFdT().\n"
				 << "\tThe error code returned was \""
				 << errorCode << "\".\n"
				 << "\tThe associated error message is:\n\t\""
				 << wrapper::es1errmsg( errorCode ) 
				 << "\"\n";

		    // This is a fatal exception right now.  It might
		    // be useful to throw a specific exception that is 
		    // derived from EospacException.  The host code
		    // could theoretically catch such an exception,
		    // fix the problem and then continue.
		    throw EospacException( outputString.str() );
		}

 	    // copy the results back into a STL vector
 	    std::vector< double > dFdT( vtemperature.size() );
	    
	    // dF/dT values are the last 1/3 of returnVals.
	    std::copy( 
		// returnVals+2*vtemperature.size(),
		returnVals.begin()+2*vtemperature.size(),
		       returnVals.end(),
		       dFdT.begin() );
	    
	    // the values in "returnVals" are actually log(vals) so
	    // that dF/dT = returnVal[2]/T, dF/drho =
	    // returnVal[1]/rho.  returnVal[0] is the EoS value (not a 
	    // derivative value.
	    for ( size_t i=0; i<vtemperature.size(); ++i )
// Return units for temperature should be keV not Kelvin.
// Heat Capacity will be kJ/g/keV!
// 		dFdT[i] = dFdT[i]/vtempsKelvin[i];
		dFdT[i] = dFdT[i]/vtemperature[i];

	    return dFdT;
	}

    /*!
     * \brief This member function examines the contents of
     *        the data member "SesTabs" and then calls the EOSPAC
     *        routine to load the required EoS Tables.
     */
    void Eospac::expandEosTable() const
	{
	    // loop over all possible EOSPAC data types.  If a matid
	    // has been assigned to a table then add this information
	    // to the vectors returnTypes[] and matIDs[] which are
	    // used by EOSPAC.	  

	    for ( int i=0; i<SesTabs.getNumReturnTypes(); ++i )
		if ( SesTabs.returnTypes( i ) != ES4null )
		    {
			// MatIDs[] and returnTypes[] are a tuple.

			returnTypes.insert( returnTypes.begin(),
					    SesTabs.returnTypes( i ) ); 
			matIDs.insert( matIDs.begin(),
				       SesTabs.matID(
					   SesTabs.returnTypes( i ) ) );
		    }

	    // Allocate eosTable.  The length and location of eosTable 
	    // will be modified by es1tabs() as needed.	  
	    eosTable = new V_FLOAT [ eosTableLength ];

	    // Initialize eosTable and find it's required length
	    int errorCode 
		= wrapper::es1tabs( numRegions, returnTypes.size(),
				    returnTypes, matIDs,
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

		    // Clean up temporaries before we throw the exception.
		    delete [] eosTable;

		    // This is a fatal exception right now.  It might
		    // be useful to throw a specific exception that is 
		    // derived from EospacException.  The host code
		    // could theoretically catch such an exception,
		    // fix the problem and then continue.
		    throw EospacException( outputString.str() );
		}
	    
	    // We don't delete eosTable until ~Eospac() is called.
	}

    /*!
     * \brief Returns true if the EoS data associated with
     *        "returnType" has been loaded.
     */
    bool Eospac::typeFound( ES4DataType returnType ) const
	{
	    // Loop over all available types.  If the requested
	    // type id matches on in the list then return true.
	    // If we reach the end of the list without a match return
	    // false. 
	    
	    for ( size_t i=0; i<returnTypes.size(); ++i )
		if ( returnType == returnTypes[i] ) return true;
	    return false;
	}

    
    /*!
     * \brief Converts a double to a length one vector.
     */
    std::vector< double > Eospac::dbl_v1( const double dbl ) const
	{
	    std::vector< double > vec(1);
	    vec[0] = dbl;
	    return vec;	    
	}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
// end of Eospac.cc
//---------------------------------------------------------------------------//
