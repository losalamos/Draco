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

#include "Eospac.hh"
#include "EospacException.hh"

// Other Draco dependencies
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"

// C++ standard library dependencies
#include <string>
#include <sstream>
#include <iostream>

namespace rtt_cdi_eospac
{

// ------------ //
// Constructors //
// ------------ //
    
//---------------------------------------------------------------------------//
/*!
 * \brief The constructor for Eospac.
 *
 * \sa The definition of rtt_cdi_eospac::SesameTables.
 *
 * \param SesTabs A rtt_cdi_eospac::SesameTables object that defines what data
 * tables will be available for queries from the Eospac object.
 *
 */
Eospac::Eospac( SesameTables const & in_SesTabs )
    : SesTabs( in_SesTabs ),
      infoItems( initializeInfoItems() ),
      infoItemDescriptions( initializeInfoItemDescriptions() )
{
    // Eospac can only be instantiated if SesameTables is provided.  If
    // SesameTables is invalid this will be caught in expandEosTable();

    // PreCache the default data type
    expandEosTable();

} // end Eospac::Eospac()

//---------------------------------------------------------------------------//
/*! 
 * \brief Construct an Eospac by unpacking a vector<char> stream
 */
Eospac::Eospac( std::vector<char> const & packed )
    : SesTabs( SesameTables( packed ) ),
      infoItems( initializeInfoItems() ),
      infoItemDescriptions( initializeInfoItemDescriptions() )
{
    // Use the initializer list to build the needed SesTabs.  Now initialize
    // libeospac using the data found in SesTabs just like in the default
    // ctor.     
    expandEosTable();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Default Eospac() destructor.
 *
 * This is required to correctly release memeroy when an Eospac object is
 * destroyed.  We define the destructor in the implementation file to avoid
 * including the unnecessary header files.
 */
Eospac::~Eospac()
{
    // Destroy all data objects:
    EOS_INTEGER errorCode( EOS_OK );
    eos_DestroyAll( &errorCode );
    if( errorCode != EOS_OK )
    {
        std::ostringstream outputString;
        for( size_t i=0; i<returnTypes.size(); ++i )
        {
            EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
            EOS_INTEGER tableHandleErrorCode = EOS_OK;
            eos_GetErrorCode( &tableHandles[i], &tableHandleErrorCode );
            eos_GetErrorMessage( &tableHandleErrorCode, errorMessage);
            
            outputString
                << "\n\tAn unsuccessful request was made to destroy the "
                << "EOSPAC table area by ~Eospac().\n"
                << "\tThe error code returned by eos_DestroyAll(...) was \""
                << tableHandleErrorCode << "\".\n"
                << "\tThe associated error message is:\n\t\""
                << errorMessage << "\"\n";
        }
        // Never throw an exception from the destructor.  This can cause
        // confusion during stack unwiding.
        // throw EospacException( outputString.str() );
        std::cerr << outputString.str() << std::endl;
    }
}
    
// --------- //
// Accessors //
// --------- //

void Eospac::printTableInformation( EOS_INTEGER const tableType,
                                    std::ostream & out ) const
{
    // Obtain the table handle for this type
    EOS_INTEGER tableHandle( tableHandles[ tableIndex( tableType )]);
    
    EOS_INTEGER one(1);
    EOS_INTEGER errorCode(EOS_OK);
    // std::vector< EOS_REAL > infoVals( infoItems.size() );
    EOS_REAL infoVal;

    out << "\nEOSPAC information for Table " << SesTabs.tableName[tableType]
        << " (" << SesTabs.tableDescription[tableType] << ")\n"
        << "-------------------------------------------------------------"
        << "-------------------------\n";
    
    for( size_t i=0; i<infoItems.size(); ++i )
    {
        eos_GetTableInfo( &tableHandle, &one, &infoItems[i], &infoVal, 
                          &errorCode );
        if( errorCode == EOS_OK )
        {
            out << std::setiosflags(std::ios::fixed)
                << std::setw(70) << std::left << infoItemDescriptions[i] << ": "
                << std::setprecision(6) << std::setiosflags(std::ios::fixed)
                << std::setw(13) << std::right << infoVal << std::endl;
        }
        else if( errorCode != EOS_INVALID_INFO_FLAG )
        {
            std::ostringstream outputString;
            EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
            // Ignore EOS_INVALID_INFO_FLAG since not all infoItems are
            // currently applicable to a specific tableHandle. 
            eos_GetErrorMessage( &errorCode, errorMessage );
            outputString << "\n\tAn unsuccessful request for EOSPAC table information "
                     << "was made by eos_GetTableInfo().\n"
                     << "\tThe requested infoType was \"" 
                     << infoItems[i] << "\" (see eos_Interface.h for type)\n"
                     << "\tThe error code returned was \""
                     << errorCode << "\".\n"
                     << "\tThe associated error message is:\n\t\""
                     << errorMessage
                     << "\"\n";
            throw EospacException( outputString.str() );
        }
    }
    return;
}

//---------------------------------------------------------------------------//
double Eospac::getSpecificElectronInternalEnergy(
    double temperature, double density ) const
{
    EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc;
    return getF( dbl_v1(temperature), dbl_v1(density),
                 returnType, ETDD_VALUE )[0];
}

//---------------------------------------------------------------------------//
std::vector< double > Eospac::getSpecificElectronInternalEnergy(
    std::vector< double > const & vtemperature,
    std::vector< double > const & vdensity ) const
{
    EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc; 
    return getF( vtemperature, vdensity, returnType, ETDD_VALUE );
}

//---------------------------------------------------------------------------//
double Eospac::getElectronHeatCapacity(
    double temperature, double density ) const
{
    // specific Heat capacity is dE/dT at constant pressure.  To obtain the
    // specific electron heat capacity we load the specific electron internal
    // energy (E) and it's first derivative w.r.t temperature.

    EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc
    return getF( dbl_v1(temperature), dbl_v1(density),
                 returnType, ETDD_DFDY )[0];
}

//---------------------------------------------------------------------------//
std::vector< double > Eospac::getElectronHeatCapacity(
    std::vector< double > const & vtemperature, 
    std::vector< double > const & vdensity ) const
{
    // specific Heat capacity is dE/dT at constant pressure.  To obtain the
    // specific electron heat capacity we load the specific electron internal
    // energy (E) and it's first derivative w.r.t temperature.

    EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc;
    return getF( vtemperature, vdensity, returnType, ETDD_DFDY );
}

//---------------------------------------------------------------------------//
double Eospac::getSpecificIonInternalEnergy(
    double temperature, double density ) const
{
    EOS_INTEGER const returnType = EOS_Uic_DT; // ES4enion;
    return getF( dbl_v1(temperature), dbl_v1(density),
                 returnType, ETDD_VALUE )[0];
}

//---------------------------------------------------------------------------//
std::vector< double > Eospac::getSpecificIonInternalEnergy(
    std::vector< double > const & vtemperature, 
    std::vector< double > const & vdensity ) const
{
    EOS_INTEGER const returnType = EOS_Uic_DT; //ES4enion;
    return getF( vtemperature, vdensity, returnType, ETDD_VALUE );
}

//---------------------------------------------------------------------------//
double Eospac::getIonHeatCapacity(
    double temperature, double density ) const
{
    // specific Heat capacity is dE/dT at constant pressure.  To obtain the
    // specific electron heat capacity we load the specific electron internal
    // energy (E) and it's first derivative w.r.t temperature.

    EOS_INTEGER const returnType = EOS_Uic_DT; //ES4enion; 
    return getF( dbl_v1(temperature), dbl_v1(density), 
                 returnType, ETDD_DFDY  )[0];
}

//---------------------------------------------------------------------------//
std::vector< double > Eospac::getIonHeatCapacity(
    std::vector< double > const & vtemperature,
    std::vector< double > const & vdensity ) const
{
    // specific Heat capacity is dE/dT at constant pressure.  To obtain the
    // specific electron heat capacity we load the specific electron internal
    // energy (E) and it's first derivative w.r.t temperature.

    EOS_INTEGER const returnType = EOS_Uic_DT; //ES4enion;
    return getF( vtemperature, vdensity, returnType, ETDD_DFDY  );
}
    
//---------------------------------------------------------------------------//
double Eospac::getNumFreeElectronsPerIon(
    double temperature, double density ) const
{
    EOS_INTEGER const returnType = EOS_Zfc_DT; // ES4zfree3; // (zfree3)
    return getF( dbl_v1(temperature), dbl_v1(density), 
                 returnType, ETDD_VALUE )[0];
}
    
//---------------------------------------------------------------------------//
std::vector< double > Eospac::getNumFreeElectronsPerIon(
    std::vector< double > const & vtemperature,
    std::vector< double > const & vdensity ) const
{
    EOS_INTEGER const returnType = EOS_Zfc_DT; //ES4zfree3; // (zfree3)
    return getF( vtemperature, vdensity, returnType, ETDD_VALUE );
}
    
//---------------------------------------------------------------------------//
double Eospac::getElectronThermalConductivity(
    double temperature, double density ) const
{
    EOS_INTEGER const returnType = EOS_Ktc_DT; //ES4tconde; // (tconde)
    return getF( dbl_v1(temperature), dbl_v1(density), 
                 returnType, ETDD_VALUE )[0];
}
    
//---------------------------------------------------------------------------//
std::vector< double > Eospac::getElectronThermalConductivity(
    std::vector< double > const & vtemperature, 
    std::vector< double > const & vdensity ) const
{
    EOS_INTEGER const returnType = EOS_Ktc_DT; //ES4tconde; // (tconde)
    return getF( vtemperature, vdensity, returnType, ETDD_VALUE );
}

// ------- //
// Packing //
// ------- //

//---------------------------------------------------------------------------//
/*!
 * Pack the Eospac state into a char string represented by a
 * vector<char>. This can be used for persistence, communication, etc. by
 * accessing the char * under the vector (required by implication by the
 * standard) with the syntax &char_string[0]. Note, it is unsafe to use
 * iterators because they are \b not required to be char *.
 *
 * Eospac has no state of its own to pack.  However, this function does call
 * the SesameTable pack() because this data is required to rebuild Eospac.
 */
std::vector<char> Eospac::pack() const
{
    return SesTabs.pack();
}

// -------------- //
// Implementation //
// -------------- //

//---------------------------------------------------------------------------//
/*!
 * \brief Retrieves the EoS data associated with the returnType specified and
 *        the given (density, temperature) tuples.
 *
 * Each of the public access functions calls either getF() or getdFdT() after
 * assigning the correct value to "returnType".
 *
 * \param vdensity A vector of density values (g/cm^3).
 * \param vtemperature A vector of temperature values (keV).
 * \param returnType The integer index that corresponds to the type of data
 *                 being retrieved from the EoS tables. 
 */
std::vector< double > Eospac::getF( 
    std::vector< double > const & vtemperature,
    std::vector< double > const & vdensity,
    EOS_INTEGER const returnType,
    EosTableDataDerivative const etdd ) const
{
    // The density and vector parameters must be a tuple.
    Require( vdensity.size() == vtemperature.size() );

    unsigned returnTypeTableIndex( tableIndex( returnType ));

    // Convert temperatures from keV to degrees Kelvin.
    std::vector< double > vtempsKelvin = vtemperature;
    std::transform( vtemperature.begin(), vtemperature.end(),
                    vtempsKelvin.begin(), keV2K );                    
	    
    // There is one piece of returned information for each (density,
    // temperature) tuple.
    int returnSize = vdensity.size();

    std::vector< double > returnVals( returnSize );
    std::vector< double > dFx( returnSize );
    std::vector< double > dFy( returnSize );
    int errorCode =0;
    std::vector< double > nc_vdensity( vdensity );

    eos_Interpolate( &tableHandles[returnTypeTableIndex],
                     &returnSize, &nc_vdensity[0], &vtempsKelvin[0],
                     &returnVals[0], &dFx[0], &dFy[0], &errorCode );
    	    
    if ( errorCode != 0 )
    {
        std::ostringstream outputString;
        EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
        eos_GetErrorMessage( &errorCode, errorMessage);

        outputString << "\n\tAn unsuccessful request for EOSPAC data "
                     << "was made by eos_Interpolate() from within getF().\n"
                     << "\tThe requested returnType was \"" 
                     << returnType << "\" (see eos_Interface.h for type)\n"
                     << "\tThe error code returned was \""
                     << errorCode << "\".\n"
                     << "\tThe associated error message is:\n\t\""
                     << errorMessage
                     << "\"\n";

    // If the EOS_INTERP_EXTRAPOLATED error code is returned by either
    // eos_Interpolate or eos_Mix, then the eos_CheckExtrap routine allows the
    // user to determine which (x,y) pairs caused extrapolation and in which
    // direction (high or low), it occurred. The units of the xVals, and yVals
    // arguments listed below are determined by the units listed for each
    // tableType in APPENDIX B and APPENDIX C.
    
    // eos_CheckExtrap( &tableHandles[returnTypeTableIndex],
    //                  &returnSize, &nc_vdensity[0], &vtempsKelvin[0],
    //                  &xyBounds, &errorCode );

        
        // This is a fatal exception right now.  It might be useful to throw a
        // specific exception that is derived from EospacException.  The host
        // code could theoretically catch such an exception, fix the problem
        // and then continue.
        throw EospacException( outputString.str() );
    }

    switch( etdd )
    {
        case ETDD_VALUE: return returnVals; break;
        case ETDD_DFDX:
        {
            return dFx;
            break;
        }
        case ETDD_DFDY:
        {
            // Convert back to temperature units of keV (kJ/g/K -> kJ/g/keV):
            std::transform( dFy.begin(), dFy.end(), dFy.begin(), keV2K );
            return dFy;
            break;
        }
        default:
            Insist( etdd == ETDD_VALUE || etdd == ETDD_DFDX || etdd == ETDD_DFDY,
                    "Bad value for EosTableDataDerivative." );
    }
    // should never get here.
    return returnVals;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This member function examines the contents of the data member
 *        "SesTabs" and then calls the EOSPAC routine to load the required EoS
 *        Tables.
 */
void Eospac::expandEosTable() const
{
    // loop over all possible EOSPAC data types.  If a matid has been assigned
    // to a table then add this information to the vectors returnTypes[] and
    // matIDs[] which are used by EOSPAC.

    for( size_t i=0; i<SesTabs.getNumReturnTypes(); ++i )
    {
        if ( SesTabs.returnTypes( i ) != EOS_NullTable )
        {
            // MatIDs[] and returnTypes[] are a tuple.

            returnTypes.insert( returnTypes.begin(),
                                SesTabs.returnTypes( i ) ); 
            matIDs.insert( matIDs.begin(),
                           SesTabs.matID(
                               SesTabs.returnTypes( i ) ) );
        }
    }

    // Allocate eosTable.  The length and location of eosTable will be
    // modified by es1tabs() as needed.
    for( size_t i=0; i<returnTypes.size(); ++i )
        tableHandles.push_back(EOS_NullTable);

    // Initialize eosTable and find it's required length

    int errorCode(0);
    int nTables( returnTypes.size() );
    eos_CreateTables( &nTables, &returnTypes[0], &matIDs[0], &tableHandles[0],
                      &errorCode );

    // Check for errors
    if ( errorCode != EOS_OK )
    {
        std::ostringstream outputString;
        for( size_t i=0; i<returnTypes.size(); ++i )
        {
            EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
            int tableHandleErrorCode( EOS_OK );
            eos_GetErrorCode(    &tableHandles[i], &tableHandleErrorCode );
            eos_GetErrorMessage( &tableHandles[i], errorMessage );
        
            outputString
                << "\n\tAn unsuccessful request was made to initialize the "
                << "EOSPAC table area by expandEosTable().\n"
                << "\tThe error code returned by eos_CreateTables(...) was \""
                << tableHandleErrorCode << "\".\n"
                << "\tThe associated error message is:\n\t\""
                << errorMessage << "\"\n";
        }
            
        // Clean up temporaries before we throw the exception.
        // delete [] eosTable;

        // This is a fatal exception right now.  It might be useful to throw a
        // specific exception that is derived from EospacException.  The host
        // code could theoretically catch such an exception, fix the problem
        // and then continue.
        throw EospacException( outputString.str() );
    }

    // Set options:
    // for( size_t i=0; i <nTables; i++)
    // {
    //     /* enable smoothing */
    //     eos_SetOption (&tableHandle[i], &EOS_SMOOTH, EOS_NullPtr, &errorCode);
    //     if (errorCode != EOS_OK) {
    //         eos_GetErrorMessage (&errorCode, errorMessage);
    //         cout << "eos_SetOption ERROR " << errorCode << ": " << errorMessage 
    //              << '\n';
    //     }
    // }
    
    // Load data into table data objects
    
    eos_LoadTables( &nTables, &tableHandles[0], &errorCode );

    if (errorCode != EOS_OK)
    {
        std::ostringstream outputString;
        for( size_t i=0; i<returnTypes.size(); ++i )
        {
            EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
            EOS_INTEGER tableHandleErrorCode( EOS_OK );
            eos_GetErrorCode(    &tableHandles[i], &tableHandleErrorCode );
            eos_GetErrorMessage( &tableHandles[i], errorMessage );
        
            outputString
                << "\n\tAn unsuccessful request was made to initialize the "
                << "EOSPAC table area by expandEosTable().\n"
                << "\tThe error code returned by eos_LoadTables(...) was \""
                << tableHandleErrorCode << "\".\n"
                << "\tThe associated error message is:\n\t\""
                << errorMessage << "\"\n";
        }
        throw EospacException( outputString.str() );        
    }
    
    // We don't delete eosTable until ~Eospac() is called.
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns true if the EoS data associated with "returnType" has been
 *        loaded.
 */
bool Eospac::typeFound( EOS_INTEGER returnType ) const
{
    // Loop over all available types.  If the requested type id matches on in
    // the list then return true.  If we reach the end of the list without a
    // match return false.
	    
    for ( size_t i=0; i<returnTypes.size(); ++i )
        if ( returnType == returnTypes[i] ) return true;
    return false;
}

//---------------------------------------------------------------------------//
unsigned Eospac::tableIndex( EOS_INTEGER returnType ) const
{
    // Loop over all available types.  If the requested type id matches on in
    // the list then return true.  If we reach the end of the list without a
    // match return false.

    // Throw an exception if the required return type has not been loaded by
    // Eospac.
    if ( ! typeFound( returnType ) ) 
    {
        std::ostringstream outputString;
        outputString << "\n\tA request was made for data by getF() "
                     << "for which EOSPAC does not have an\n"
                     << "\tassociated material identifier.\n"
                     << "\tRequested returnType = \"" 
                     << SesTabs.tableName[returnType]
                     << " (" << SesTabs.tableDescription[returnType] << ")\"\n";
        throw EospacUnknownDataType( outputString.str() );
    }
    
    for ( size_t i=0; i<returnTypes.size(); ++i )
        if ( returnType == returnTypes[i] )
            return i;
    
    return 0;
}

//---------------------------------------------------------------------------//
std::vector< EOS_INTEGER > Eospac::initializeInfoItems(void)
{
    std::vector< EOS_INTEGER > ii;
    ii.push_back( EOS_Cmnt_Len );
    ii.push_back( EOS_Exchange_Coeff );
    ii.push_back( EOS_F_Convert_Factor );
    ii.push_back( EOS_Log_Val );
    ii.push_back( EOS_Material_ID );
    ii.push_back( EOS_Mean_Atomic_Mass );
    ii.push_back( EOS_Mean_Atomic_Num );
    ii.push_back( EOS_Modulus );
    ii.push_back( EOS_Normal_Density );
    ii.push_back( EOS_Table_Type );
    ii.push_back( EOS_X_Convert_Factor );
    ii.push_back( EOS_Y_Convert_Factor );
    return ii;
}

//---------------------------------------------------------------------------//
std::vector< std::string > Eospac::initializeInfoItemDescriptions(void)
{
    using std::string;
    
    std::vector< std::string > iid;
     iid.push_back(
         string("The length in chars of the comments for the specified data table") );
     iid.push_back(
         string("The exchange coefficient") );
     iid.push_back(
         string("The conversion factor corresponding to the dependent variable, F(x,y)") );
     iid.push_back(
         string("Non-zero if the data table is in a log10 format") );
     iid.push_back(
         string("The SESAME material identification number") );
     iid.push_back(
         string("The mean atomic mass") );
     iid.push_back(
         string("The mean atomic number") );
     iid.push_back(
         string("The solid bulk modulus") );
     iid.push_back(
         string("The normal density") );
     iid.push_back(
         string("The type of data table. See APPENDIX B and APPENDIX C") );
     iid.push_back(
         string("The conv. factor corresponding to the primary indep. variable, x") );
     iid.push_back(
         string("The conv. factor corresponding to the secondary indep. variable, y") );
     return iid;
}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
// end of Eospac.cc
//---------------------------------------------------------------------------//
