//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/IpcressGrayOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 22 14:11:10 2001
 * \brief  IpcressGrayOpacity templated class implementation file.
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IpcressGrayOpacity.hh"
#include "IpcressFile.hh"       // we have smart pointers to
                                // IpcressFile objects.
#include "IpcressDataTable.hh"  // we have a smart pointer to a
                                // IpcressDataTable object.
#include "ds++/Assert.hh"       // we make use of Require()
#include "ds++/Packing_Utils.hh"
#include <cmath>                // we need to define log(double) and exp(double)

namespace rtt_cdi_ipcress
{

// ------------ //
// Constructors //
// ------------ //
    
/*!
 * \brief Constructor for IpcressGrayOpacity object.
 * 
 * See IpcressGrayOpacity.hh for details.
 */
IpcressGrayOpacity::IpcressGrayOpacity( 
    rtt_dsxx::SP< IpcressFile const > const & in_spIpcressFile,
    size_t            in_materialID,
    rtt_cdi::Model    in_opacityModel,
    rtt_cdi::Reaction in_opacityReaction )
    : spIpcressFile( in_spIpcressFile ),
      materialID(    in_materialID ),
      numKeys(       0 ),
      vKnownKeys(),
      opacityModel(    in_opacityModel ),
      opacityReaction( in_opacityReaction ),
      energyPolicyDescriptor( "gray" ),
      spIpcressDataTable()
{
    // Verify that the requested material ID is available in the
    // specified IPCRESS file.
    Insist( spIpcressFile->materialFound( materialID ),
            "The requested material ID is not available in the specified Ipcress file."); 
	    
    // Retrieve keys available for this material from the IPCRESS
    // file.  wgkeys() returns vKnownKeys, numKeys and errorCode.
    int errorCode = -99;

    // spIpcressFile->getKeys();
    
    vKnownKeys = spIpcressFile->listDataFieldNames( materialID );
    for( size_t i=0; i<vKnownKeys.size(); ++i )
        std::cout << "vKnownKeys[" << i << "] = " << vKnownKeys[i] << std::endl;

    // // find 'comp' entry for materialID
    // std::string keyword_comp("comp");
    // size_t idx=0;
    // for( size_t i=0; i<vKnownKeys.size(); i+=3 )
    //     if( atoi(vKnownKeys[i].c_str()) == static_cast<int>(materialID) &&
    //         keyword_comp == vKnownKeys[i+1].substr(0,keyword_comp.size()) )
    //             idx=i;
    // Check(idx>0);
    // std::cout << "idx = " << idx << std::endl; /// idx = 60
    // size_t comp_size = spIpcressFile->getArraySize(idx/3);
    // size_t comp_addr = spIpcressFile->getArrayAddr(idx/3);


    
    
    //     wrapper::wgkeys( spIpcressFile->getDataFilename(),
    //     		 materialID, vKnownKeys,
    //     		 wrapper::maxKeys, numKeys );
    // if ( errorCode !=0 ) throw gkeysException( errorCode );
    Check(errorCode>0);
	    
    // Create the data table object and fill it with the table
    // data from the IPCRESS file.
    spIpcressDataTable = new IpcressDataTable(
	energyPolicyDescriptor,
	opacityModel, 
	opacityReaction,
	vKnownKeys,
	materialID, 
	spIpcressFile );
	    
} // end of IpcressData constructor
    
/*!
 * \brief Unpacking constructor for IpcressGrayOpacity object.
 * 
 * See IpcressGrayOpacity.hh for details.
 */
IpcressGrayOpacity::IpcressGrayOpacity( std::vector<char> const & packed )
    : spIpcressFile(),
      materialID( 0 ),
      numKeys( 0 ),
      vKnownKeys(),
      opacityModel(),
      opacityReaction(),
      energyPolicyDescriptor( "gray" ),
      spIpcressDataTable()
{
    Require (packed.size() >= 5 * sizeof(int));

    // make an unpacker
    rtt_dsxx::Unpacker unpacker;
    unpacker.set_buffer(packed.size(), &packed[0]);

    // unpack and check the descriptor
    int         packed_descriptor_size = 0;
    unpacker >> packed_descriptor_size;
    Check (packed_descriptor_size > 0);

    // make a vector<char> for the packed descriptor
    std::vector<char> packed_descriptor(packed_descriptor_size);

    // unpack it
    std::string descriptor;
    for (int i = 0; i < packed_descriptor_size; ++i)
	unpacker >> packed_descriptor[i];
    rtt_dsxx::unpack_data(descriptor, packed_descriptor);

    // make sure it is "gray"
    Insist (descriptor == "gray", 
	    "Tried to unpack a non-gray opacity in IpcressGrayOpacity.");

    // unpack the size of the packed filename
    int packed_filename_size(0);
    unpacker >> packed_filename_size;

    // make a vector<char> for the packed filename
    std::vector<char> packed_filename(packed_filename_size);

    // unpack it
    std::string filename;
    for (int i = 0; i < packed_filename_size; ++i)
	unpacker >> packed_filename[i];
    rtt_dsxx::unpack_data(filename, packed_filename);

    // unpack the material id
    // unpacker >> materialID;
    int itmp(0);
    unpacker >> itmp;
    materialID = static_cast<size_t>(itmp);

    // unpack the model and reaction
    int model    = 0;
    int reaction = 0;
    unpacker >> model >> reaction;

    opacityModel    = static_cast<rtt_cdi::Model>(model);
    opacityReaction = static_cast<rtt_cdi::Reaction>(reaction);

    Ensure (unpacker.get_ptr() == &packed[0] + packed.size());
    
    // build a new IpcressFile
    spIpcressFile = new IpcressFile(filename);
    Check (spIpcressFile);

    // Verify that the requested material ID is available in the
    // specified IPCRESS file.
    Insist( spIpcressFile->materialFound( materialID ),
        "Requested material ID is not found in the specified Ipcress file.");
	    
    // Retrieve keys available for this material from the IPCRESS
    // file.  wgkeys() returns vKnownKeys, numKeys and errorCode.
    int errorCode = -99;
    //     wrapper::wgkeys( spIpcressFile->getDataFilename(),
    //     		 materialID, vKnownKeys,
    //     		 wrapper::maxKeys, numKeys );
    // if ( errorCode !=0 ) throw gkeysException( errorCode );
    Check(errorCode>0);
	    
    // Create the data table object and fill it with the table
    // data from the IPCRESS file.
    spIpcressDataTable = new IpcressDataTable(
	energyPolicyDescriptor,
	opacityModel, 
	opacityReaction,
	vKnownKeys,
	materialID, 
	spIpcressFile );

    Ensure (spIpcressFile);
    Ensure (spIpcressDataTable);
} 

/*!
 * \brief Desctructor for IpcressGrayOpacity class.
 */ 
IpcressGrayOpacity::~IpcressGrayOpacity()
{
    // empty
}
    
// --------- //
// Accessors //
// --------- //
    
/*!
 * \brief Returns a "plain English" description of the opacity
 *     data that this class references. (e.g. "Gray Rosseland
 *     Scattering".) 
 *
 * The definition of this function is not included here to prevent 
 *     the inclusion of the IpcressFile.hh definitions within this 
 *     header file.
 */
std::string IpcressGrayOpacity::getDataDescriptor() const 
{
    // call the correct function from the IpcressDataTable
    // object. 
    return spIpcressDataTable->getDataDescriptor(); 
}
    
/*!
 * \brief Returns the name of the associated IPCRESS file.
 *
 *     The definition of this function is not included here to
 *     prevent the inclusion of the IpcressFile.hh definitions
 *     within this header file.
 */
std::string IpcressGrayOpacity::getDataFilename() const 
{ 
    return spIpcressFile->getDataFilename(); 
}
 
/*!
 * \brief Opacity accessor that returns a single opacity (or a
 *     vector of opacities for the multigroup EnergyPolicy) that 
 *     corresponds to the provided temperature and density.
 */
double IpcressGrayOpacity::getOpacity(
    double targetTemperature,
    double targetDensity ) const
{ 
    double opacity=-99.0*targetTemperature*targetDensity;
    // logorithmic interpolation:
    // opacity =
    //     wrapper::wgintgrlog( spIpcressDataTable->getLogTemperatures(),
    //     		     spIpcressDataTable->getNumTemperatures(), 
    //     		     spIpcressDataTable->getLogDensities(), 
    //     		     spIpcressDataTable->getNumDensities(),
    //     		     spIpcressDataTable->getLogOpacities(), 
    //     		     spIpcressDataTable->getNumOpacities(),
    //     		     std::log(targetTemperature),
    //     		     std::log(targetDensity) );
    Check(opacity>0.0);
    return opacity;
}
    
/*!
 * \brief Opacity accessor that returns a vector of opacities (or a
 *     vector of vectors of opacities for the multigroup
 *     EnergyPolicy) that correspond to the provided vector of
 *     temperatures and a single density value.
 */
std::vector< double > IpcressGrayOpacity::getOpacity(
    std::vector<double> const & targetTemperature,
    double                      targetDensity ) const
{ 
    std::vector< double > opacity( targetTemperature.size(),-99.0*targetTemperature[0]*targetDensity );
    for ( size_t i=0; i<targetTemperature.size(); ++i )
	// logorithmic interpolation:
	// opacity[i] = 
	//     wrapper::wgintgrlog( spIpcressDataTable->getLogTemperatures(),
	// 			 spIpcressDataTable->getNumTemperatures(), 
	// 			 spIpcressDataTable->getLogDensities(), 
	// 			 spIpcressDataTable->getNumDensities(),
	// 			 spIpcressDataTable->getLogOpacities(), 
	// 			 spIpcressDataTable->getNumOpacities(),
	// 			 std::log(targetTemperature[i]),
	// 			 std::log(targetDensity) );
        Check(false);
    return opacity;
}

/*!
 * \brief Opacity accessor that returns a vector of opacities (or a
 *     vector of vectors of opacities for the multigroup
 *     EnergyPolicy) that correspond to the provided vector of
 *     densities and a single temperature value.
 */
std::vector< double > IpcressGrayOpacity::getOpacity(
    double                      targetTemperature,
    std::vector<double> const & targetDensity ) const
{ 
    std::vector< double > opacity( targetDensity.size(), -99*targetTemperature );
    for ( size_t i=0; i<targetDensity.size(); ++i )
	// logorithmic interpolation:
	// opacity[i] = 
	//     wrapper::wgintgrlog( spIpcressDataTable->getLogTemperatures(),
	// 			 spIpcressDataTable->getNumTemperatures(), 
	// 			 spIpcressDataTable->getLogDensities(), 
	// 			 spIpcressDataTable->getNumDensities(),
	// 			 spIpcressDataTable->getLogOpacities(), 
	// 			 spIpcressDataTable->getNumOpacities(),
	// 			 std::log(targetTemperature),
	// 			 std::log(targetDensity[i]) );
        Check(false);
    return opacity;
}
    
/*!
 * \brief Returns a vector of temperatures that define the cached
 *     opacity data table.
 */
std::vector< double > IpcressGrayOpacity::getTemperatureGrid() const
{
    return spIpcressDataTable->getTemperatures();
}
    
/*!
 * \brief Returns the size of the temperature grid.
 */
size_t IpcressGrayOpacity::getNumTemperatures() const
{
    return spIpcressDataTable->getNumTemperatures();
}
    
/*!
 * \brief Returns a vector of densities that define the cached opacity data
 *     table.
 */
std::vector<double> IpcressGrayOpacity::getDensityGrid() const
{
    return spIpcressDataTable->getDensities();
}
    
//!  Returns the size of the density grid.
size_t IpcressGrayOpacity::getNumDensities() const
{
    return spIpcressDataTable->getNumDensities();
}
     
// ------- //
// Packing //
// ------- //

/*!
 * Pack the IpcressGrayOpacity state into a char string represented by a \c
 * vector<char>. This can be used for persistence, communication, etc. by
 * accessing the \c char* under the vector (required by implication by the
 * standard) with the syntax \c &char_string[0]. Note, it is unsafe to use
 * iterators because they are \b not required to be \c char*.
 */
std::vector<char> IpcressGrayOpacity::pack() const
{
    using std::vector;
    using std::string;

    Require (spIpcressFile);

    // pack up the energy policy descriptor
    vector<char> packed_descriptor;
    rtt_dsxx::pack_data(energyPolicyDescriptor, packed_descriptor);

    // pack up the ipcress file name
    string       filename = spIpcressFile->getDataFilename();
    vector<char> packed_filename;
    rtt_dsxx::pack_data(filename, packed_filename);

    // determine the total size: 3 ints (reaction, model, material id) + 2
    // ints for packed_filename size and packed_descriptor size + char in
    // packed_filename and packed_descriptor
    size_t size = 5 * sizeof(int) + packed_filename.size() + 
                  packed_descriptor.size();

    // make a container to hold packed data
    vector<char> packed(size);

    // make a packer and set it
    rtt_dsxx::Packer packer;
    packer.set_buffer(size, &packed[0]);

    // pack the descriptor
    packer << static_cast<int>(packed_descriptor.size());
    for (size_t i = 0; i < packed_descriptor.size(); i++)
	packer << packed_descriptor[i];

    // pack the filename (size and elements)
    packer << static_cast<int>(packed_filename.size());
    for (size_t i = 0; i < packed_filename.size(); i++)
	packer << packed_filename[i];

    // pack the material id
    packer << static_cast<int>(materialID);

    // pack the model and reaction
    packer << static_cast<int>(opacityModel) 
	   << static_cast<int>(opacityReaction);

    Ensure (packer.get_ptr() == &packed[0] + size);
    return packed;
}

} // end namespace rtt_cdi_ipcress

//---------------------------------------------------------------------------//
// end of IpcressGrayOpacity.cc
//---------------------------------------------------------------------------//
