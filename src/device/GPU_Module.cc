//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/GPU_Module.cc
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 20 15:28:48 2011
 * \brief  
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "device/config.h"
#include "GPU_Module.hh"
#include "GPU_Device.hh"
#include "ds++/path.hh"
#include <sstream>

namespace rtt_device
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * Create a GPU_Module object.
 */
GPU_Module::GPU_Module( std::string const & myPtxFile ) 
    : ptxFile( findPtxFile( myPtxFile ) )
{
    // load the module
    cudaError_enum err = cuModuleLoad(&cuModule, ptxFile.c_str());
    rtt_device::GPU_Device::checkForCudaError( err );
}

//---------------------------------------------------------------------------//
/*!
 * \brief destructor
 *
 * Free the the loaded modules.
 */
GPU_Module::~GPU_Module()
{
    // Unload the module
    cudaError_enum err = cuModuleUnload(cuModule);
    rtt_device::GPU_Device::checkForCudaError( err );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief findPtxFile
 * 
 * \param myPtxFile filename or path to the ptx file that is to be loaded.
 * \return fully qualified path to the disired ptx file
 */
std::string GPU_Module::findPtxFile( std::string const & myPtxFile )
{
    // Location of GPU ptx files - read from config.h
    std::string const testDir( rtt_device::test_ppe_bindir ); 
    // return value
    std::string ptxFile;
    
    // Find the ptx file
    if( rtt_dsxx::fileExists( myPtxFile ) )
        ptxFile = myPtxFile;
    else if( rtt_dsxx::fileExists( testDir + std::string("/") + myPtxFile ) )
        ptxFile = testDir + std::string("/") + myPtxFile;
    else
        Insist( rtt_dsxx::fileExists( ptxFile ),
                ( std::string("Cannot find requested file: ")
                  + ptxFile).c_str() );
    return ptxFile;
}

} // end namespace rtt_device

//---------------------------------------------------------------------------//
// end of GPU_Module.cc
//---------------------------------------------------------------------------//
