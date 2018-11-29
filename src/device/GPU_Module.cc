//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/GPU_Module.cc
 * \brief  
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "GPU_Module.hh"
#include "GPU_Device.hh"
#include "device/config.h"
#include "ds++/path.hh"
#include <sstream>

namespace rtt_device {

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 * \arg myPtxFile - the name of a ptx or cubin file.  Ptx files will be
 * compiled at runtime into cubins.  Prefer the use of cubin to ruduce runtime.
 *
 * Create a GPU_Module object.
 */
GPU_Module::GPU_Module(std::string const &myPtxFile)
    : ptxFile(findPtxFile(myPtxFile)) {
  // load the module
  cudaError_enum err = cuModuleLoad(&cuModule, ptxFile.c_str());
  rtt_device::GPU_Device::checkForCudaError(err);
}

//---------------------------------------------------------------------------//
/*!
 * \brief destructor
 *
 * Free the the loaded modules.
 */
GPU_Module::~GPU_Module() {
  // Unload the module
  cudaError_enum err = cuModuleUnload(cuModule);
  rtt_device::GPU_Device::checkForCudaError(err);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief findPtxFile
 * 
 * \param myPtxFile filename or path to the ptx file that is to be loaded.
 * \return fully qualified path to the disired ptx file
 */
std::string GPU_Module::findPtxFile(std::string const &myPtxFile) {
  Require(myPtxFile.length() > 0);

  // Location of GPU ptx files - read from config.h
  std::string const testDir(rtt_device::test_kernel_bindir);
  // return value
  std::string ptxFile;

  // std::cout << "Looking at:\n"
  //           << myPtxFile << "\n"
  //           << testDir + std::string("/") + myPtxFile << std::endl;

  // Find the ptx file
  if (rtt_dsxx::fileExists(myPtxFile))
    ptxFile = myPtxFile;
  else if (rtt_dsxx::fileExists(std::string("../") + myPtxFile))
    ptxFile = std::string("../") + myPtxFile;
  else if (rtt_dsxx::fileExists(testDir + std::string("/") + myPtxFile))
    ptxFile = testDir + std::string("/") + myPtxFile;

  Insist(rtt_dsxx::fileExists(ptxFile),
         (std::string("Cannot find requested file: ") + myPtxFile).c_str());

  return ptxFile;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Find a function in the current module and return a handle.
 * 
 * \param functionName the name of the CUDA function difined in the ptx
 *                     module. 
 * \return a CUfunction handle that points to the requested function.
 */
CUfunction
GPU_Module::getModuleFunction(std::string const &functionName) const {
  // Load the kernel from the module
  CUfunction kernel;
  cudaError_enum err =
      cuModuleGetFunction(&kernel, cuModule, functionName.c_str());
  GPU_Device::checkForCudaError(err);
  return kernel;
}

} // end namespace rtt_device

//---------------------------------------------------------------------------//
// end of GPU_Module.cc
//---------------------------------------------------------------------------//
