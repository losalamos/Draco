//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/GPU_Module.hh
 * \author Kelly (KT) Thompson
 * \brief  Define class GPU_Module
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef device_GPU_Module_hh
#define device_GPU_Module_hh

#include "device_cuda.h"
#include <string>

namespace rtt_device {

//===========================================================================//
/*!
 * \class GPU_Module
 * \brief
 *
 * CUDA modules are dynamically loadable packages of device code and data
 * (akin to dynamic shared libraries)
 *
 * \sa GPU_Module.cc for detailed descriptions.
 *
 * \par Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 *
 * \sa device/test/gpu_hello_driver_api.cc
 */
//===========================================================================//

class GPU_Module {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Default constructors.
  GPU_Module(std::string const &myPtxFile);

  //! Copy constructor (the long doxygen description is in the .cc file).
  // GPU_Module(const GPU_Module &rhs);

  //! Destructor.
  ~GPU_Module();

  // MANIPULATORS

  //! Assignment operator for GPU_Module.
  // GPU_Module& operator=(const GPU_Module &rhs);

  // ACCESSORS
  CUmodule handle(void) { return cuModule; }

  // SERVICES
  CUfunction getModuleFunction(std::string const &functionName) const;

  // IMPLEMENTATION
  static std::string findPtxFile(std::string const &myPtxFile);

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  // DATA
  std::string const ptxFile;
  CUmodule cuModule;
};

} // end namespace rtt_device

#endif // device_GPU_Module_hh

//---------------------------------------------------------------------------//
// end of device/GPU_Module.hh
//---------------------------------------------------------------------------//
