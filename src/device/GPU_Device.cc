//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/GPU_Device.cc
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 20 15:28:48 2011
 * \brief  
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GPU_Device.hh"
#include "ds++/Assert.hh"

namespace rtt_device
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
GPU_Device::GPU_Device(void) // int /*argc*/, char */*argv*/[] )
    : deviceCount(0),
      computeCapability(),
      deviceName()
{
    // Initialize the library
    int err = cuInit(0); // currently must be 0.
    Check( err == CUDA_SUCCESS );
    
    // Get a device count, determine compute capability
    cuDeviceGetCount(&deviceCount);
    Insist( deviceCount>0, "No GPU devices found!");

    // Collect information about each GPU device found
    computeCapability.resize(deviceCount);

    for (int device = 0; device < deviceCount; device++)
    {
        CUdevice cuDevice;
        cuDeviceGet(&cuDevice, device);

        // Compute capability revision
        int major = 0;
        int minor = 0;
        cuDeviceComputeCapability(&major, &minor, cuDevice);
        computeCapability[device].push_back( major );
        computeCapability[device].push_back( minor );

        // Device name
        char name[200];
        cuDeviceGetName( name, 200, cuDevice );
        deviceName.push_back( std::string( name ) );

        // Query and archive device properties.
        CUdevprop_st properties;
        cuDeviceGetProperties( &properties, cuDevice );
        deviceProperties.push_back( properties );
    }
}

//---------------------------------------------------------------------------//
// Print a summary of all GPU devices found
//---------------------------------------------------------------------------//
void GPU_Device::printDeviceSummary(std::ostream & out) const
{
    using std::endl;
    
    out << "GPU device summary:\n\n"
         << "   Number of devices: " << deviceCount << "\n" << endl;
    
    for( int device=0; device<deviceCount; ++device )
    {
        out << "Device: " << device
            << "\n   Name               : " << deviceName[device]
            << "\n   Compute capability : " << computeCapability[device][0]
            << "." << computeCapability[device][1]
            << "\n   maxThreadsPerBlock : " << maxThreadsPerBlock(device)
            << "\n   maxThreadsDim      : " << maxThreadsDim(0,device) << " x "
            << maxThreadsDim(1,device) << " x " <<  maxThreadsDim(2,device)
            << "\n   maxGridSize        : " << maxGridSize(0,device) << " x "
            << maxGridSize(1,device) << " x " <<  maxGridSize(2,device)
            << "\n   sharedMemPerBlock  : " << sharedMemPerBlock(device)
            << "\n   totalConstantMemory: " << totalConstantMemory(device)
            << "\n   SIMDWidth          : " << SIMDWidth(device)
            << "\n   memPitch           : " << memPitch(device)
            << "\n   regsPerBlock       : " << regsPerBlock(device)
            << "\n   clockRate          : " << clockRate(device)
            << "\n   textureAlign       : " << textureAlign(device)
            << "\n" << endl;
    }
    return;
}


} // end namespace rtt_device

//---------------------------------------------------------------------------//
// end of GPU_Device.cc
//---------------------------------------------------------------------------//
