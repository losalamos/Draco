
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/gpu_hello_driver_api.cc
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 25 15:28:48 2011
 * \brief  Simple test of the CUDA Driver API.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../GPU_Device.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/path.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <vector>

//---------------------------------------------------------------------------//
// Helpers
//---------------------------------------------------------------------------//

int align(int offset, int alignment)
{
    // std::cout << "(offset,alignment) = " << offset << ", " << alignment
    //           << std::endl;
    // std::cout << "~(alignment-1) = " << ~(alignment - 1) << std::endl;
    // int foo = (offset + alignment - 1) & ~(alignment - 1);
    // std::cout << "result = " << foo << std::endl; 
    return (offset + alignment - 1) & ~(alignment - 1);
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

void simple_add( rtt_dsxx::ScalarUnitTest & ut )
{
    using namespace std;

    // Create a GPU_Device object.
    rtt_device::GPU_Device gpu;

    // Print a summary of the devices found.
    gpu.printDeviceSummary();

    // Only use device 0 (if there are more than 1 GPUs available)
    int iDev(0);

    // Create a context
    cout << "Get device handle" << endl;
    CUdevice device;
    cuDeviceGet(&device, iDev);

    cout << "Create a context" << endl;
    CUcontext ctx;
    cuCtxCreate(&ctx, iDev, device);

    // Load the module, must compile the kernel with nvcc -ptx -m32 kernel.cu
    string ptxFile("gpu_kernel.ptx");
    cout << "Load module " << ptxFile << endl;
    Insist( rtt_dsxx::fileExists( ptxFile ),
            (string("Cannot find requested file: kernel.ptx.")
             + ptxFile).c_str());
    CUmodule module;
    cuModuleLoad(&module, ptxFile.c_str());

    // Load the kernel from the module
    cout << "Load kernel \"sum\" from the module." << endl;
    CUfunction kernel;
    cuModuleGetFunction(&kernel,module,"sum");

    // Allocate some memory for the result
    cout << "Allocate memory on the device." << endl;
    CUdeviceptr dest;
    cuMemAlloc(&dest, sizeof(int));
	
    // Setup kernel parameters
    int offset(0);
    offset = align(offset, __alignof(CUdeviceptr));

    // cuParamSetv is used for pointers...
    cuParamSetv(kernel, offset, &dest, sizeof(CUdeviceptr));
    offset += sizeof(CUdeviceptr);

    offset = align(offset, __alignof(int));
    cuParamSeti(kernel, offset, 4);    // cuParamSeti is used for integers.
    offset += sizeof(int);
    offset = align(offset, __alignof(int));
    cuParamSeti(kernel, offset, 34);
    offset += sizeof(int);
    cuParamSetSize(kernel, offset);
    
    // Launch the grid
    cout << "Launch the grid" << endl;
    cuFuncSetBlockShape(kernel, 1, 1, 1);
    cuLaunchGrid(kernel, 1, 1);

    // Read the result off of the GPU
    cout << "Read the result" << endl;
    int result = 0;
    cuMemcpyDtoH(&result, dest, sizeof(int));

    cout << "Sum of 4 and 34 is " << result << endl;

    if( result == 38 )
        ut.passes("Sum of 4 and 34 is 38.");
    else
        ut.failure("Sum of 4 and 34 was incorrect.");

// deallocate memory, free the context.
    cout << "deallocate device memory." << endl;
    cuMemFree(dest);
    cout << "unload module" << endl;
    cuModuleUnload(module);
    cout <<  "destroy context" << endl;
    cuCtxDestroy(ctx);
    
    return;
}

//---------------------------------------------------------------------------//
// Main
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using namespace std;
    
    rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        simple_add(ut);
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing gpu_hello_driver_api, " << err.what()
             << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing gpu_hello_driver_api, " 
             << "An unknown exception was thrown." << endl;
        ut.numFails++;
    }
    return ut.numFails;
} 

//---------------------------------------------------------------------------//
// end of gpu_hello_driver_api.cc
//---------------------------------------------------------------------------//
