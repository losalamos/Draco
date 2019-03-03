//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/gpu_hello_driver_api.cc
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 25 15:28:48 2011
 * \brief  Simple test of the CUDA Driver API.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "device/GPU_Device.hh"
#include "device/GPU_Module.hh"
#include "device/config.h"
#include "ds++/Assert.hh"
#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/SystemCall.hh"
#include "ds++/path.hh"
#include <cstdlib> // RAND_MAX
#include <iostream>
#include <sstream>
#include <time.h>
#include <vector>

//---------------------------------------------------------------------------//
// Helpers
//---------------------------------------------------------------------------//

void genTestData(std::vector<double> &a, std::vector<double> &b,
                 std::vector<double> &ref) {
  // Initialize the random seed
  srand(time(NULL));

  // Fill arrays
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = static_cast<double>(rand() % 1000);
    b[i] = static_cast<double>(rand() % 1000);
    ref[i] = a[i] + b[i];
  }

  return;
}

//---------------------------------------------------------------------------//
// query_device
//---------------------------------------------------------------------------//

void query_device(rtt_dsxx::ScalarUnitTest &ut) {
  using namespace std;

  cout << "Starting gpu_hello_driver_api::query_device()...\n" << endl;

  // Create a GPU_Device object.
  // Initialize the CUDA library and sets device and context handles.
  rtt_device::GPU_Device gpu;

  // Create and then print a summary of the devices found.
  std::ostringstream out;
  size_t const numDev(gpu.numDevicesAvailable());
  out << "GPU device summary:\n\n"
      << "   Number of devices found: " << numDev << "\n"
      << endl;
  for (size_t device = 0; device < numDev; ++device)
    gpu.printDeviceSummary(device, out);

  // Print the message to stdout
  cout << out.str();

  // Parse the output
  bool verbose(false);
  std::map<std::string, unsigned> wordCount =
      rtt_dsxx::get_word_count(out, verbose);

  if (wordCount[string("Device")] == numDev)
    ut.passes("Found a report for each available device.");
  else
    ut.failure("Did not find a report for each available device.");

  return;
}

//---------------------------------------------------------------------------//
// Test: simple_add
//---------------------------------------------------------------------------//

void simple_add(rtt_dsxx::ScalarUnitTest &ut) {
  using namespace std;

  cout << "\nStarting gpu_hello_driver_api::simple_add()...\n" << endl;

  // Where are we?
  cout << "Paths:"
       << "\n   Current working dir = " << rtt_dsxx::draco_getcwd()
       << "\n   GPU kernel files at = " << rtt_device::test_kernel_bindir
       << endl;

  // Create a GPU_Device object.
  // Initialize the CUDA library and sets device and context handles.
  rtt_device::GPU_Device gpu;

  // Load the module, must compile the kernel with nvcc -ptx -m32 kernel.cu
  rtt_device::GPU_Module myModule("gpu_kernel.cubin");

  // Load the kernel from the module
  cout << "Load kernel \"sum\" from the module." << endl;
  CUfunction kernel;
  cudaError_enum err = cuModuleGetFunction(&kernel, myModule.handle(), "sum");
  gpu.checkForCudaError(err);

  // Allocate some memory for the result
  cout << "Allocate memory on the device." << endl;
  CUdeviceptr dest;
  err = cuMemAlloc(&dest, sizeof(int));
  gpu.checkForCudaError(err);

  // Setup kernel parameters
  int offset(0);
  offset = gpu.align(offset, __alignof(CUdeviceptr));

  // cuParamSetv is used for pointers...
  err = cuParamSetv(kernel, offset, &dest, sizeof(CUdeviceptr));
  gpu.checkForCudaError(err);
  offset += sizeof(CUdeviceptr);

  offset = gpu.align(offset, __alignof(int));
  err = cuParamSeti(kernel, offset, 4); // cuParamSeti is used for integers.
  gpu.checkForCudaError(err);
  offset += sizeof(int);
  offset = gpu.align(offset, __alignof(int));
  err = cuParamSeti(kernel, offset, 34);
  gpu.checkForCudaError(err);
  offset += sizeof(int);
  err = cuParamSetSize(kernel, offset);
  gpu.checkForCudaError(err);

  // Launch the grid
  cout << "Launch the grid" << endl;
  err = cuFuncSetBlockShape(kernel, 1, 1, 1);
  gpu.checkForCudaError(err);
  err = cuLaunchGrid(kernel, 1, 1);
  gpu.checkForCudaError(err);

  // Read the result off of the GPU
  cout << "Read the result" << endl;
  int result = 0;
  err = cuMemcpyDtoH(&result, dest, sizeof(int));
  gpu.checkForCudaError(err);

  cout << "Sum of 4 and 34 is " << result << endl;

  if (result == 38)
    ut.passes("Sum of 4 and 34 is 38.");
  else
    ut.failure("Sum of 4 and 34 was incorrect.");

  // deallocate memory, free the context.
  cout << "deallocate device memory." << endl;
  err = cuMemFree(dest);
  gpu.checkForCudaError(err);

  return;
}

//---------------------------------------------------------------------------//
// vector_add
//---------------------------------------------------------------------------//

void vector_add(rtt_dsxx::ScalarUnitTest &ut) {
  using namespace std;

  cout << "\nStarting gpu_hello_driver_api::vector_add()...\n" << endl;

  // Create a GPU_Device object.
  // Initialize the CUDA library and sets device and context handles.
  rtt_device::GPU_Device gpu;

  // Load the module, must compile the kernel with nvcc -ptx -m32 kernel.cu
  rtt_device::GPU_Module myModule("vector_add.cubin");

  // Host data
  size_t len(1024);
  size_t const threadsPerBlock(gpu.maxThreadsPerBlock());
  size_t const blocksPerGrid = (len + threadsPerBlock - 1) / threadsPerBlock;
  vector<double> aH(len);
  vector<double> bH(len);
  vector<double> cH(len, 0.0);
  vector<double> refH(len);
  genTestData(aH, bH, refH);

  // Load the kernel from the module
  CUfunction kernel;
  cudaError_enum err =
      cuModuleGetFunction(&kernel, myModule.handle(), "vector_add");
  gpu.checkForCudaError(err);

  // Allocate some memory for the result
  CUdeviceptr d_A, d_B, d_C;
  err = cuMemAlloc(&d_A, len * sizeof(double));
  gpu.checkForCudaError(err);
  err = cuMemAlloc(&d_B, len * sizeof(double));
  gpu.checkForCudaError(err);
  err = cuMemAlloc(&d_C, len * sizeof(double));
  gpu.checkForCudaError(err);

  // Copy host data to device
  err = cuMemcpyHtoD(d_A, &aH[0], len * sizeof(double));
  gpu.checkForCudaError(err);
  err = cuMemcpyHtoD(d_B, &bH[0], len * sizeof(double));
  gpu.checkForCudaError(err);

  // This is the function signature
  void *args[] = {&d_A, &d_B, &d_C, &len};

  // Execute the kernel
  err = cuLaunchKernel(kernel, blocksPerGrid, 1, 1, threadsPerBlock, 1, 1, 0, 0,
                       args, 0);
  gpu.checkForCudaError(err);

  // Copy result from device to host
  err = cuMemcpyDtoH((void *)(&cH[0]), d_C, len * sizeof(double));
  gpu.checkForCudaError(err);

  // Free device memory
  err = cuMemFree(d_A);
  gpu.checkForCudaError(err);
  err = cuMemFree(d_B);
  gpu.checkForCudaError(err);
  err = cuMemFree(d_C);
  gpu.checkForCudaError(err);

  // Check the result
  if (rtt_dsxx::soft_equiv(cH.begin(), cH.end(), refH.begin(), refH.end()))
    ut.passes("vector_add worked!");
  else
    ut.failure("vector_add failed.");

  return;
}

//---------------------------------------------------------------------------//
// vector_add_using_wrappers
//---------------------------------------------------------------------------//

void vector_add_using_wrappers(rtt_dsxx::ScalarUnitTest &ut) {
  using namespace std;

  cout << "\nStarting gpu_hello_driver_api::vector_add_using_wrappers()...\n"
       << endl;

  // Create a GPU_Device object.
  // Initialize the CUDA library and sets device and context handles.
  rtt_device::GPU_Device gpu;

  // Load the module, must compile the kernel with nvcc -ptx -m32 kernel.cu
  rtt_device::GPU_Module myModule("vector_add.cubin");

  // Host data
  size_t len(1024);
  size_t const threadsPerBlock(gpu.maxThreadsPerBlock());
  size_t const blocksPerGrid = (len + threadsPerBlock - 1) / threadsPerBlock;
  vector<double> aH(len);
  vector<double> bH(len);
  vector<double> cH(len, 0.0);
  vector<double> refH(len);
  genTestData(aH, bH, refH);

  // Load the kernel from the module
  CUfunction kernel = myModule.getModuleFunction("vector_add");

  // Allocate some memory for the result
  unsigned const nbytes = len * sizeof(double);
  CUdeviceptr d_A = gpu.MemAlloc(nbytes);
  CUdeviceptr d_B = gpu.MemAlloc(nbytes);
  CUdeviceptr d_C = gpu.MemAlloc(nbytes);

  // Copy host data to device
  gpu.MemcpyHtoD(d_A, &aH[0], nbytes);
  gpu.MemcpyHtoD(d_B, &bH[0], nbytes);

  // This is the function signature
  void *args[] = {&d_A, &d_B, &d_C, &len};

  // Execute the kernel
  cudaError_enum err = cuLaunchKernel(kernel, blocksPerGrid, 1, 1,
                                      threadsPerBlock, 1, 1, 0, 0, args, 0);
  gpu.checkForCudaError(err);

  // Copy result from device to host
  gpu.MemcpyDtoH((void *)(&cH[0]), d_C, nbytes);

  // Free device memory
  gpu.MemFree(d_A);
  gpu.MemFree(d_B);
  gpu.MemFree(d_C);

  // Check the result
  if (rtt_dsxx::soft_equiv(cH.begin(), cH.end(), refH.begin(), refH.end()))
    ut.passes("vector_add worked!");
  else
    ut.failure("vector_add failed.");

  return;
}

//---------------------------------------------------------------------------//
// Main
//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  using namespace std;

  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    query_device(ut);
    simple_add(ut);
    vector_add(ut);
    vector_add_using_wrappers(ut);
  } catch (exception &err) {
    cout << "ERROR: While testing gpu_hello_driver_api, " << err.what() << endl;
    ut.numFails++;
  } catch (...) {
    cout << "ERROR: While testing gpu_hello_driver_api, "
         << "An unknown exception was thrown." << endl;
    ut.numFails++;
  }
  return ut.numFails;
}

//---------------------------------------------------------------------------//
// end of gpu_hello_driver_api.cc
//---------------------------------------------------------------------------//
