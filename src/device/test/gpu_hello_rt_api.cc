//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/gpu_hello_rt_api.cc
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 25 15:28:48 2011
 * \brief  Simple test of the CUDA Runtime API.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "gpu_hello_rt_api.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <iostream>
#include <string>
#include <vector>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

void copy_and_retrieve(rtt_dsxx::ScalarUnitTest &ut) {
  using namespace std;

  // Size of array
  size_t const count(100);

  // Host arrays
  vector<float> inH(count, -1.0);
  vector<float> outH(count, -1.0);

  // Allocate device array
  float *inD = NULL; // inD = NULL;
  cudaError_t ret = cudaMalloc((void **)&inD, count * sizeof(float));
  if (ret != cudaSuccess)
    ut.failure(string("cudaMalloc returned \"") +
               string(cudaGetErrorString(ret)) + string("\""));
  else
    ut.passes("Allocated inD array on device.");

  // Initialize host data
  for (size_t i = 0; i < count; ++i)
    inH[i] = i;

  // Copy to device
  ret = cudaMemcpy(inD, &inH[0], count * sizeof(float), cudaMemcpyHostToDevice);
  if (ret != cudaSuccess)
    ut.failure(string("cudaMemcpy returned \"") +
               string(cudaGetErrorString(ret)) + string("\""));
  else
    ut.passes("Copied array from host to device.");

  // Copy back to host
  ret =
      cudaMemcpy(&outH[0], inD, count * sizeof(float), cudaMemcpyDeviceToHost);
  if (ret != cudaSuccess)
    ut.failure(string("cudaMemcpy returned \"") +
               string(cudaGetErrorString(ret)) + string("\""));
  else
    ut.passes("Copied array from device to host.");

  // Compare original and result arrays.
  size_t diff(0);
  for (size_t i = 0; i < count; ++i) {
    if (inH[i] != outH[i]) {
      diff = 1;
      cout << "diff: inH[" << i << "] = " << inH[i] << "; outH[" << i
           << "] = " << outH[i] << endl;
    }
  }
  if (diff != 0)
    ut.failure("inH != outH");
  else
    ut.passes("inH == outH");

  // Clean up
  ret = cudaFree(inD);
  if (ret != cudaSuccess)
    ut.failure(string("cudaFree returned \"") +
               string(cudaGetErrorString(ret)) + string("\""));
  else
    ut.passes("Freed allocated memory on device.");

  return;
}

//---------------------------------------------------------------------------//
// Main
//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  using namespace std;

  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    copy_and_retrieve(ut);
  } catch (exception &err) {
    cout << "ERROR: While testing gpu_hello_rt_api, " << err.what() << endl;
    ut.numFails++;
  } catch (...) {
    cout << "ERROR: While testing gpu_hello_rt_api, "
         << "An unknown exception was thrown." << endl;
    ut.numFails++;
  }
  return ut.numFails;
}

//---------------------------------------------------------------------------//
// end of gpu_hello_rt_api.cc
//---------------------------------------------------------------------------//
