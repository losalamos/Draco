//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   config/query_gpu.cu
 * \author Alex Long
 * \brief  Small CUDA code that prints the architecture version, used by CMake
 * \date   Thu Mat 21 15:53:51 2019
 * \note   Copyright (C) 2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

// NOTE: This code is from
// wagonhelm.github.io/articles/2018-03/detecting-cuda-capability-with-cmake

#include <stdio.h>

int main(int argc, char **argv) {
  cudaDeviceProp dP;
  float min_cc = 3.0;

  int rc = cudaGetDeviceProperties(&dP, 0);
  if (rc != cudaSuccess) {
    cudaError_t error = cudaGetLastError();
    printf("CUDA error: %s", cudaGetErrorString(error));
    return rc; /* Failure */
  }
  if ((dP.major + (dP.minor / 10)) < min_cc) {
    printf("Min Compute Capability of %2.1f required:  %d.%d found\n Not "
           "Building CUDA Code",
           min_cc, dP.major, dP.minor);
    return 1; /* Failure */
  } else {
    printf("-arch=sm_%d%d", dP.major, dP.minor);
    return 0; /* Success */
  }
}
//---------------------------------------------------------------------------//
//  end of query_gpu.cu
//---------------------------------------------------------------------------//
