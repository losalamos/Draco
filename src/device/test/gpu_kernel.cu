//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/gpu_kernel.cu
 * \author Kelly Thompson
 * \date
 * \brief  Small kernel code for testing GPU Device framework.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

extern "C" {
__global__ void sum(int *dest, int a, int b) {
  // Assuming a single thread, 1x1x1 block, 1x1 grid
  *dest = a + b;
}
}

//---------------------------------------------------------------------------//
// end of gpu_kernel.cu
//---------------------------------------------------------------------------//
