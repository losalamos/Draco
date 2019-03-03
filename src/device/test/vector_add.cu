//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/vector_add.cu
 * \author Kelly Thompson
 * \date
 * \brief  Small kernel code for testing GPU Device framework.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

__global__ void vector_add(double const *A_dev, double const *B_dev,
                           double *C_dev, int const N) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  /* if(i%512==0)
         *   printf("index %d\n",i); */
  if (i < N)
    C_dev[i] = A_dev[i] + B_dev[i];
}

//----------------------------------------------------------------------------//
// end of vector_add.cu
//----------------------------------------------------------------------------//
