//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/vector_add.cu
 * \author Kelly Thompson
 * \date   
 * \brief  Small kernel code for testing GPU Device framework.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

extern "C"
{
    __global__ void
    vector_add(double const *A, double const *B, double *C, int const N)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x; 
        // if(i%512==0)
        //     printf("index %d\n",i);
        if (i < N) 
            C[i] = A[i] + B[i];
    }
}

//---------------------------------------------------------------------------//
// end of vector_add.cu
//---------------------------------------------------------------------------//
