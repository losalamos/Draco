//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   lapack_wrap/Blas_Prototypes.hh
 * \author Thomas M. Evans
 * \date   Thu Aug 29 11:23:27 2002
 * \brief  Header declaring BLAS prototypes
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __lapack_wrap_Blas_Prototypes_hh__
#define __lapack_wrap_Blas_Prototypes_hh__

#include <lapack_wrap/config.h>

extern "C"
{
    // >>> LEVEL 1 BLAS

    // y <- x
    void SCOPY(int *, float *, int *, float *, int *);
    void DCOPY(int *, double *, int *, double *, int *);

    // x <- ax
    void SSCAL(int *, float *, float *, int *);
    void DSCAL(int *, double *, double *, int *);

    // dot <- x^T y
    float  SDOT(int *, float *, int *, float *, int *);
    double DDOT(int *, double *, int *, double *, int *);

    // y <- ax + y
    void SAXPY(int *, float *, float *, int *, float *, int *);
    void DAXPY(int *, double *, double *, int *, double *, int *);

    // nrm2 <- ||x||_2
    float SNRM2(int *, float *, int *);
    double DNRM2(int *, double *, int *);

} // end of extern "C"

#endif                          // __lapack_wrap_Blas_Prototypes_hh__

//---------------------------------------------------------------------------//
//                              end of lapack_wrap/Blas_Prototypes.hh
//---------------------------------------------------------------------------//
