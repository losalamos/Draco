//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   lapack_wrap/Release.hh
 * \author Thomas M. Evans
 * \date   Thu Aug 29 11:06:46 2002
 * \brief  Release function for the lapack_wrap library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __lapack_wrap_Release_hh__
#define __lapack_wrap_Release_hh__

//===========================================================================//
/*!
 * \page lapack_wrap_overview Overview of the lapack_wrap package
 * \version 1_0_0
 * \author  Tom Evans
 * 
 * The lapack_wrap package contains a C++ functional interface to BLAS and
 * LAPACK functions. 
 *
 * BLAS functions are included in the headers Blas_Level_1.hh (level 1 BLAS),
 * Blas_Level_2.hh (level 2 BLAS), and Blas_Level_3.hh (level 3 BLAS).  To
 * include the C++ wrapped BLAS include the header <Blas.hh>.  This header
 * includes all three Blas_Level headers.
 *
 * LAPACK functionality has not been wrapped yet.
 * 
 */
//===========================================================================//
/*!
 * \namespace rtt_lapack_wrap
 *
 * \brief Namespace that contains the lapack_wrap package classes and
 * variables.
 *
 */
//===========================================================================//

#include <string>

namespace rtt_lapack_wrap 
{
    const std::string release();
}

#endif                          // __lapack_wrap_Release_hh__

//---------------------------------------------------------------------------//
//                           end of lapack_wrap/Release.hh
//---------------------------------------------------------------------------//
