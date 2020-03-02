//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   VendorChecks/test/superlu-dist-wrapper.h
 * \date   Wednesday, May 01, 2019, 09:52 am
 * \brief  Wrapper for superlu_ddefs.h to allow warning suppression
 * \note   Copyright (C) 2019, Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC system_header
#endif
#include <superlu_ddefs.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

//----------------------------------------------------------------------------//
// End superlu-dist-wrapper.h
//----------------------------------------------------------------------------//
