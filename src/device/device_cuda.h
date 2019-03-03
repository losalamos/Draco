//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/device_cuda.h
 * \author Kelly (KT) Thompson
 * \brief  Wrap the cuda.h header while preventing comiler warnings about
 *         vendor code.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#ifndef device_device_cuda_h
#define device_device_cuda_h

// All this garbage suppresses warnings found in "cuda.h".
// http://wiki.services.openoffice.org/wiki/Writing_warning-free_code#When_all_else_fails
#if defined __GNUC__
#pragma GCC system_header
// Intel defines __GNUC__ by default
#ifdef __INTEL_COMPILER
#pragma warning push
#endif
#elif defined __SUNPRO_CC
#pragma disable_warn
#elif defined _MSC_VER
#pragma warning(push, 1)
#endif

#include <cuda.h>

#if defined __GNUC__
#pragma GCC system_header
#ifdef __INTEL_COMPILER
#pragma warning pop
#endif
#elif defined __SUNPRO_CC
#pragma enable_warn
#elif defined _MSC_VER
#pragma warning(pop)
#endif

#endif // device_device_cuda_h

//---------------------------------------------------------------------------//
// end of device_cuda.h
//---------------------------------------------------------------------------//
