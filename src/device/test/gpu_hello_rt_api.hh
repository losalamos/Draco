//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/gpu_hello_rt_api.hh
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 25 15:28:48 2011
 * \brief  Wrap the cuda_runtime_api.h header while preventing comiler
 *         warnings about vendor code.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

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

#include <cuda_runtime_api.h>

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

// #define Error(format, args...) (error_and_exit)("%s:%d: " format, __FILE__, __LINE__, ##args)

//---------------------------------------------------------------------------//
// end of gpu_hello_rt_api.hh
//---------------------------------------------------------------------------//
