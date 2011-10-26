//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/gpu_hello_rt_api.hh
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 25 15:28:48 2011
 * \brief  
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

// All this garbage suppresses warnings found in "cuda.h".
// http://wiki.services.openoffice.org/wiki/Writing_warning-free_code#When_all_else_fails
#if defined __GNUC__
#pragma GCC system_header
#elif defined __SUNPRO_CC
#pragma disable_warn
#elif defined _MSC_VER
#pragma warning(push, 1)
#endif
#include <cuda_runtime_api.h>
#if defined __SUNPRO_CC
#pragma enable_warn
#elif defined _MSC_VER
#pragma warning(pop)
#endif

// #define Error(format, args...) (error_and_exit)("%s:%d: " format, __FILE__, __LINE__, ##args)


