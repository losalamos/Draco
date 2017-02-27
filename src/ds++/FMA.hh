//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    ds++/FMA.hh
 * \author  Kelly Thompson
 * \date    Thursday, Feb 09, 2017, 11:22 am
 * \brief   Provide extra control for FMA operations.
 * \note    Copyright (C) 2017 Los Alamos National Security, LLC.
 *          All rights reserved.
 *
 * Intel Haswell and later (and also modern AMD cpus), have hardware FMA
 * features. Using \c fma(c,b,a) is more accurate than \c a*b+c because there is
 * less roundoff error. This difference in accuracy causes issues for IMC
 * because solutions are not consistent between older and newer hardware.
 *
 * To improve solution consistency across platforms, we can choose to use a
 * software-based infinite precision fma call.
 *
 * \note Intel's code for detecting FMA availability on hardware.
 * https://software.intel.com/en-us/node/405250?language=es&wapkw=avx2+cpuid
 */
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_FMA_hh
#define rtt_dsxx_FMA_hh

#include "ds++/config.h"
#include <cmath>

//---------------------------------------------------------------------------//
// HAVE_HARDWARE_FMA is set by config/platform_checks.cmake.

#ifdef HAVE_HARDWARE_FMA

/* This will be the default for new machines. Both versions should be equally
 * accurate, choose whatever is faster. FP_FAST_FMA is set by the compiler
 *
 * \ref http://en.cppreference.com/w/c/numeric/math/fma
 */
#ifdef FP_FAST_FMA
#define FMA(a, b, c) fma(a, b, c)
#else
#define FMA(a, b, c) ((a) * (b) + c)
#endif

#else /* HAVE_HARDWARE_FMA is false */

/*! For older hardware that does not support FMA natively, provide a macro that
 *  will call an infinite precision fma function so that numerical
 *  reproducibility is enhanced (if FP_ACCURATE_FMA is defined). This allows
 *  solutions generated on older machines to match solutions generated on newer
 *  machines. */
#ifndef FP_ACCURATE_FMA
#error "Must set FP_ACCURATE_FMA to 0 or 1."
#endif

#if FP_ACCURATE_FMA > 0
#define FMA(a, b, c) fma(a, b, c)
#else
#define FMA(a, b, c) ((a) * (b) + c)
#endif

#endif /* HAVE_HARDWARE_FMA */

#endif // rtt_dsxx_FMA_hh

//---------------------------------------------------------------------------//
// end of ds++/FMA.hh
//---------------------------------------------------------------------------//
