/*-----------------------------------*-C-*-----------------------------------*/
/* file  : c4_omp.h
 * author: Kelly Thompson
 * date  : Wednesday, Nov 23, 2011, 08:42 am
 * brief : wrapper for omp.h
 * note  : Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
/*---------------------------------------------------------------------------*/

#ifndef rtt_c4_c4_omp_h
#define rtt_c4_c4_omp_h

#include "c4/config.h"

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

//---------------------------------------------------------------------------*/
// Wrap OpenMP API calls
//---------------------------------------------------------------------------*/
namespace rtt_c4 {

#ifdef OPENMP_FOUND

/*! Set the number of OpenMP threads
 *\param n: desired number of threads
 */
void set_omp_num_threads(int n) { omp_set_num_threads(n); }

/*! Get the number of OpenMP threads in current section
 *\return: current number of threads
 */
int get_omp_num_threads() { return omp_get_num_threads(); }

/*! Get the maximum number of OpenMP threads
 *\return: maximum number of threads
 */
int get_omp_max_threads() { return omp_get_max_threads(); }

#else

/*! Bypass setting the number of OpenMP threads */
void set_omp_num_threads(int /*n*/) { return; }

/*! Bypass getting the current number of OpenMP threads
 *\return: 1 (always)
 */
int get_omp_num_threads() { return 1; }

/*! Bypass getting the maximum number of OpenMP threads
 *\return: 1 (always)
 */
int get_omp_max_threads() { return 1; }

#endif

} // namespace rtt_c4

#endif /* rtt_c4_c4_omp_h */

/*---------------------------------------------------------------------------*/
/* end of c4_omp.h */
/*---------------------------------------------------------------------------*/
