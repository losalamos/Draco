/*-----------------------------------*-C-*-----------------------------------*/
/* c4_mpi.h */
/* Thomas M. Evans */
/* Fri Jan  8 15:06:30 1999 */
/*---------------------------------------------------------------------------*/
/* @> put the right includes for MPI header files                            */
/*---------------------------------------------------------------------------*/

#ifndef rtt_c4_c4_mpi_h
#define rtt_c4_c4_mpi_h

#include "c4/config.h"

/* defined in ac_vendors.m4, location of <mpi.h> */
#include <mpi.h>

#ifdef MPI_MAX_PROCESSOR_NAME
#define DRACO_MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#else

#endif

#endif                          /* __c4_c4_mpi_h__ */

/*---------------------------------------------------------------------------*/
/*                              end of c4_mpi.h */
/*---------------------------------------------------------------------------*/
