/*-----------------------------------*-C-*-----------------------------------*/
/*!
 * \file   c4_mpi.h 
 * \author Thomas M. Evans 
 * \date   Fri Jan  8 15:06:30 1999
 * \brief  put the right includes for MPI header files
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
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

#endif /* __c4_c4_mpi_h__ */

/*---------------------------------------------------------------------------*/
/* end of c4_mpi.h */
/*---------------------------------------------------------------------------*/
