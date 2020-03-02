//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/C4_MPI_gather_scatter_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:44:54 2002
 * \brief  C4 MPI non-blocking send/recv instantiations.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "c4/config.h"

#ifdef C4_MPI
#include "C4_MPI.t.hh"
#else
#include "C4_Serial.t.hh"
#endif

namespace rtt_c4 {

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF GATHER/SCATTER
//----------------------------------------------------------------------------//

template DLL_PUBLIC_c4 int gather<unsigned>(unsigned *send_buffer,
                                            unsigned *receive_buffer, int size);

template DLL_PUBLIC_c4 int gather<int>(int *send_buffer, int *receive_buffer,
                                       int size);

template DLL_PUBLIC_c4 int gather<char>(char *send_buffer, char *receive_buffer,
                                        int size);

//----------------------------------------------------------------------------//
template DLL_PUBLIC_c4 int allgather<int>(int *send_buffer, int *receive_buffer,
                                          int size);

//----------------------------------------------------------------------------//
template DLL_PUBLIC_c4 int gatherv<unsigned>(unsigned *send_buffer,
                                             int send_size,
                                             unsigned *receive_buffer,
                                             int *receive_sizes,
                                             int *receive_displs);

template DLL_PUBLIC_c4 int gatherv<int>(int *send_buffer, int send_size,
                                        int *receive_buffer, int *receive_sizes,
                                        int *receive_displs);

template DLL_PUBLIC_c4 int gatherv<double>(double *send_buffer, int send_size,
                                           double *receive_buffer,
                                           int *receive_sizes,
                                           int *receive_displs);
template DLL_PUBLIC_c4 int gatherv<char>(char *send_buffer, int send_size,
                                         char *receive_buffer,
                                         int *receive_sizes,
                                         int *receive_displs);

//----------------------------------------------------------------------------//
template DLL_PUBLIC_c4 int
scatter<unsigned>(unsigned *send_buffer, unsigned *receive_buffer, int size);

template DLL_PUBLIC_c4 int scatter<int>(int *send_buffer, int *receive_buffer,
                                        int size);

template DLL_PUBLIC_c4 int scatterv<unsigned>(unsigned *send_buffer,
                                              int *send_sizes, int *send_displs,
                                              unsigned *receive_buffer,
                                              int receive_size);

template DLL_PUBLIC_c4 int scatterv<int>(int *send_buffer, int *send_sizes,
                                         int *send_displs, int *receive_buffer,
                                         int receive_size);

template DLL_PUBLIC_c4 int scatterv<double>(double *send_buffer,
                                            int *send_sizes, int *send_displs,
                                            double *receive_buffer,
                                            int receive_size);

} // end namespace rtt_c4

//----------------------------------------------------------------------------//
// end of C4_MPI_gather_scatter_pt.cc
//----------------------------------------------------------------------------//
