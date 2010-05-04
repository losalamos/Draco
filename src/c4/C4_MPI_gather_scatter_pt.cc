//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI_gather_scatter_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:44:54 2002
 * \brief  C4 MPI non-blocking send/recv instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>

#ifdef C4_MPI

#include "C4_MPI.t.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF GATHER/SCATTER
//---------------------------------------------------------------------------//

template int gather(int *send_buffer,
                    int *receive_buffer,
                    int size);

template int gather(unsigned *send_buffer,
                    unsigned *receive_buffer,
                    int size);

template int scatter(unsigned *send_buffer,
                     unsigned *receive_buffer,
                     int size);
template
int gatherv(unsigned *send_buffer,
            int send_size,
            unsigned *receive_buffer,
            int *receive_sizes,
            int *receive_displs);

} // end namespace rtt_c4

#endif // C4_MPI

//---------------------------------------------------------------------------//
//                              end of C4_MPI_gather_scatter_pt.cc
//---------------------------------------------------------------------------//
