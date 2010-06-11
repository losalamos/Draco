//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI_scatterv_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:44:54 2002
 * \brief  C4 MPI determinate and indeterminate scatterv instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>
#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "scatterv.t.hh"

namespace rtt_c4
{
using std::vector;

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF NON-BLOCKING SEND/RECEIVE
//---------------------------------------------------------------------------//

template
void indeterminate_scatterv(std::vector<std::vector<unsigned> > &outgoing_data,
                            std::vector<unsigned> &incoming_data);

template
void determinate_scatterv(std::vector<std::vector<unsigned> > &outgoing_data,
                          std::vector<unsigned> &incoming_data);

template
void determinate_scatterv(std::vector<std::vector<int> > &outgoing_data,
                          std::vector<int> &incoming_data);

template
void determinate_scatterv(std::vector<std::vector<double> > &outgoing_data,
                          std::vector<double> &incoming_data);

} // end namespace rtt_c4


//---------------------------------------------------------------------------//
//                              end of C4_MPI_scatterv_pt.cc
//---------------------------------------------------------------------------//
