//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Processor_Group.i.hh
 * \author Kent Budge
 * \date   Fri Oct 20 13:49:10 2006
 * \brief  Template method definitions of class Processor_Group
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_Processor_Group_i_hh
#define c4_Processor_Group_i_hh

#include "ds++/Assert.hh"
#include "MPI_Traits.hh"
#ifdef C4_MPI
namespace rtt_c4
{
//---------------------------------------------------------------------------//
template<class T>
void Processor_Group::sum(std::vector<T> &x)
{
    // copy data into send buffer
    std::vector<T> y = x;
    
    // do global MPI reduction (result is on all processors) into x
    int status = MPI_Allreduce(&y[0],
                               &x[0],
                               y.size(),
                               rtt_c4::MPI_Traits<T>::element_type(),
                               MPI_SUM,
                               comm_);
    Insist(status==0, "MPI_Allreduce failed");
}

} // end namespace rtt_c4

#endif  // C4_MPI
#endif // c4_Processor_Group_i_hh

//---------------------------------------------------------------------------//
//              end of c4/Processor_Group.i.hh
//---------------------------------------------------------------------------//
