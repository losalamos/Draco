//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/swap.t.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI template implementation.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_swap_t_hh
#define c4_swap_t_hh

#include "c4/config.h"
#include "ds++/Assert.hh"
#include "swap.hh"
#include "C4_Functions.hh"

namespace rtt_c4
{
using std::vector;

//---------------------------------------------------------------------------//
// EXCHANGE
//---------------------------------------------------------------------------//

template<class T>
void determinate_swap(vector<unsigned>   const &outgoing_pid,
                      vector<vector<T> > const &outgoing_data,
                      vector<unsigned>   const &incoming_pid,
                      vector<vector<T> >       &incoming_data,
                      int tag = C4_Traits<T*>::tag )
{
    Require(outgoing_pid.size()==outgoing_data.size());
    Require(incoming_pid.size()==incoming_data.size());

    unsigned incoming_processor_count = incoming_pid.size();
    unsigned outgoing_processor_count = outgoing_pid.size();

#ifdef C4_MPI
    { // This block is a no-op for with-c4=scalar 
        
        // Post the asynchronous sends.
        vector<C4_Req> outgoing_C4_Req(outgoing_processor_count);
        for (unsigned p=0; p<outgoing_processor_count; ++p)
        {
            outgoing_C4_Req[p] =
                rtt_c4::send_async(&outgoing_data[p][0],
                                   outgoing_data[p].size(),
                                   outgoing_pid[p],
                                   tag);
        }
        
        // Post the asynchronous receives
        vector<C4_Req> incoming_C4_Req(incoming_processor_count);
        for (unsigned p=0; p<incoming_processor_count; ++p)
        {
            incoming_C4_Req[p] =
                receive_async(&incoming_data[p][0],
                              incoming_data[p].size(),
                              incoming_pid[p],
                              tag);
        }
        
        // Wait for all the receives to complete.
        
        wait_all(incoming_processor_count, &incoming_C4_Req[0]);
        
        // Wait until all the posted sends have completed.
        
        wait_all(outgoing_processor_count, &outgoing_C4_Req[0]);
        
    }
#endif // C4_MPI

    return;
}

} // end namespace rtt_c4

#endif // c4_swap_t_hh

//---------------------------------------------------------------------------//
//                              end of c4/swap.t.hh
//---------------------------------------------------------------------------//
