//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/gatherv.t.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI template implementation.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_gatherv_t_hh
#define c4_gatherv_t_hh

#include "c4/config.h"
#include "gatherv.hh"
#include "C4_Functions.hh"
#include "ds++/Assert.hh"
#include <algorithm>

namespace rtt_c4
{
using std::vector;
using std::copy;

//---------------------------------------------------------------------------//
// EXCHANGE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template<class T>
void indeterminate_gatherv(vector<T> &outgoing_data,
                           vector<vector<T> >       &incoming_data)
{
#ifdef C4_MPI
    { // This block is a no-op for with-c4=scalar 

        unsigned const N = rtt_c4::nodes();

        incoming_data.resize(N);

        int count = outgoing_data.size();
        if (rtt_c4::node()==0)
        {
            vector<int> counts(N), displs(N);
            gather(&count, &counts[0], 1);
            unsigned total_count = 0;
            for (unsigned p=0; p<N; ++p)
            {
                displs[p] = total_count;
                total_count += counts[p];
            }
            
            vector<T> recbuf(total_count);
            rtt_c4::gatherv(&outgoing_data[0],
                            outgoing_data.size(),
                            &recbuf[0],
                            &counts[0],
                            &displs[0]);
            
            for (unsigned p=0; p<N; ++p)
            {
                incoming_data[p].resize(counts[p]);
                copy(recbuf.begin()+displs[p],
                     recbuf.begin()+displs[p]+counts[p],
                     incoming_data[p].begin());
            }

        }
        else
        {
            gather(&count, static_cast<int*>(NULL), 1);
            gatherv(&outgoing_data[0],
                    outgoing_data.size(),
                    static_cast<T*>(NULL),
                    0,
                    static_cast<int*>(NULL));
        }
    }
#else
    {
        // Only need to copy outgoing to incoming
        incoming_data.resize(0);
        incoming_data.resize(1, outgoing_data);
    }
#endif // C4_MPI

    return;
}

} // end namespace rtt_c4

#endif // c4_gatherv_t_hh

//---------------------------------------------------------------------------//
//                              end of c4/gatherv.t.hh
//---------------------------------------------------------------------------//
