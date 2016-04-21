//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/accumulatev.t.hh
 * \author Kelly Thompson, Thomas M. Evans, Bob Webster
 * \date   Monday, Nov 05, 2012, 13:41 pm
 * \brief  C4 accumulatev template implementation.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: accumulatev.t.hh 6288 2011-12-04 03:43:52Z kellyt $
//---------------------------------------------------------------------------//

#ifndef c4_accumulatev_t_hh
#define c4_accumulatev_t_hh

#include "C4_Functions.hh"
#include "gatherv.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// ACCUMULATE
//---------------------------------------------------------------------------//
#ifdef C4_MPI

template<typename T, typename Tciter, typename BinaryOp>
void accumulatev(Tciter   localBegin,
                 Tciter   localEnd,
                 T        init,
                 BinaryOp op)
{
    // one processor - nothing to do.
    if( rtt_c4::nodes() == 1 ) return;

    size_t const ndata( std::distance(localBegin,localEnd) );

    std::vector<T> data(ndata);
    std::copy(localBegin,localEnd,data.begin());

    std::vector< std::vector< T > > incoming_data(rtt_c4::nodes());

    // Bring everything to PE0
    indeterminate_gatherv( data, incoming_data );

    // accumulate
    if( rtt_c4::node() == 0 )
    {
        for( size_t i=0; i<ndata; ++i )
        {
            data[i]=init;
            for( int rank=0; rank<rtt_c4::nodes(); ++rank )
                data[i]=op(data[i],incoming_data[rank][i]);
        }
    }

    // Send the data to all procs
    broadcast( data.begin(), data.end(), localBegin );
    // Also update localBegin for processor 0
    if( rtt_c4::node()==0 )
        std::copy( data.begin(), data.end(), localBegin );

    return;
}

#else // NOT C4_MPI

template<typename T, typename Tciter, typename BinaryOp>
void accumulatev(Tciter   ,
                 Tciter   ,
                 T        ,
                 BinaryOp )
{
    // Nothing to do.
    return;
}

#endif // C4_MPI

} // end namespace rtt_c4

#endif // c4_accumulatev_t_hh

//---------------------------------------------------------------------------//
// end of c4/accumulatev.t.hh
//---------------------------------------------------------------------------//
