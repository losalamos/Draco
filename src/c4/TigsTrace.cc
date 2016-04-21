//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/TigsTrace.cc
 * \author Thomas M. Evans
 * \date   Fri Sep 30 12:53:25 2005
 * \brief  TigsTrace function definitions.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Tigs.cc 6056 2012-06-19 19:05:27Z kellyt $
//---------------------------------------------------------------------------//

#include "TigsTrace.hh"
#include "accumulatev.hh"
#include <cstdlib>
#include <functional>

namespace rtt_c4
{

//===========================================================================//
// COMM::TRACE CLASS MEMBER DEFINITIONS
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTORS AND ASSIGNMENT
//---------------------------------------------------------------------------//
/*!
 * \brief TigsTrace constructor.
 *
 * A "trace" is the object that handles hidden communication associated
 * with an indirection map into some global range of indices.
 *
 * The trace is defined by a map JMap:I->J given to the
 * trace constructor as a container of ints that with a length equal
 * to the local length of the domain (I) space. Each of those ints is
 * the Global index of the range (j) location to which the local I
 * location is mapped.
 *
 * Data stored on the range of the map is assumed to be layed out
 * in processor rank ascending order. If J# is the number of range
 * values stored on processor #, then the global index is an index in
 * the collation of the local processor storage of the range, that is,
 * it indexes into the sequence as
 *
 * \verbatim
 * D[0,0] ... D[J0-1,0],D[0,1] ... D[J1-0,1], ... ,D[0,N] ... D[JN-1,N]
 * \endverbatim
 *
 * where D[0,n] is the first element of the container D stored on
 * processor n, and Jn is the length of that conatainer.
 *
 *    \param JMap The values of the above indirection for the local domain
 *    \param J    the value of Jn for the current process.
 */
TigsTrace::TigsTrace( std::vector<int> const & M, unsigned const J )
    : onProcDomain( M.size() ),
      onProcRange( J ),
      IM(),
      IMV(),
      counts(),
      IsideBufferSize(0),
      JsideBufferSize(0),
      IsideConnects(),
      IsideIndirect(), 
      JsideConnects(),
      JsideIndirect(),
      BmapList(),
      countsList() 
{
    if( rtt_c4::nodes()==1 )
    {
        // form inverse map
        for( size_t i=0; i<onProcDomain; i++ )
        {
            int j=M[i];
            TigsComm_map::iterator itm=IM.find(j);
            if( itm != IM.end() )
            {
                // Append to entry
                (itm->second).push_back(i);
            }
            else
            {
                // Add a new entry...
                std::vector <int> ilist(1,i);
                IM.insert(std::make_pair(j,   ilist));
            }
        }
        //
        // Express the inverse in a CSR like fashion
        size_t total_count=0;
        counts.resize(onProcRange,0);
        for( size_t j=0; j<onProcRange; j++ )
        {
            TigsComm_map::iterator itm=IM.find(j);
            if( itm != IM.end() )
            {
                std::vector<int> &V(itm->second);
                counts[j]=V.size();
                total_count+=V.size();
            }
        }
        size_t k=0;
        IMV.resize(total_count);
        for( size_t j=0; j<onProcRange; j++ )
        {
            TigsComm_map::iterator itm=IM.find(j);
            if( itm != IM.end() )
            {
                std::vector<int> &V(itm->second);
                for( size_t l=0; l<V.size(); l++ )
                {
                    IMV[k++]= V[l];
                }
            }
        }
        
        JsideBufferSize=total_count;
    }
    else
    {
#ifdef C4_MPI
        // first get all processors the information needed to map global
        // indices to (processor, localindex) pairs. Jstar[pe] is going to
        // contain the CSN of the first element of the range array on process
        // pe.
        std::vector<int> Jstarb(rtt_c4::nodes(),0);
        std::vector<int> Jstar(rtt_c4::nodes()+1,0);
        Jstarb[rtt_c4::node()]=J;
        int const initial_value(0);
        rtt_c4::accumulatev( Jstarb.begin(), Jstarb.end(),
                             initial_value,  std::plus<int>());
        for( int pe=1; pe<rtt_c4::nodes()+1; pe++ )
            Jstar[pe]=Jstar[pe-1] + Jstarb[pe-1];

        // Now form the maps from local index "i" of the map M processor and
        // local index on that processor.
        std::vector<int> Pbar(M.size(),-1);
        std::vector<int> jbar(M.size(),-1);
        for( size_t i=0; i<M.size(); i++ )
        {
            for( int pe=0; pe<rtt_c4::nodes(); pe++ )
            {
                if( M[i] >= Jstar[pe] )
                {
                    Pbar[i]=pe;
                    jbar[i]=M[i]-Jstar[pe];
                }
                else
                {
                    pe=rtt_c4::nodes()+1;
                }
            }
        }

        // Form maps of data movement to target processor Itilde and Jtilde
        // are maps of std::vectors. For the same processor number, these
        // std::vectors are the same length and store the local index of the
        // src data and the local index (ont the target processor) of the
        // destination data.
        TigsComm_map Itilde;
        TigsComm_map Jtilde;
        for( size_t i=0; i<M.size(); i++ )
        {    
            // make sure the processor is actually communicating with another
            // pe
            if( Pbar[i] >= 0 )
            {
                Itilde[Pbar[i]].push_back(i); 
                Jtilde[Pbar[i]].push_back(jbar[i]);
            }
        }

        // recv flag to eliminate any send-to-selfs
        size_t recv_from_self = 0;

        // now check to see if we have a send-to-self in Itilde
        if( Itilde.count(rtt_c4::node()) == 1 ) recv_from_self = 1;
        
        // All processors now need to know how many messages they should
        // expect during a gather_list or scatter_list operation. gathers and
        // scatters are implemented by post processing those more basic data
        // movements.
        std::vector<int> numSrcProcs(rtt_c4::nodes(),0);
        TigsComm_map::iterator itb;
        TigsComm_map::iterator jtb;
        size_t numSend=0;
        for( itb=Itilde.begin(); itb!=Itilde.end() ; itb++ )
        {
            numSrcProcs[itb->first]++;
            numSend++;
        }
        // save the source side indirection maps as std::vectors (maps are
        // overkill only really useful in creating the connections)
        IsideConnects.resize(numSend);
        IsideIndirect.resize(numSend);
        std::vector< std::vector<int> > JsideIndirectA(numSend);
        itb=Itilde.begin();
        for(  size_t s=0; s<numSend ; itb++, s++)
        {
            jtb=Jtilde.find(itb->first);
            IsideConnects[s]  = itb->first;
            IsideIndirect[s]  = itb->second;
            JsideIndirectA[s] = jtb->second;
        }

        rtt_c4::accumulatev(numSrcProcs.begin(), numSrcProcs.end(),
                            initial_value, std::plus<int>() );
        rtt_c4::global_barrier();
        
        // Now we have enough data to make the maps, we need to get the Jside
        // information over to the processors that have the Jside data.
        // Currently the data is stored in JsideIndirectA, but that will be
        // moved to the processor local to the Jside data and stored as
        // JsideIndirect there.
        
        // Post the receives so that the processes can trade buffer size
        // information, buf size will receive 2 ints, the first will be the
        // source process rank, the second will be the amount of data coming
        // from the source process.        
        size_t const numRecv=numSrcProcs[rtt_c4::node()];
        std::vector< std::vector<int> > rbufsiz(numRecv,std::vector<int>(2));
        std::vector< std::vector<int> > sbufsiz(numSend,std::vector<int>(2));

        // add one to sizes so that dereferencing the std::vectors does not result
        // in an out-of-bounds read (ie. &rreqs[0] will have a valid pointer
        // to send)
        std::vector< rtt_c4::C4_Req > reqs_recv;
        for( size_t s=0; s<numSend; s++ )
        {
            sbufsiz[s][0]=rtt_c4::node();
            sbufsiz[s][1]=JsideIndirectA[s].size();
        }

        // the number of posted receives = numRecv - recv_from_self such that
        // numPostRecv == numSend || numPostRecv == numSend - 1
        Check( numRecv >= recv_from_self );
        size_t const numPostRecv = numRecv - recv_from_self;
        for( size_t r=0; r<numPostRecv ; r++ )
        {
            reqs_recv.push_back( rtt_c4::receive_async(
                                 &(rbufsiz[r][0]), 2, MPI_ANY_SOURCE ) );
        }
        for( size_t s=0; s<numSend ; s++ )
        {
            // only send if not send-to-self
            if( IsideConnects[s] != rtt_c4::node() )
            {
                rtt_c4::C4_Req c4req_send = rtt_c4::send_async(
                    &(sbufsiz[s][0]), 2, IsideConnects[s] );
                c4req_send.wait();
            }
            // otherwise, dump into the receive buffer directly
            else
            {
                Check( numPostRecv + 1 == numRecv );
                Check( Itilde.count(rtt_c4::node()) == 1 );

                // add sbufsiz[s] to the last rbufsize
                rbufsiz[numPostRecv] = sbufsiz[s];
            }
        }

        // wait on all communication
        for( size_t s=0; s<reqs_recv.size(); ++s )
            reqs_recv[s].wait();
        
        // wait for communications to finish
        rtt_c4::global_barrier();

        // dummy line to make sure compiler doesn't release sbufsiz early
        // if( numSend<0)
        //  Check((sbufsiz[0][0] > -10));

        // Now we have the data needed to allocate the tables for indirections
        // on both sides of the communication. Allocate the space and use a
        // second communication to move over the indirection lists...  Let's
        // be smart and order the receive buffers in increasing processor
        // rank. This makes testing easier.

        std::map<int, int> sortRcv;
        JsideConnects.resize(numRecv);
        JsideIndirect.resize(numRecv);
        for( size_t r=0; r<numRecv ; r++ )
        {
            sortRcv[rbufsiz[r][0]]=rbufsiz[r][1];
        }
        std::map<int, int>::iterator its=sortRcv.begin();
      
        for( int r=0; its!=sortRcv.end(); its++, r++ )
        {
            JsideConnects[r] =  its->first;;
            JsideIndirect[r].resize(its->second);
        }

        // store the recv-from-self index
        int recv_from_self_index = -1;
        Remember( size_t np = 0; );
        for( size_t r=0; r<numRecv ; r++ )
        {
            // only post receives if not recv-from-self
            if( JsideConnects[r] != rtt_c4::node() )
            {
                reqs_recv.push_back(
                    rtt_c4::receive_async( &(JsideIndirect[r][0]),
                                           JsideIndirect[r].size(),
                                           JsideConnects[r] ) );
                Remember( np++; );
            }
            // store the recv_from_self_index for future reference on the send
            else
            {
                Check (JsideConnects[r] == rtt_c4::node());
                recv_from_self_index = r;
            }
        }
        Check( np == numPostRecv );
        
        itb=Jtilde.begin();
        for( size_t s=0; s<numSend ; s++, itb++ )
        {
            // only send if not send-to-self
            if( IsideConnects[s] != rtt_c4::node() )
            {
                rtt_c4::C4_Req c4req_send = rtt_c4::send_async(
                    &(JsideIndirectA[s][0]), JsideIndirectA[s].size(),
                    IsideConnects[s] );
                c4req_send.wait();
            }
            // for send-to-self write directly into receive buffer
            else
            {
                Check( recv_from_self_index >= 0 );
                Check( recv_from_self_index < static_cast<int>(numRecv) );

                // put send buffer into receive buffer directly
                JsideIndirect[recv_from_self_index] = JsideIndirectA[s];
            }
        }

        // wait on the receives
        for( size_t s=0; s<reqs_recv.size(); ++s )
            // if( JsideConnects[s] != rtt_c4::node() )
                reqs_recv[s].wait();

        // wait for communications to finish
        rtt_c4::global_barrier();

        // finally, before leaving lets set the total buffer size needed for
        // the sends and receives
        IsideBufferSize=0;
        for( size_t s=0; s<numSend; s++ ) 
            IsideBufferSize+=IsideIndirect[s].size();
        JsideBufferSize=0;
        for( size_t r=0; r<numRecv; r++ ) 
            JsideBufferSize+=JsideIndirect[r].size();

        Check( JsideConnects.size() == numRecv );
        Check( JsideConnects.size() == JsideIndirect.size() );

        // Lastly save the mappings for the scatterList and gatherList
        // operations
        std::vector< std::vector<int> > Bmap(onProcRange,std::vector<int>(0));
        for( size_t r=0, k=0; r<numRecv; r++ )
        {
            for( int l=0, numBuf = JsideIndirect[r].size(); l<numBuf; l++ )
            {
                Bmap[JsideIndirect[r][l]].push_back(k++);
            }
        }
          
        // flatten the loops
        countsList.resize(onProcRange);
        BmapList.resize(JsideBufferSize);
        for( size_t j=0, k=0; j<onProcRange; j++ )
        {
            countsList[j]=Bmap[j].size();
            for( size_t l=0; l<Bmap[j].size(); l++ )
            {
                BmapList[k++]=Bmap[j][l];
            }
        }
#endif // C4_MPI
    }
    return;
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// end of TigsTrace.cc
//---------------------------------------------------------------------------//
