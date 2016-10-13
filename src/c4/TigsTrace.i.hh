//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/TigsTrace.i.hh
 * \author Bob Webster, Thomas M. Evans
 * \date   Fri Sep 30 12:53:25 2005
 * \brief  Member definitions of class TigsTrace
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Tigs.i.hh 6056 2012-06-19 19:05:27Z kellyt $
//---------------------------------------------------------------------------//

#ifndef c4_TigsTrace_i_hh
#define c4_TigsTrace_i_hh

#include "C4_Functions.hh"

namespace rtt_c4 {

//---------------------------------------------------------------------------//
// COMM::TRACE PUBLIC TEMPLATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Copy data defined on the domain (I) to the range (J).
 *
 * This is the inverse operation of gatherList
 *
 * Recall that the trace is defined by a map JMap:I->J given to the trace
 * constructor as a container of ints that with a length equal to the local
 * length of the domain (I) space. Each of those ints is the Global index of
 * the range (J) location to which the local I location is mapped.
 *
 * This function takes a container A[I] (specified by its begining and ending
 * iterators) and copies the data to a container B in the following CSR like
 * manner.
 *
 * \verbatim
 *  B= { (A[i1[1]],...,A[iN[1]]), (A[i1[2]], ... , A[iN[2])
 *
 *     ... (A[i1[j]], ... ,A[iN[j]]), ... , (A[i1[J]], ... ,A[iN[J] }
 * \endverbatim
 *
 * where i1[j] is the first value of i (globally) such that JMap[i1]=j.  The
 * number of i's that map to each j (N above) is returned in the container Cnt
 * (again specified through iterators.)
 *
 * \param Afirst    beginning of the source data
 * \param Alast     end of the source data
 * \param Cntfirst  beginning of the count list
 * \param Cntlast   end of the count list
 * \param Bfirst    beginning of the destination data
 * \param Blast     end of the destination data
 */
template <typename iterA, typename iterC, typename iterB>
void TigsTrace::scatterList(iterA Afirst, iterA Remember(Alast), iterC Cntfirst,
                            iterC Remember(Cntlast), iterB Bfirst,
                            iterB Remember(Blast)) const {
  Require(JsideIndirect.size() == JsideConnects.size());
  Require(IsideIndirect.size() == IsideConnects.size());

  Require(std::distance(Afirst, Alast) == onProcDomain);
  Require(std::distance(Cntfirst, Cntlast) == onProcRange);
  Require(std::distance(Bfirst, Blast) == JsideBufferSize);

  if (rtt_c4::nodes() == 1) {
    for (size_t j = 0, k = 0; j < onProcRange; j++) {
      *(Cntfirst + j) = counts[j];
      for (size_t l = 0; l < counts[j]; l++, k++)
        *(Bfirst + k) = *(Afirst + IMV[k]);
    }
  } else {
#ifdef C4_MPI
    // Scatter always moves data from the I side to the Jside
    // The traditional scatter involves a reduction
    // Scatter_list form a CSR like packed representation
    // of the moved data.
    typedef typename std::iterator_traits<iterB>::value_type T;
    size_t const numSend(IsideIndirect.size());
    size_t const numRecv(JsideIndirect.size());
    std::vector<T> sbuffer(IsideBufferSize);
    std::vector<T> rbuffer(JsideBufferSize);
    std::vector<rtt_c4::C4_Req> reqs_recv;

    // fill the send buffer on the Iside
    size_t k = 0;
    for (size_t s = 0; s < IsideConnects.size(); s++) {
      for (size_t l = 0; l < IsideIndirect[s].size(); l++) {
        sbuffer[k++] = *(Afirst + IsideIndirect[s][l]);
      }
    }

    size_t krcv = 0;
    int recv_from_self_index = -1;
    for (size_t r = 0; r < numRecv; r++) {
      // do not post receives-from-self
      if (JsideConnects[r] != rtt_c4::node()) {
        reqs_recv.push_back(rtt_c4::receive_async(
            &rbuffer[krcv], JsideIndirect[r].size(), JsideConnects[r]));
      }
      // store the receive buffer index for use later
      else {
        Check(JsideConnects[r] == rtt_c4::node());
        recv_from_self_index = krcv;
      }
      krcv += JsideIndirect[r].size();
    }
    size_t ksend = 0;
    for (size_t s = 0; s < numSend; s++) {
      // do not send-to-self
      if (IsideConnects[s] != rtt_c4::node()) {
        rtt_c4::C4_Req c4req_send = rtt_c4::send_async(
            &sbuffer[ksend], IsideIndirect[s].size(), IsideConnects[s]);
        c4req_send.wait();
      }
      // otherwise put directly into receive buffer
      else {
        Check(recv_from_self_index >= 0);
        Check(recv_from_self_index < static_cast<int>(JsideBufferSize));
        std::copy(&sbuffer[ksend], &sbuffer[ksend] + IsideIndirect[s].size(),
                  &rbuffer[recv_from_self_index]);
      }
      ksend += IsideIndirect[s].size();
    }

    for (size_t s = 0; s < reqs_recv.size(); ++s)
      reqs_recv[s].wait();

    rtt_c4::global_barrier();

    // unload the receive buffer Here is where we form the
    // CSR packing...
    for (size_t j = 0; j < onProcRange; j++)
      *(Cntfirst + j) = countsList[j];
    for (size_t k = 0; k < JsideBufferSize; k++)
      *(Bfirst + k) = rbuffer[BmapList[k]];

#endif // C4_MPI
  }
  return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy data from the range to the domain.
 * 
 * Recall that the trace is defined by a map JMap:I->J given to the trace
 * constructor as a container of ints that with a length equal to the local
 * length of the domain (I) space. Each of those ints is the Global index of
 * the range (J) location to which the local I location is mapped.
 *
 * This function takes a copies data from the container B (given by its
 * beginning and ending iterators) to a container A using the relation
 * A[i]=B[Jmap[i]] for all local i. Because the map JMap must be many to one
 * (or none) this operation is well defined.
 *
 * \param Bfirst    beginning of the source data
 * \param Blast     end of the source data
 * \param Afirst    beginning of the destination data
 * \param Alast     end of the destination data
 */
template <typename iterA, typename iterB>
void TigsTrace::gather(iterB Bfirst, iterB Remember(Blast), iterA Afirst,
                       iterA Remember(Alast)) {
  Require(JsideIndirect.size() == JsideConnects.size());
  Require(IsideIndirect.size() == IsideConnects.size());

  Require(std::distance(Afirst, Alast) == onProcDomain);
  Require(std::distance(Bfirst, Blast) == onProcRange);
  if (rtt_c4::nodes() == 1) {
    for (size_t j = 0, k = 0; j < onProcRange; j++)
      for (size_t l = 0; l < counts[j]; l++)
        *(Afirst + IMV[k++]) = *(Bfirst + j);
  } else {
#ifdef C4_MPI
    typedef typename std::iterator_traits<iterB>::value_type T;
    size_t const numSend = JsideIndirect.size();
    size_t const numRecv = IsideIndirect.size();
    std::vector<T> sbuffer(JsideBufferSize);
    std::vector<T> rbuffer(IsideBufferSize);
    std::vector<rtt_c4::C4_Req> reqs_recv;

    // fill the send buffer on the Jside
    size_t k = 0;
    for (size_t s = 0; s < numSend; s++) {
      for (size_t l = 0; l < JsideIndirect[s].size(); l++) {
        sbuffer[k++] = *(Bfirst + JsideIndirect[s][l]);
      }
    }

    size_t krcv = 0;
    int recv_from_self_index = -1;
    for (size_t r = 0; r < numRecv; r++) {
      // do not post receives-from-self
      if (IsideConnects[r] != rtt_c4::node()) {
        reqs_recv.push_back(rtt_c4::receive_async(
            &rbuffer[krcv], IsideIndirect[r].size(), IsideConnects[r]));

      }
      // store the receive buffer index for use later
      else {
        Check(IsideConnects[r] == rtt_c4::node());
        recv_from_self_index = krcv;
      }
      krcv += IsideIndirect[r].size();
    }
    size_t ksend = 0;
    for (size_t s = 0; s < numSend; s++) {
      // do not send-to-self
      if (JsideConnects[s] != rtt_c4::node()) {
        rtt_c4::C4_Req c4req_send = rtt_c4::send_async(
            &sbuffer[ksend], JsideIndirect[s].size(), JsideConnects[s]);
        c4req_send.wait();
      }
      // otherwise put directly into receive buffer
      else {
        Check(recv_from_self_index >= 0);
        Check(recv_from_self_index < static_cast<int>(IsideBufferSize));
        std::copy(&sbuffer[ksend], &sbuffer[ksend] + JsideIndirect[s].size(),
                  &rbuffer[recv_from_self_index]);
      }
      ksend += JsideIndirect[s].size();
    }

    for (size_t s = 0; s < reqs_recv.size(); ++s)
      reqs_recv[s].wait();
    rtt_c4::global_barrier();

    // unload the receive buffer
    k = 0;
    for (size_t r = 0; r < numRecv; r++) {
      for (size_t l = 0, numBuf = IsideIndirect[r].size(); l < numBuf; l++) {
        *(Afirst + IsideIndirect[r][l]) = rbuffer[k++];
      }
    }

#endif // C4_MPI
  }
  return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy with reduction of data from the domain to the range.
 *
 * Recall that the trace is defined by a map JMap:I->J given to the trace
 * constructor as a container of ints that with a length equal to the local
 * length of the domain (I) space. Each of those ints is the Global index of
 * the range (J) location to which the local I location is mapped.
 *
 * Because the map JMap is many to one (or none), many locations in the domain
 * may be mapped to the same location in the range. The scatter operation
 * loads the container B (specified by beginning and ending iterators) with
 * the value
 *
 * \verbatim
 * B[j]= B[j] op A[i1[j]] op ... op A[iN[j]]
 * \endverbatim
 * 
 * where i1[j] is the first value of i (globally) such that JMap[i1]=j, and
 * iN[j] is the last such i.  The global number of i's that map to j is N. The
 * operand op must be commutative and associative.
 *
 * \param Afirst    beginning of the source data
 * \param Alast     end of the source data
 * \param Bfirst    beginning of the reduced data
 * \param Blast     end of the reduced data
 * \param op        the binary reduction operation
 */
template <typename iterA, typename iterB, typename BinaryOp>
void TigsTrace::scatter(iterA Afirst, iterA Alast, iterB Bfirst,
                        iterB Remember(Blast), BinaryOp op) {
  Check(std::distance(Afirst, Alast) == onProcDomain);
  Check(std::distance(Bfirst, Blast) == onProcRange);
  if (rtt_c4::nodes() == 1) {
    for (size_t j = 0, k = 0; j < onProcRange; j++)
      for (size_t l = 0; l < counts[j]; l++, k++)
        *(Bfirst + j) = op(*(Bfirst + j), *(Afirst + IMV[k]));
  } else {
    typedef typename std::iterator_traits<iterB>::value_type T;
    std::vector<unsigned> counts_ret(onProcRange);
    std::vector<T> Btmp(JsideBufferSize);
    scatterList(Afirst, Alast, counts_ret.begin(), counts_ret.end(),
                Btmp.begin(), Btmp.end());
    for (size_t j = 0, k = 0; j < onProcRange; j++)
      for (size_t l = 0; l < counts_ret[j]; l++, k++)
        *(Bfirst + j) = op(*(Bfirst + j), Btmp[k]);
  }
  return;
}

} // end namespace rtt_c4

#endif // c4_TigsTrace_i_hh

//---------------------------------------------------------------------------//
// end of c4/Tigs.i.hh
//---------------------------------------------------------------------------//
