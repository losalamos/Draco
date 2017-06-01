//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Serial.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:06:25 2002
 * \brief  Serial implementation of C4.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_C4_Serial_hh__
#define __c4_C4_Serial_hh__

#include "c4/config.h"

#ifdef C4_SCALAR

#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "C4_Tags.hh"
#include "ds++/Assert.hh"
#include <algorithm>

namespace rtt_c4 {

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

template <class Comm> void inherit(const Comm & /*comm*/) {}

template <class T>
int create_vector_type(unsigned /*count*/, unsigned /*blocklength*/,
                       unsigned /*stride*/, C4_Datatype & /*new_type*/) {
  return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template <class T>
int send(const T * /* buffer*/, int /* size */, int /* destination */,
         int /* tag */) {
  return C4_SUCCESS;
}

template <class T>
int send(const T * /* buffer*/, int /* size */, int /* destination */,
         C4_Datatype & /*data_type*/, int /* tag */) {
  return C4_SUCCESS;
}

template <class T>
DLL_PUBLIC_c4 void send_is(C4_Req & /*request*/, const T * /*buffer*/,
                           int /* size*/, int /*destination*/, int /*tag*/) {
  Insist(false, "send_is is not support for C4_SCALAR builds.");
}

template <typename T>
DLL_PUBLIC_c4 int send_udt(const T * /*buffer*/, int /*size*/,
                           int /*destination*/, C4_Datatype & /*data_type*/,
                           int /*tag*/) {
  return C4_SUCCESS;
}

template <typename TS, typename TR>
DLL_PUBLIC_c4 int send_receive(TS * /*sendbuf*/, int /*sendcount*/,
                               int /*destination*/, TR * /*recvbuf*/,
                               int /*recvcount*/, int /*source*/,
                               int /*sendtag*/, int /*recvtag*/) {
  Insist(false, "send_receive is not support for C4_SCALAR builds.");
  return 1;
}

//---------------------------------------------------------------------------//

template <class T>
int receive(T * /* buffer */, int /* size */, int /* source */, int /* tag */) {
  return C4_SUCCESS;
}

template <class T>
int receive(T * /* buffer */, int /* size */, int /* source */,
            C4_Datatype & /*data_type*/, int /* tag */) {
  return C4_SUCCESS;
}

template <typename T>
DLL_PUBLIC_c4 int receive_udt(T * /*buffer*/, int /*size*/, int /*destination*/,
                              C4_Datatype & /*data_type*/, int /*tag*/) {
  return C4_SUCCESS;
}

template <typename T> DLL_PUBLIC_c4 T prefix_sum(T &node_value) {
  return node_value;
}

//---------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template <class T>
C4_Req send_async(T const * /* buffer */, int /* size */, int /* destination */,
                  int /* tag */) {
  // make a c4 request handle
  C4_Req request;
  return request;
}

//---------------------------------------------------------------------------//

template <class T>
DLL_PUBLIC_c4 void send_async(C4_Req &Remember(request), T const * /* buffer */,
                              int /* size   */, int /* destination */,
                              int /* tag    */) {
  Require(!request.inuse());
}

//---------------------------------------------------------------------------//

template <class T>
C4_Req receive_async(T * /*buffer*/, int /*size  */, int /*source*/,
                     int /*tag   */) {
  // make a c4 request handle
  C4_Req request;
  return request;
}

//---------------------------------------------------------------------------//

template <class T>
void receive_async(C4_Req &Remember(request), T * /* buffer  */,
                   int /* size    */, int /* source  */, int /* tag     */) {
  Require(!request.inuse());
}

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template <class T>
int broadcast(T * /* buffer */, int /* size   */, int /* root   */) {
  return C4_SUCCESS;
}

template <class ForwardIterator, class OutputIterator>
void broadcast(ForwardIterator /* first  */, ForwardIterator /* last   */,
               OutputIterator /* result */) {
  // No communication needed for Serial use.
  return;
}

// safer version of broadcast using stl ranges
template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator /*first*/, ForwardIterator /*last*/,
               OutputIterator /*result*/, OutputIterator /*result_end*/) {
  // No communication needed for Serial use.
  return;
}

//---------------------------------------------------------------------------//
// GATHER/SCATTER
//---------------------------------------------------------------------------//

template <class T> int gather(T *send_buffer, T *receive_buffer, int size) {
  std::copy(send_buffer, send_buffer + size, receive_buffer);
  return C4_SUCCESS;
}

template <class T> int allgather(T *send_buffer, T *receive_buffer, int size) {
  std::copy(send_buffer, send_buffer + size, receive_buffer);
  return C4_SUCCESS;
}

template <class T> int scatter(T *send_buffer, T *receive_buffer, int size) {
  std::copy(send_buffer, send_buffer + size, receive_buffer);
  return C4_SUCCESS;
}

template <class T>
int gatherv(T *send_buffer, int send_size, T *receive_buffer,
            int *receive_sizes, int *receive_displs) {
  std::copy(send_buffer, send_buffer + send_size,
            receive_buffer + receive_displs[0]);

  return C4_SUCCESS;
}

template <class T>
int scatterv(T *send_buffer, int *send_sizes, int *send_displs,
             T *receive_buffer, int receive_size) {
  std::copy(send_buffer + send_displs[0],
            send_buffer + send_displs[0] + send_sizes[0], receive_buffer);

  return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template <class T> DLL_PUBLIC_c4 void global_sum(T & /*x*/) { /* empty */
}

//---------------------------------------------------------------------------//

template <class T> DLL_PUBLIC_c4 void global_prod(T & /*x*/) { /* empty */
}

//---------------------------------------------------------------------------//

template <class T> DLL_PUBLIC_c4 void global_min(T & /*x*/) { /* empty */
}

//---------------------------------------------------------------------------//

template <class T> DLL_PUBLIC_c4 void global_max(T & /*x*/) { /* empty */
}

//---------------------------------------------------------------------------//

template <class T>
DLL_PUBLIC_c4 void global_sum(T * /*x*/, int /*n*/) { /* empty */
}

//---------------------------------------------------------------------------//

template <class T>
DLL_PUBLIC_c4 void global_prod(T * /*x*/, int /*n*/) { /* empty */
}

//---------------------------------------------------------------------------//

template <class T>
DLL_PUBLIC_c4 void global_min(T * /*x*/, int /*n*/) { /* empty */
}

//---------------------------------------------------------------------------//

template <class T>
DLL_PUBLIC_c4 void global_max(T * /*x*/, int /*n*/) { /* empty */
}

} // end namespace rtt_c4

#endif // C4_SCALAR

#endif // __c4_C4_Serial_hh__

//---------------------------------------------------------------------------//
// end of c4/C4_Serial.hh
//---------------------------------------------------------------------------//
