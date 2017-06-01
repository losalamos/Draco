//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Functions.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:42:03 2002
 * \brief  C4 Communication Functions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * This file contains the declarations for communication functions provided by
 * C4. This file allows the client to include the message passing services
 * provided by C4.  The function declarations and class definitions are
 * contained in the rtt_c4 namespace.  For backwards compatibility, the
 * old-style C4 functions and classes are declared in the C4 namespace.
 */
//---------------------------------------------------------------------------//

#ifndef c4_C4_Functions_hh
#define c4_C4_Functions_hh

#include "C4_Datatype.hh"
#include "C4_Req.hh"
#include "C4_Status.hh"
#include "C4_Traits.hh"
#include "C4_sys_times.h"

namespace rtt_c4 {

//---------------------------------------------------------------------------//
/*!
 * C4 unit tests.
 */
/*! \example c4/test/tstAbort.cc
 * Example of MPI abort functions.
 */
/*! \example c4/test/tstBroadcast.cc
 * Example of MPI broadcast-like functions
 */
/*! \example c4/test/tstComm_Dup.cc
 * Example
 */
/*! \example c4/test/tstPingPong.cc
 * Example of point-to-point communications
 */
/*! \example c4/test/tstReduction.cc
 * Example of many-to-one communications
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// GLOBAL CONSTANTS
//---------------------------------------------------------------------------//

//! Any source rank
DLL_PUBLIC_c4 extern const int any_source;
//! Null source/destination rank
DLL_PUBLIC_c4 extern const int proc_null;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a parallel job.
 */
DLL_PUBLIC_c4 int initialize(int &argc, char **&argv,
                             int required = DRACO_MPI_THREAD_SINGLE);

//---------------------------------------------------------------------------//
/*!
 * \brief Finish a parallel job.
 */
DLL_PUBLIC_c4 void finalize();

//---------------------------------------------------------------------------//
/*!
 * \brief Inherit a communicator from another application.
 */
template <class Comm> void inherit(const Comm &);

//---------------------------------------------------------------------------//
/*!
 * \brief Free an inherited communicator from another application.
 */
DLL_PUBLIC_c4 void free_inherited_comm();

//---------------------------------------------------------------------------//
/*!
 * \brief Create up a new vector type.
 *
 * \param count Number of blocks in the data type
 * \param blocklength Length of each block (in units of base type)
 * \param stride Spacing between start of each block (in units of base type)
 * \param new_type On return, contains the new type descriptor.
 */
template <class T>
int create_vector_type(unsigned count, unsigned blocklength, unsigned stride,
                       C4_Datatype &new_type);

//---------------------------------------------------------------------------//
//! Free a user defined type, such as a vector type.

DLL_PUBLIC_c4 void type_free(C4_Datatype &old_type);

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Get the node (rank) of the current processor.
 *
 * The rank is determined by the current communicator.
 */
DLL_PUBLIC_c4 int node();

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of processors used for this job.
 *
 * The number of nodes is determined by the current communicator.
 */
DLL_PUBLIC_c4 int nodes();

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set a global barrier for the communicator.
 */
DLL_PUBLIC_c4 void global_barrier();

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//
//! Do a point-to-point, blocking send.
template <typename T>
DLL_PUBLIC_c4 int send(const T *buffer, int size, int destination,
                       int tag = C4_Traits<T *>::tag);

//---------------------------------------------------------------------------//
//! Do a point-to-point, blocking receive.
template <typename T>
DLL_PUBLIC_c4 int receive(T *buffer, int size, int source,
                          int tag = C4_Traits<T *>::tag);

//---------------------------------------------------------------------------//
//! Do a point-to-point, blocking send of a user-defined type.
template <typename T>
DLL_PUBLIC_c4 int send_udt(const T *buffer, int size, int destination,
                           C4_Datatype &, int tag = C4_Traits<T *>::tag);

//---------------------------------------------------------------------------//
//! Do a point-to-point, blocking send-receive.
template <typename TS, typename TR>
DLL_PUBLIC_c4 int send_receive(TS *sendbuf, int sendcount, int destination,
                               TR *recvbuf, int recvcount, int source,
                               int sendtag = C4_Traits<TS *>::tag,
                               int recvtag = C4_Traits<TR *>::tag);

//---------------------------------------------------------------------------//
//! Do a point-to-point, blocking receive of a user-defined type.
template <typename T>
DLL_PUBLIC_c4 int receive_udt(T *buffer, int size, int source, C4_Datatype &,
                              int tag = C4_Traits<T *>::tag);

//---------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 *
 * \return C4_Req object to handle communciation requests
 */
template <typename T>
DLL_PUBLIC_c4 C4_Req send_async(T const *buffer, int size, int destination,
                                int tag);

template <typename T>
C4_Req send_async(T const *buffer, int size, int destination) {
  int tag = C4_Traits<T *>::tag;
  return send_async(buffer, size, destination, tag);
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced
//    function template is nonstandard

// template<typename T>
// C4_Req send_async( T const * buffer,
//                    int       size,
//                    int       destination,
//                    int       tag = C4_Traits<T*>::tag );

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 */
template <typename T>
DLL_PUBLIC_c4 void send_async(C4_Req &request, T const *buffer, int size,
                              int destination, int tag);
template <typename T>
void send_async(C4_Req &request, T const *buffer, int size, int destination) {
  int tag = C4_Traits<T *>::tag;
  send_async(request, buffer, size, destination, tag);
  return;
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced
//    function template is nonstandard

// template<typename T>
// void send_async( C4_Req       & request,
//                  T      const * buffer,
//                  int            size,
//                  int            destination,
//                  int            tag = C4_Traits<T*>::tag );

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking synchronous send.
 */
template <typename T>
DLL_PUBLIC_c4 void send_is(C4_Req &request, T const *buffer, int size,
                           int destination, int tag);
template <typename T>
void send_is(C4_Req &request, T const *buffer, int size, int destination) {
  int tag = C4_Traits<T *>::tag;
  send_is(request, buffer, size, destination, tag);
  return;
}

// [2011-05-11 GMR] This declaration should replace the two preceeding ones.
// However, I expect that PGI-10 doesn't like this syntax for the same reason
// that it doesn't like the syntax for send_async above.

// template<typename T>
// void send_is(C4_Req       & request,
//              T      const * buffer,
//              int            size,
//              int            destination,
//              int            tag = C4_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 *
 * \return C4_Req object to handle communciation requests
 */

template <typename T>
C4_Req receive_async(T *buffer, int size, int source, int tag);
template <typename T> C4_Req receive_async(T *buffer, int size, int source) {
  int tag = C4_Traits<T *>::tag;
  return receive_async<T>(buffer, size, source, tag);
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced
//    function template is nonstandard

// template<typename T>
// C4_Req receive_async( T   * buffer,
//                       int   size,
//                       int   source,
//                       int   tag = C4_Traits<T*>::tag );

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 */
template <typename T>
DLL_PUBLIC_c4 void receive_async(C4_Req &request, T *buffer, int size,
                                 int source, int tag);
template <typename T>
void receive_async(C4_Req &request, T *buffer, int size, int source) {
  int tag = C4_Traits<T *>::tag;
  receive_async(request, buffer, size, source, tag);
  return;
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced
//    function template is nonstandard

// template<typename T>
// void receive_async( C4_Req & request,
//                     T      * buffer,
//                     int      size,
//                     int      source,
//                     int      tag = C4_Traits<T*>::tag );
//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template <typename T>
DLL_PUBLIC_c4 int broadcast(T *buffer, int size, int root);

/*---------------------------------------------------------------------------*/
/*!
 * \brief Send data from processor 0 to all other processors.
 */
template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator first, ForwardIterator last,
               OutputIterator result);

//---------------------------------------------------------------------------//
// GATHER/SCATTER
//---------------------------------------------------------------------------//

template <typename T>
DLL_PUBLIC_c4 int gather(T *send_buffer, T *receive_buffer, int size);

template <typename T>
DLL_PUBLIC_c4 int allgather(T *send_buffer, T *receive_buffer, int size);

template <typename T>
int gatherv(T *send_buffer, int send_size, T *receive_buffer,
            int *receive_sizes, int *receive_displs);

template <typename T>
DLL_PUBLIC_c4 int scatter(T *send_buffer, T *receive_buffer, int size);

template <typename T>
int scatterv(T *send_buffer, int *send_sizes, int *send_displs,
             T *receive_buffer, int receive_size);

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a global sum of a scalar variable.
 */
template <typename T> DLL_PUBLIC_c4 void global_sum(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global product of a scalar variable.
 */
template <typename T> DLL_PUBLIC_c4 void global_prod(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global minimum of a scalar variable.
 */
template <typename T> DLL_PUBLIC_c4 void global_min(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global maximum of a scalar variable.
 */
template <typename T> DLL_PUBLIC_c4 void global_max(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global sum of an array.
 */
template <typename T> DLL_PUBLIC_c4 void global_sum(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global product of an array.
 */
template <typename T> DLL_PUBLIC_c4 void global_prod(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global minimum of an array.
 */
template <typename T> DLL_PUBLIC_c4 void global_min(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global maximum of an array.
 */
template <typename T> DLL_PUBLIC_c4 void global_max(T *x, int n);

//---------------------------------------------------------------------------//
// TIMING FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the wall-clock time in seconds.
 */
DLL_PUBLIC_c4 double wall_clock_time();
double wall_clock_time(DRACO_TIME_TYPE &now);

//---------------------------------------------------------------------------//
/*!
 * \brief Return the resolution of wall_clock_time.
 */
DLL_PUBLIC_c4 double wall_clock_resolution();

//---------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief See if a message is pending.
 *
 * \param source
 * Processor from which a message may be pending.
 * \param tag
 * Tag for pending message.
 * \param message_size
 * On return, size of the pending message in bytes.
 * \return \c true if a message from the specified processor with the
 * specified tag is pending; \c false otherwise.
 */
DLL_PUBLIC_c4 bool probe(int source, int tag, int &message_size);

//---------------------------------------------------------------------------//
/*!
 * \brief Wait until a message (of unknown size) is pending.
 *
 * \param source
 * Processor from which a message of unknown size is expected.
 * \param tag
 * Tag for pending message.
 * \param message_size
 * On return, size of the pending message in bytes.
 */
DLL_PUBLIC_c4 void blocking_probe(int source, int tag, int &message_size);

//---------------------------------------------------------------------------//
/*!
 * \brief Wait until every one of a set of posted sends/receives is complete.
 *
 * This version returns no status information.
 *
 * \param count
 * Size of the set of requests to wait on.
 * \param requests
 * Set of requests to wait on.
 */
DLL_PUBLIC_c4 void wait_all(unsigned count, C4_Req *requests);

//---------------------------------------------------------------------------//
/*!
 * \brief Wait until one of a set of posted sends/receives is complete.
 *
 * \param count
 * Size of the set of requests to wait on.
 * \param requests
 * Set of requests to wait on.
 * \return The request that completed.
 */
DLL_PUBLIC_c4 unsigned wait_any(unsigned count, C4_Req *requests);

//---------------------------------------------------------------------------//
// ABORT
//---------------------------------------------------------------------------//
/*!
 * \brief Abort across all processors.
 *
 * \param error suggested return error, defaults to 1
 */
DLL_PUBLIC_c4 int abort(int error = 1);

//---------------------------------------------------------------------------//
// isScalar
//---------------------------------------------------------------------------//
/*!
 * \brief Is C4 executing in scalar-only mode?
 */
DLL_PUBLIC_c4 bool isScalar();

//---------------------------------------------------------------------------//
// get_processor_name
//---------------------------------------------------------------------------//
//! Return the processor name for each rank.
DLL_PUBLIC_c4 std::string get_processor_name();

//---------------------------------------------------------------------------//
// prefix_sum
//---------------------------------------------------------------------------//
/*!
 * \brief Return the value of the prefix sum at this processor.
 *
 * \param node_value Current node's value of variable to be prefix summed
 * \return Sum of value over nodes up to and including this node.
 */
template <typename T> DLL_PUBLIC_c4 T prefix_sum(const T node_value);

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// Include the appropriate header for an underlying message passing
// implementation.

#ifdef C4_SCALAR
#include "C4_Serial.hh"
#endif

#ifdef C4_MPI
#include "C4_MPI.hh"
#endif

#endif // c4_C4_Functions_hh

//---------------------------------------------------------------------------//
// end of c4/C4_Functions.hh
//---------------------------------------------------------------------------//
