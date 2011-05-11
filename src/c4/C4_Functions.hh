//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Functions.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:42:03 2002
 * \brief  C4 Communication Functions.
 *
 * This file contains the declarations for communication functions provided
 * by C4.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_C4_Functions_hh
#define c4_C4_Functions_hh

#include "ds++/config.h"
#include "C4_sys_times.h"
#include "C4_Traits.hh"
#include "C4_Req.hh"

namespace rtt_c4
{

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
 
// Forward declarations.
class C4_Req;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a parallel job.
 */
void initialize(int &argc, char **&argv);

//---------------------------------------------------------------------------//
/*!
 * \brief Finish a parallel job.
 */
void finalize();

//---------------------------------------------------------------------------//
/*!
 * \brief Inherit a communicator from another application.
 */
template<class Comm>
void inherit(const Comm &);

//---------------------------------------------------------------------------//
/*!
 * \brief Free an inherited communicator from another application.
 */
void free_inherited_comm();

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Get the node (rank) of the current processor.
 *
 * The rank is determined by the current communicator.
 */
int node();

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of processors used for this job.
 *
 * The number of nodes is determined by the current communicator.
 */
int nodes();

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set a global barrier for the communicator.
 */
void global_barrier();

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, blocking send.
 */
template<class T>
int send(const T *buffer, int size, int destination,
	 int tag = C4_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, blocking receive.
 */
template<class T>
int receive(T *buffer, int size, int source, int tag = C4_Traits<T*>::tag);

//---------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 *
 * \return C4_Req object to handle communciation requests
 */
template<class T>
C4_Req send_async( T const * buffer,
                   int       size,
                   int       destination,
                   int       tag );

template<class T>
C4_Req send_async( T const * buffer,
                   int       size,
                   int       destination )
{
    int tag =  C4_Traits<T*>::tag;
    return send_async( buffer, size, destination, tag );
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced 
//    function template is nonstandard 

// template<class T>
// C4_Req send_async( T const * buffer,
//                    int       size,
//                    int       destination,
//                    int       tag = C4_Traits<T*>::tag );

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 */
template<class T>
void send_async( C4_Req       & request,
                 T      const * buffer,
                 int            size,
                 int            destination,
                 int            tag );
template<class T>
void send_async( C4_Req       & request,
                 T      const * buffer,
                 int            size,
                 int            destination )
{
    int tag = C4_Traits<T*>::tag;
    send_async( request, buffer, size, destination, tag );
    return;
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced 
//    function template is nonstandard 

// template<class T>
// void send_async( C4_Req       & request,
//                  T      const * buffer,
//                  int            size,
//                  int            destination,
//                  int            tag = C4_Traits<T*>::tag );

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking synchronous send.
 */
template<class T>
void send_is(C4_Req       & request,
             T      const * buffer,
             int            size,
             int            destination,
             int            tag );
template<class T>
void send_is(C4_Req       & request,
             T      const * buffer,
             int            size,
             int            destination )
{
    int tag = C4_Traits<T*>::tag;
    send_is( request, buffer, size, destination, tag );
    return;
}

// [2011-05-11 GMR] This declaration should replace the two preceeding ones.
// However, I expect that PGI-10 doesn't like this syntax for the same reason
// that it doesn't like the syntax for send_async above.

// template<class T>
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

template<class T>
C4_Req receive_async( T   * buffer,
                      int   size,
                      int   source,
                      int   tag );
template<class T>
C4_Req receive_async( T   * buffer,
                      int   size,
                      int   source )
{
    int tag = C4_Traits<T*>::tag;
    return receive_async<T>( buffer, size, source, tag );
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced 
//    function template is nonstandard 

// template<class T>
// C4_Req receive_async( T   * buffer,
//                       int   size,
//                       int   source,
//                       int   tag = C4_Traits<T*>::tag );

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 */
template<class T>
void receive_async( C4_Req & request,
                    T      * buffer,
                    int      size,
                    int      source,
                    int      tag );
template<class T>
void receive_async( C4_Req & request,
                    T      * buffer,
                    int      size,
                    int      source)
{
    int tag = C4_Traits<T*>::tag;
    receive_async( request, buffer, size, source, tag );
    return;
}

// [2010-07-22 KT] This declaration should replace the two preceeding ones.
// However, PGI-10 doesn't like this syntax and issues the warning:
//    error: specifying a default argument when redeclaring an unreferenced 
//    function template is nonstandard 

// template<class T>
// void receive_async( C4_Req & request,
//                     T      * buffer,
//                     int      size,
//                     int      source,
//                     int      tag = C4_Traits<T*>::tag );
//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template<class T>
int broadcast(T *buffer, int size, int root);


/*---------------------------------------------------------------------------*/
/*! 
 * \brief Send data from processor 0 to all other processors.
 */
template<class ForwardIterator, class OutputIterator>
void broadcast(ForwardIterator first,
	       ForwardIterator last,
	       OutputIterator result);

//---------------------------------------------------------------------------//
// GATHER/SCATTER
//---------------------------------------------------------------------------//

template<class T>
int gather(T *send_buffer, T *receive_buffer, int size);

template<class T>
int gatherv(T *send_buffer,
            int send_size,
            T *receive_buffer,
            int *receive_sizes,
            int *receive_displs);

template<class T>
int scatter(T *send_buffer, T *receive_buffer, int size);

template<class T>
int scatterv(T *send_buffer,
             int *send_sizes,
             int *send_displs,
             T *receive_buffer,
             int receive_size);

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a global sum of a scalar variable.
 */
template<class T> 
void global_sum(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global product of a scalar variable.
 */
template<class T>
void global_prod(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global minimum of a scalar variable.
 */
template<class T> 
void global_min(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global maximum of a scalar variable.
 */
template<class T> 
void global_max(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global sum of an array.
 */
template<class T> 
void global_sum(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global product of an array.
 */
template<class T>
void global_prod(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global minimum of an array.
 */
template<class T> 
void global_min(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global maximum of an array.
 */
template<class T> 
void global_max(T *x, int n);

//---------------------------------------------------------------------------//
// TIMING FUNCTIONS
//---------------------------------------------------------------------------//
/*! 
 * \brief Return the wall-clock time in seconds.
 */
double wall_clock_time();
double wall_clock_time( DRACO_TIME_TYPE & now );


//---------------------------------------------------------------------------//
/*!
 * \brief Return the resolution of wall_clock_time.
 */
double wall_clock_resolution();

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
bool probe(int source, int tag, int &message_size);

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
void blocking_probe(int source, int tag, int &message_size);

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
void wait_all(int count, C4_Req *requests);

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
unsigned wait_any(int count, C4_Req *requests);

//---------------------------------------------------------------------------//
// ABORT
//---------------------------------------------------------------------------//
/*!
 * \brief Abort across all processors.
 *
 * \param error suggested return error, defaults to 1
 */
int abort(int error = 1);

//---------------------------------------------------------------------------//
// isScalar
//---------------------------------------------------------------------------//
/*!
 * \brief Is C4 executing in scalar-only mode?
 */
bool isScalar();

} // end namespace rtt_c4

#endif                          // c4_C4_Functions_hh

//---------------------------------------------------------------------------//
//                              end of c4/C4_Functions.hh
//---------------------------------------------------------------------------//
