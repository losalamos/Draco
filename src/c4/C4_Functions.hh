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

#include <unistd.h>
#include <sys/times.h>  // defines the struct tms.
#include "C4_Traits.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
/*!
 * C4 unit tests.
 *
 * \example c4/test/tstAbort.cc
 * \example c4/test/tstBroadcast.cc
 * \example c4/test/tstComm_Dup.cc
 * \example c4/test/tstPingPong.cc
 * \example c4/test/tstReduction.cc
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
C4_Req send_async(const T *buffer, int size, int destination,
		  int tag = C4_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 */
template<class T>
void send_async(C4_Req &request, const T *buffer, int size, int destination,
		int tag = C4_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 *
 * \return C4_Req object to handle communciation requests
 */
template<class T>
C4_Req receive_async(T *buffer, int size, int source,
		     int tag = C4_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 */
template<class T>
void receive_async(C4_Req& request, T *buffer, int size, int source,
		   int tag = C4_Traits<T*>::tag);

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
double wall_clock_time( tms & now );

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
// ABORT
//---------------------------------------------------------------------------//
/*!
 * \brief Abort across all processors.
 *
 * \param error suggested return error, defaults to 1
 */
int abort(int error = 1);

} // end namespace rtt_c4

#endif                          // c4_C4_Functions_hh

//---------------------------------------------------------------------------//
//                              end of c4/C4_Functions.hh
//---------------------------------------------------------------------------//
