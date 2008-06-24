//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Serial.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:06:25 2002
 * \brief  Serial implementation of C4.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_C4_Serial_hh__
#define __c4_C4_Serial_hh__

#include <c4/config.h>

#ifdef C4_SCALAR

#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "C4_Tags.hh"
#include "ds++/Assert.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// Null source/destination rank
//---------------------------------------------------------------------------//

extern const int proc_null;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

template<class Comm>
void inherit(const Comm &comm)
{
}

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
int send(const T *buffer, 
	 int      size,
	 int      destination, 
	 int      tag)
{
    return C4_SUCCESS;
}

//---------------------------------------------------------------------------//

template<class T>
int receive(T   *buffer, 
	    int  size, 
	    int  source, 
	    int  tag)
{
    return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
C4_Req send_async(const T *buffer, 
		  int      size, 
		  int      destination, 
		  int      tag)
{
    // make a c4 request handle
    C4_Req request;
    return request;
}

//---------------------------------------------------------------------------//

template<class T>
void send_async(C4_Req  &request, 
		const T *buffer, 
		int      size, 
		int      destination,
		int      tag)
{
    Require (!request.inuse());
}

//---------------------------------------------------------------------------//

template<class T>
C4_Req receive_async(T   *buffer, 
		     int  size, 
		     int  source, 
		     int  tag)
{
    // make a c4 request handle
    C4_Req request;
    return request;
}

//---------------------------------------------------------------------------//

template<class T>
void receive_async(C4_Req &request, 
		   T      *buffer, 
		   int     size, 
		   int     source, 
		   int     tag)
{
    Require (!request.inuse());
}

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template<class T>
int broadcast(T *buffer, 
	      int      size,
	      int      root)
{
    return C4_SUCCESS;
}

template<class ForwardIterator, class OutputIterator>
void broadcast(ForwardIterator first,
	       ForwardIterator last,
	       OutputIterator  result)
{
    // No communication needed for Serial use.    
    return;
}

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template<class T> 
void global_sum(T &x)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_prod(T &x)
{
}

//---------------------------------------------------------------------------//

template<class T> 
void global_min(T &x)
{
}

//---------------------------------------------------------------------------//

template<class T> 
void global_max(T &x)
{   
}

//---------------------------------------------------------------------------//

template<class T> 
void global_sum(T *x, int n)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_prod(T *x, int n)
{
}

//---------------------------------------------------------------------------//

template<class T> 
void global_min(T *x, int n)
{
}

//---------------------------------------------------------------------------//

template<class T> 
void global_max(T *x, int n)
{
}

} // end namespace rtt_c4

#endif // C4_SCALAR

#endif  // __c4_C4_Serial_hh__

//---------------------------------------------------------------------------//
//                              end of c4/C4_Serial.hh
//---------------------------------------------------------------------------//
