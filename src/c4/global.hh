//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/global.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 10:56:16 2002
 * \brief  C4 function declarations and class * definitions.
 * \note   Copyright (C) 2002-2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * This file allows the client to include the message passing services
 * provided by C4.  The function declarations and class definitions are
 * contained in the rtt_c4 namespace.  For backwards compatibility, the
 * old-style C4 functions and classes are declared in the C4 namespace.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_c4_global_hh
#define rtt_c4_global_hh

// C4 Message Passing Functions
#include "C4_Functions.hh" 

//---------------------------------------------------------------------------//
// Include the appropriate header for an underlying message passing
// implementation.  This allows the definition of inline functions declared
// in C4_Functions.hh.

#ifdef C4_SCALAR
#include "C4_Serial.hh"
#endif

#ifdef C4_MPI
#include "C4_MPI.hh"
#endif

//===========================================================================//
/*!
 * \namespace C4
 *
 * \brief Deprecated namespace that contains the C4 package.
 *
 * This namespace contains all of C4 <= 1_6_0.  c4-2_0_0 is written in the
 * rtt_c4 namespace.  However, for porting purposes, the functionality from
 * c4-1_6_0 and back is preserved in the C4 namespace.  This is accomplished
 * by introducing the new c4 names into the old namespace.  Only functions
 * and services from c4-1_6_0 and earlier are in this deprecated namespace.
 * New services are defined only in the rtt_c4 namespace.
 */
//===========================================================================//

namespace C4
{

using rtt_c4::node;
using rtt_c4::nodes;
using rtt_c4::C4_Req;

//---------------------------------------------------------------------------//
// Backwards compatibility functions

inline void Init(int &argc, char **&argv)
{
    rtt_c4::initialize(argc, argv);
}

inline void Finalize()
{
    rtt_c4::finalize();
}

inline void gsync()
{
    rtt_c4::global_barrier();
}

template<class T>
inline int Send(const T *buf, int nels, int dest,
		int tag = rtt_c4::C4_Traits<T*>::tag, int /*group*/ = 0)
{
    return rtt_c4::send(buf, nels, dest, tag);
}

template<class T>
inline int Recv(T *buf, int nels, int source,
		int tag = rtt_c4::C4_Traits<T*>::tag, int /*group*/ = 0)
{
    return rtt_c4::receive(buf, nels, source, tag);
}

// When using either of these functions, you may need to specify the template
// type explicitly.  For example:
//
//     int myInt(rtt_c4::node()+10);
//     if( rtt_c4::node() == 0 )  C4::Send<int>(myInt,1);
//     if( rtt_c4::node() == 1 )  C4::Recv<int>(myInt,0);

// template<class T>
// inline int Send(T data, int destination,
// 		int tag = rtt_c4::C4_Traits<T>::tag, int group = 0)
// {
//     return rtt_c4::send(&data, 1, destination, tag);
// }

// template<class T>
// inline int Recv(T &data, int source,
// 		int tag = rtt_c4::C4_Traits<T>::tag, int group = 0)
// {
//     return rtt_c4::receive(&data, 1, source, tag);
// }

template<class T>
inline void SendAsync(C4_Req& r, const T *buf, int nels, int dest,
		      int tag = rtt_c4::C4_Traits<T*>::tag, int /*group*/ = 0)
{
    rtt_c4::send_async(r, buf, nels, dest, tag);
}

template<class T>
inline void RecvAsync(C4_Req& r, T *buf, int nels, int source,
		      int tag = rtt_c4::C4_Traits<T*>::tag, int /*group*/ = 0)
{
    rtt_c4::receive_async(r, buf, nels, source, tag);
}

template<class T>
inline C4_Req SendAsync(const T *buf, int nels, int dest,
			int tag = rtt_c4::C4_Traits<T*>::tag, int /*group*/ = 0)
{
    return rtt_c4::send_async(buf, nels, dest, tag); 
}

template<class T>
inline C4_Req RecvAsync(T *buf, int nels, int source,
			int tag = rtt_c4::C4_Traits<T*>::tag, int /*group*/ = 0)
{
    return rtt_c4::receive_async(buf, nels, source, tag);
}

template<class T>
inline void gsum(T &x)
{
    rtt_c4::global_sum(x);
    return;
}

template<class T>
inline void gprod(T &x)
{
    rtt_c4::global_prod(x);
    return;
}

template<class T>
inline void gmin(T &x)
{
    rtt_c4::global_min(x);
    return;
}

template<class T>
inline void gmax(T &x)
{
    rtt_c4::global_max(x);
    return;
}

template<class T>
inline void gsum(T *x, int n)
{
    rtt_c4::global_sum(x, n);
    return;
}

template<class T>
inline void gprod(T *x, int n)
{
    rtt_c4::global_prod(x, n);
    return;
}

template<class T>
inline void gmin(T *x, int n)
{
    rtt_c4::global_min(x, n);
    return;
}

template<class T>
inline void gmax(T *x, int n)
{
    rtt_c4::global_max(x, n);
    return;
}

inline double Wtime()
{
    return rtt_c4::wall_clock_time();
}

inline double Wtick()
{
    return rtt_c4::wall_clock_resolution();
}

} // end of namespace C4

#endif                          // rtt_c4_global_hh

//---------------------------------------------------------------------------//
//                              end of c4/global.hh
//---------------------------------------------------------------------------//
