//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI_blocking_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:41:05 2002
 * \brief  C4 MPI Blocking Send/Recv instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>

#ifdef C4_MPI

#include "C4_MPI.t.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF BLOCKING SEND/RECEIVE AND BROADCAST
//---------------------------------------------------------------------------//

template int send(const char *, int, int, int);
template int send(const unsigned char *, int, int, int);
template int send(const short *, int, int, int);
template int send(const unsigned short *, int, int, int);
template int send(const int *, int, int, int);
template int send(const unsigned int *, int, int, int);
template int send(const long *, int, int, int);
template int send(const unsigned long *, int, int, int);
template int send(const float *, int, int, int);
template int send(const double *, int, int, int); 
template int send(const long double *, int, int, int);

template int receive(char *, int, int, int);
template int receive(unsigned char *, int, int, int);
template int receive(short *, int, int, int);
template int receive(unsigned short *, int, int, int);
template int receive(int *, int, int, int);
template int receive(unsigned int *, int, int, int);
template int receive(long *, int, int, int);
template int receive(unsigned long *, int, int, int);
template int receive(float *, int, int, int);
template int receive(double *, int, int, int); 
template int receive(long double *, int, int, int);

template int broadcast(const char *, int, int);
template int broadcast(const unsigned char *, int, int);
template int broadcast(const short *, int, int);
template int broadcast(const unsigned short *, int, int);
template int broadcast(const int *, int, int);
template int broadcast(const unsigned int *, int, int);
template int broadcast(const long *, int, int);
template int broadcast(const unsigned long *, int, int);
template int broadcast(const float *, int, int);
template int broadcast(const double *, int, int); 
template int broadcast(const long double *, int, int);

} // end namespace rtt_c4

#endif // C4_MPI

//---------------------------------------------------------------------------//
//                              end of C4_MPI_blocking_pt.cc
//---------------------------------------------------------------------------//
