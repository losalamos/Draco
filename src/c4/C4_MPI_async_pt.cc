//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI_async_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:44:54 2002
 * \brief  C4 MPI non-blocking send/recv instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>

#ifdef C4_MPI

#include "C4_MPI.t.hh"

#else

#include "C4_Serial.hh"

#endif

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF NON-BLOCKING SEND/RECEIVE
//---------------------------------------------------------------------------//

template C4_Req send_async(const char *, int, int, int);
template C4_Req send_async(const unsigned char *, int, int, int);
template C4_Req send_async(const short *, int, int, int);
template C4_Req send_async(const unsigned short *, int, int, int);
template C4_Req send_async(const int *, int, int, int);
template C4_Req send_async(const unsigned int *, int, int, int);
template C4_Req send_async(const long *, int, int, int);
template C4_Req send_async(const unsigned long *, int, int, int);
template C4_Req send_async(const float *, int, int, int);
template C4_Req send_async(const double *, int, int, int);
template C4_Req send_async(const long double *, int, int, int);

template C4_Req receive_async(char *, int, int, int);
template C4_Req receive_async(unsigned char *, int, int, int);
template C4_Req receive_async(short *, int, int, int);
template C4_Req receive_async(unsigned short *, int, int, int);
template C4_Req receive_async(int *, int, int, int);
template C4_Req receive_async(unsigned int *, int, int, int);
template C4_Req receive_async(long *, int, int, int);
template C4_Req receive_async(unsigned long *, int, int, int);
template C4_Req receive_async(float *, int, int, int);
template C4_Req receive_async(double *, int, int, int);
template C4_Req receive_async(long double *, int, int, int);

template void send_async(C4_Req &, const char *, int, int, int);
template void send_async(C4_Req &, const unsigned char *, int, int, int);
template void send_async(C4_Req &, const short *, int, int, int);
template void send_async(C4_Req &, const unsigned short *, int, int, int);
template void send_async(C4_Req &, const int *, int, int, int);
template void send_async(C4_Req &, const unsigned int *, int, int, int);
template void send_async(C4_Req &, const long *, int, int, int);
template void send_async(C4_Req &, const unsigned long *, int, int, int);
template void send_async(C4_Req &, const float *, int, int, int);
template void send_async(C4_Req &, const double *, int, int, int);
template void send_async(C4_Req &, const long double *, int, int, int);

template void receive_async(C4_Req &, char *, int, int, int);
template void receive_async(C4_Req &, unsigned char *, int, int, int);
template void receive_async(C4_Req &, short *, int, int, int);
template void receive_async(C4_Req &, unsigned short *, int, int, int);
template void receive_async(C4_Req &, int *, int, int, int);
template void receive_async(C4_Req &, unsigned int *, int, int, int);
template void receive_async(C4_Req &, long *, int, int, int);
template void receive_async(C4_Req &, unsigned long *, int, int, int);
template void receive_async(C4_Req &, float *, int, int, int);
template void receive_async(C4_Req &, double *, int, int, int);
template void receive_async(C4_Req &, long double *, int, int, int);

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                              end of C4_MPI_async_pt.cc
//---------------------------------------------------------------------------//
