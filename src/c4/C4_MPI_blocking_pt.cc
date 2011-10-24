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

template int send<char>(const char *, int, int, int);
template int send<unsigned char>(const unsigned char *, int, int, int);
template int send<short>(const short *, int, int, int);
template int send<unsigned short>(const unsigned short *, int, int, int);
template int send<int>(const int *, int, int, int);
template int send<unsigned int>(const unsigned int *, int, int, int);
template int send<long>(const long *, int, int, int);
template int send<unsigned long>(const unsigned long *, int, int, int);
template int send<unsigned long long>(const unsigned long long *, int, int, int);
template int send<float>(const float *, int, int, int);
template int send<double>(const double *, int, int, int); 
template int send<long double>(const long double *, int, int, int);
    

template int receive<char>(char *, int, int, int);
template int receive<unsigned char>(unsigned char *, int, int, int);
template int receive<short>(short *, int, int, int);
template int receive<unsigned short>(unsigned short *, int, int, int);
template int receive<int>(int *, int, int, int);
template int receive<unsigned int>(unsigned int *, int, int, int);
template int receive<long>(long *, int, int, int);
template int receive<unsigned long>(unsigned long *, int, int, int);
template int receive<unsigned long long>(unsigned long long *, int, int, int);
template int receive<float>(float *, int, int, int);
template int receive<double>(double *, int, int, int); 
template int receive<long double>(long double *, int, int, int);

template int broadcast<char>(char *, int, int);
template int broadcast<unsigned char>(unsigned char *, int, int);
template int broadcast<short>(short *, int, int);
template int broadcast<unsigned short>(unsigned short *, int, int);
template int broadcast<int>(int *, int, int);
template int broadcast<unsigned int>(unsigned int *, int, int);
template int broadcast<long>(long *, int, int);
template int broadcast<unsigned long>(unsigned long *, int, int);
template int broadcast<unsigned long long>(unsigned long long *, int, int);
template int broadcast<float>(float *, int, int);
template int broadcast<double>(double *, int, int); 
template int broadcast<long double>(long double *, int, int);

} // end namespace rtt_c4

#endif // C4_MPI

//---------------------------------------------------------------------------//
//                              end of C4_MPI_blocking_pt.cc
//---------------------------------------------------------------------------//
