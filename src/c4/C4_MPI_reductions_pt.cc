//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI_reductions_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 11:12:35 2002
 * \brief  C4 MPI global reduction instantiations.
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
// EXPLICIT INSTANTIATIONS OF GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template void global_sum(short &);
template void global_sum(unsigned short &);
template void global_sum(int &);
template void global_sum(unsigned int &);
template void global_sum(long &);
template void global_sum(unsigned long &);
template void global_sum(float &);
template void global_sum(double &);
template void global_sum(long double &);

template void global_prod(short &);
template void global_prod(unsigned short &);
template void global_prod(int &);
template void global_prod(unsigned int &);
template void global_prod(long &);
template void global_prod(unsigned long &);
template void global_prod(float &);
template void global_prod(double &);
template void global_prod(long double &);

template void global_max(short &);
template void global_max(unsigned short &);
template void global_max(int &);
template void global_max(unsigned int &);
template void global_max(long &);
template void global_max(unsigned long &);
template void global_max(float &);
template void global_max(double &);
template void global_max(long double &);

template void global_min(short &);
template void global_min(unsigned short &);
template void global_min(int &);
template void global_min(unsigned int &);
template void global_min(long &);
template void global_min(unsigned long &);
template void global_min(float &);
template void global_min(double &);
template void global_min(long double &);

template void global_sum(short *, int);
template void global_sum(unsigned short *, int);
template void global_sum(int *, int);
template void global_sum(unsigned int *, int);
template void global_sum(long *, int);
template void global_sum(unsigned long *, int);
template void global_sum(float *, int);
template void global_sum(double *, int);
template void global_sum(long double *, int);

template void global_prod(short *, int);
template void global_prod(unsigned short *, int);
template void global_prod(int *, int);
template void global_prod(unsigned int *, int);
template void global_prod(long *, int);
template void global_prod(unsigned long *, int);
template void global_prod(float *, int);
template void global_prod(double *, int);
template void global_prod(long double *, int);

template void global_max(short *, int);
template void global_max(unsigned short *, int);
template void global_max(int *, int);
template void global_max(unsigned int *, int);
template void global_max(long *, int);
template void global_max(unsigned long *, int);
template void global_max(float *, int);
template void global_max(double *, int);
template void global_max(long double *, int);

template void global_min(short *, int);
template void global_min(unsigned short *, int);
template void global_min(int *, int);
template void global_min(unsigned int *, int);
template void global_min(long *, int);
template void global_min(unsigned long *, int);
template void global_min(float *, int);
template void global_min(double *, int);
template void global_min(long double *, int);

} // end namespace rtt_c4

#endif // C4_MPI

//---------------------------------------------------------------------------//
//                              end of C4_MPI_reductions_pt.cc
//---------------------------------------------------------------------------//
