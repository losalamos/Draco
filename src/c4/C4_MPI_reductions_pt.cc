//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/C4_MPI_reductions_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 11:12:35 2002
 * \brief  C4 MPI global reduction instantiations.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "c4/config.h"

#ifdef C4_MPI
#include "C4_MPI.t.hh"
#else
#include "C4_Serial.t.hh"
#endif

namespace rtt_c4 {

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF GLOBAL REDUCTIONS
//----------------------------------------------------------------------------//

template DLL_PUBLIC_c4 void global_sum<short>(short &);
template DLL_PUBLIC_c4 void global_sum<unsigned short>(unsigned short &);
template DLL_PUBLIC_c4 void global_sum<int>(int &);
template DLL_PUBLIC_c4 void global_sum<unsigned int>(unsigned int &);
template DLL_PUBLIC_c4 void global_sum<long>(long &);
template DLL_PUBLIC_c4 void global_sum<unsigned long>(unsigned long &);
template DLL_PUBLIC_c4 void global_sum<float>(float &);
template DLL_PUBLIC_c4 void global_sum<double>(double &);
template DLL_PUBLIC_c4 void global_sum<long double>(long double &);
template DLL_PUBLIC_c4 void global_sum<long long>(long long &);
template DLL_PUBLIC_c4 void
global_sum<unsigned long long>(unsigned long long &);

template DLL_PUBLIC_c4 void global_isum<short>(short &, short &, C4_Req &);
template DLL_PUBLIC_c4 void
global_isum<unsigned short>(unsigned short &, unsigned short &, C4_Req &);
template DLL_PUBLIC_c4 void global_isum<int>(int &, int &, C4_Req &);
template DLL_PUBLIC_c4 void global_isum<unsigned int>(unsigned int &,
                                                      unsigned int &, C4_Req &);
template DLL_PUBLIC_c4 void global_isum<long>(long &, long &, C4_Req &);
template DLL_PUBLIC_c4 void
global_isum<unsigned long>(unsigned long &, unsigned long &, C4_Req &);
template DLL_PUBLIC_c4 void global_isum<float>(float &, float &, C4_Req &);
template DLL_PUBLIC_c4 void global_isum<double>(double &, double &, C4_Req &);
template DLL_PUBLIC_c4 void global_isum<long double>(long double &,
                                                     long double &, C4_Req &);
template DLL_PUBLIC_c4 void global_isum<long long>(long long &, long long &,
                                                   C4_Req &);
template DLL_PUBLIC_c4 void
global_isum<unsigned long long>(unsigned long long &, unsigned long long &,
                                C4_Req &);

template DLL_PUBLIC_c4 void global_prod<short>(short &);
template DLL_PUBLIC_c4 void global_prod<unsigned short>(unsigned short &);
template DLL_PUBLIC_c4 void global_prod<int>(int &);
template DLL_PUBLIC_c4 void global_prod<unsigned int>(unsigned int &);
template DLL_PUBLIC_c4 void global_prod<long>(long &);
template DLL_PUBLIC_c4 void global_prod<unsigned long>(unsigned long &);
template DLL_PUBLIC_c4 void global_prod<float>(float &);
template DLL_PUBLIC_c4 void global_prod<double>(double &);
template DLL_PUBLIC_c4 void global_prod<long double>(long double &);
template DLL_PUBLIC_c4 void global_prod<long long>(long long &);
template DLL_PUBLIC_c4 void
global_prod<unsigned long long>(unsigned long long &);

template DLL_PUBLIC_c4 void global_max<short>(short &);
template DLL_PUBLIC_c4 void global_max<unsigned short>(unsigned short &);
template DLL_PUBLIC_c4 void global_max<int>(int &);
template DLL_PUBLIC_c4 void global_max<unsigned int>(unsigned int &);
template DLL_PUBLIC_c4 void global_max<long>(long &);
template DLL_PUBLIC_c4 void global_max<unsigned long>(unsigned long &);
template DLL_PUBLIC_c4 void global_max<float>(float &);
template DLL_PUBLIC_c4 void global_max<double>(double &);
template DLL_PUBLIC_c4 void global_max<long double>(long double &);
template DLL_PUBLIC_c4 void global_max<long long>(long long &);
template DLL_PUBLIC_c4 void
global_max<unsigned long long>(unsigned long long &);

template DLL_PUBLIC_c4 void global_min<short>(short &);
template DLL_PUBLIC_c4 void global_min<unsigned short>(unsigned short &);
template DLL_PUBLIC_c4 void global_min<int>(int &);
template DLL_PUBLIC_c4 void global_min<unsigned int>(unsigned int &);
template DLL_PUBLIC_c4 void global_min<long>(long &);
template DLL_PUBLIC_c4 void global_min<unsigned long>(unsigned long &);
template DLL_PUBLIC_c4 void global_min<float>(float &);
template DLL_PUBLIC_c4 void global_min<double>(double &);
template DLL_PUBLIC_c4 void global_min<long double>(long double &);
template DLL_PUBLIC_c4 void
global_min<unsigned long long>(unsigned long long &);

template DLL_PUBLIC_c4 void global_sum<short>(short *, int);
template DLL_PUBLIC_c4 void global_sum<unsigned short>(unsigned short *, int);
template DLL_PUBLIC_c4 void global_sum<int>(int *, int);
template DLL_PUBLIC_c4 void global_sum<unsigned int>(unsigned int *, int);
template DLL_PUBLIC_c4 void global_sum<long>(long *, int);
template DLL_PUBLIC_c4 void global_sum<unsigned long>(unsigned long *, int);
template DLL_PUBLIC_c4 void global_sum<float>(float *, int);
template DLL_PUBLIC_c4 void global_sum<double>(double *, int);
template DLL_PUBLIC_c4 void global_sum<long double>(long double *, int);
template DLL_PUBLIC_c4 void global_sum<long long>(long long *, int);
template DLL_PUBLIC_c4 void global_sum<unsigned long long>(unsigned long long *,
                                                           int);

template DLL_PUBLIC_c4 void global_prod<short>(short *, int);
template DLL_PUBLIC_c4 void global_prod<unsigned short>(unsigned short *, int);
template DLL_PUBLIC_c4 void global_prod<int>(int *, int);
template DLL_PUBLIC_c4 void global_prod<unsigned int>(unsigned int *, int);
template DLL_PUBLIC_c4 void global_prod<long>(long *, int);
template DLL_PUBLIC_c4 void global_prod<unsigned long>(unsigned long *, int);
template DLL_PUBLIC_c4 void global_prod<float>(float *, int);
template DLL_PUBLIC_c4 void global_prod<double>(double *, int);
template DLL_PUBLIC_c4 void global_prod<long double>(long double *, int);
template DLL_PUBLIC_c4 void global_prod<long long>(long long *, int);
template DLL_PUBLIC_c4 void
global_prod<unsigned long long>(unsigned long long *, int);

template DLL_PUBLIC_c4 void global_max<short>(short *, int);
template DLL_PUBLIC_c4 void global_max<unsigned short>(unsigned short *, int);
template DLL_PUBLIC_c4 void global_max<int>(int *, int);
template DLL_PUBLIC_c4 void global_max<unsigned int>(unsigned int *, int);
template DLL_PUBLIC_c4 void global_max<long>(long *, int);
template DLL_PUBLIC_c4 void global_max<unsigned long>(unsigned long *, int);
template DLL_PUBLIC_c4 void global_max<float>(float *, int);
template DLL_PUBLIC_c4 void global_max<double>(double *, int);
template DLL_PUBLIC_c4 void global_max<long double>(long double *, int);
template DLL_PUBLIC_c4 void global_max<long long>(long long *, int);
template DLL_PUBLIC_c4 void global_max<unsigned long long>(unsigned long long *,
                                                           int);

template DLL_PUBLIC_c4 void global_min<short>(short *, int);
template DLL_PUBLIC_c4 void global_min<unsigned short>(unsigned short *, int);
template DLL_PUBLIC_c4 void global_min<int>(int *, int);
template DLL_PUBLIC_c4 void global_min<unsigned int>(unsigned int *, int);
template DLL_PUBLIC_c4 void global_min<long>(long *, int);
template DLL_PUBLIC_c4 void global_min<unsigned long>(unsigned long *, int);
template DLL_PUBLIC_c4 void global_min<float>(float *, int);
template DLL_PUBLIC_c4 void global_min<double>(double *, int);
template DLL_PUBLIC_c4 void global_min<long double>(long double *, int);
template DLL_PUBLIC_c4 void global_min<long long>(long long *, int);
template DLL_PUBLIC_c4 void global_min<unsigned long long>(unsigned long long *,
                                                           int);

} // end namespace rtt_c4

//----------------------------------------------------------------------------//
// end of C4_MPI_reductions_pt.cc
//----------------------------------------------------------------------------//
