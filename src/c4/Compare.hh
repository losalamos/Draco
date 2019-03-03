//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Compare.hh
 * \author Mike Buksas
 * \brief  Define class Compare
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef c4_Compare_hh
#define c4_Compare_hh

#include "ds++/config.h"

namespace rtt_c4 {

//===========================================================================//
/*!
 * \class Compare
 *
 * \example c4/test/tstCompare.cc
 */
//===========================================================================//

DLL_PUBLIC_c4 bool check_global_equiv(int local_value);
DLL_PUBLIC_c4 bool check_global_equiv(unsigned long long local_value);
DLL_PUBLIC_c4 bool check_global_equiv(unsigned long local_value);
DLL_PUBLIC_c4 bool check_global_equiv(long long local_value);
DLL_PUBLIC_c4 bool check_global_equiv(long local_value);
DLL_PUBLIC_c4 bool check_global_equiv(double local_value, double eps = 1.0e-8);

} // end namespace rtt_c4

#endif // c4_Compare_hh

//---------------------------------------------------------------------------//
// end of c4/Compare.hh
//---------------------------------------------------------------------------//
