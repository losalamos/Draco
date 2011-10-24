//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Compare.hh
 * \author Mike Buksas
 * \brief  Define class Compare
 * \note   Copyright (C) 2007 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_Compare_hh
#define c4_Compare_hh

namespace rtt_c4
{

//===========================================================================//
/*!
 * \class Compare
 * \brief
 *
 * \endcode
 */
/*! 
 * \example c4/test/tstCompare.cc
 *
 * Test of Compare.
 */
//===========================================================================//

bool check_global_equiv(int local_value);
bool check_global_equiv(unsigned long long local_value);
bool check_global_equiv(double local_value, double eps = 1.0e-8);

} // end namespace rtt_c4

#endif // c4_Compare_hh

//---------------------------------------------------------------------------//
//              end of c4/Compare.hh
//---------------------------------------------------------------------------//
