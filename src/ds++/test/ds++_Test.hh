//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/ds++_Test.hh
 * \author Thomas M. Evans
 * \date   Fri Jul 20 17:25:07 2001
 * \brief  Header for testing functionality in ds++.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ds_test_dsxx_Test_hh__
#define __ds_test_dsxx_Test_hh__

#include <iostream>

namespace rtt_dsxx_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

inline bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

} // end namespace rtt_dsxx_test

#endif                          // __ds_test_dsxx_Test_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/test/ds++_Test.hh
//---------------------------------------------------------------------------//
