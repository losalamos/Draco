//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   traits/test/Traits_Test.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:49:41 2000
 * \brief  Traits testing definitions.
 * \note   Copyright (C) 2000-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __traits_test_Traits_Test_hh__
#define __traits_test_Traits_Test_hh__

#include <string>
#include <iostream>

namespace rtt_traits_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

inline bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

} // end namespace rtt_traits_test

#endif // __traits_test_Traits_Test_hh__

//---------------------------------------------------------------------------//
// end of traits/test/Traits_Test.hh
//---------------------------------------------------------------------------//
