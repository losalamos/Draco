//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/Viz_Test.hh
 * \author Thomas M. Evans
 * \date   Mon Jan 24 11:13:39 2000
 * \brief  Viz_Test common header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __viz_test_Viz_Test_hh__
#define __viz_test_Viz_Test_hh__

#include <string>
#include <iostream>

namespace rtt_viz_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

inline bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

} // end namespace rtt_viz_test

#endif                          // __viz_test_Viz_Test_hh__

//---------------------------------------------------------------------------//
//                              end of viz/test/Viz_Test.hh
//---------------------------------------------------------------------------//
