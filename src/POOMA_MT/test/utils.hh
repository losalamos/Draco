//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/test/utils.hh
 * \author Randy M. Roberts
 * \date   Wed Oct 27 14:32:09 1999
 * \brief  utilities for testing the POOMA_MT
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __POOMA_MT_test_utils_hh__
#define __POOMA_MT_test_utils_hh__

#include "PoomaMesh_XYZFactory.hh"
#include <string>

namespace rtt_POOMA_MT_test
{

void version(const std::string &progname);

PoomaMesh_XYZFactory getMTFactory(const std::string &filename);

} // end namespace rtt_POOMA_MT_test

#endif                          // __POOMA_MT_test_utils_hh__

//---------------------------------------------------------------------------//
//                              end of POOMA_MT/test/utils.hh
//---------------------------------------------------------------------------//
