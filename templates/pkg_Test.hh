//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/<spkg>.hh
 * \author <user>
 * \date   <date>
 * \brief  <start>
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __<spkg>_hh__
#define __<spkg>_hh__

#include <iostream>

namespace rtt_<spkg>
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_<spkg>::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_<spkg>::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//
// USAGE:
// if (!condition) ITFAILS;

extern bool passed;

} // end namespace rtt_<spkg>

#define ITFAILS rtt_<spkg>::fail(__LINE__);
#define FAILURE rtt_<spkg>::fail(__LINE__, __FILE__);

#endif                          // __<spkg>_hh__

//---------------------------------------------------------------------------//
//                              end of <pkg>/<spkg>.hh
//---------------------------------------------------------------------------//
