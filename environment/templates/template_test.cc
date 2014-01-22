//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright (C) 2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace <namespace>;

#define PASSMSG(a) ut.passes(a)
#define ITFAILS    ut.failure(__LINE__)
#define FAILURE    ut.failure(__LINE__, __FILE__)
#define FAILMSG(a) ut.failure(a)

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//



//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        // >>> UNIT TESTS
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of <class>.cc
//---------------------------------------------------------------------------//
