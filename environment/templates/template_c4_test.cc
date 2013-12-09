//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright (C) 2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"

using namespace std;

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
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        // >>> UNIT TESTS
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of <class>.cc
//---------------------------------------------------------------------------//
