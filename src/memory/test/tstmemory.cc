//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   memory/test/tstmemory.cc
 * \author Kent G. Budge
 * \brief  memory test.
 * \note   Copyright (C) 2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
//! \version $Id: tstProcmon.cc 5830 2011-05-05 19:43:43Z kellyt $
//---------------------------------------------------------------------------//

#include "../memory.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <sstream>

using namespace std;
using namespace rtt_memory;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tst_memory( rtt_dsxx::UnitTest & ut )
{
    if (total_allocation()==0)
    {
        ut.passes("correct initial total allocation");
    }
    else
    {
        ut.failure("NOT correct initial total allocation");
    }
    if (peak_allocation()==0)
    {
        ut.passes("correct initial peak allocation");
    }
    else
    {
        ut.failure("NOT correct initial peak allocation");
    }
    set_memory_checking(true);
    
    double *array = new double[20];

    double *array2 = new double[30];

#if DRACO_DIAGNOSTICS & 2
    if (total_allocation()==50*sizeof(double))
    {
        ut.passes("correct total allocation");
    }
    else
    {
        ut.failure("NOT correct total allocation");
    }
    if (peak_allocation()>=50*sizeof(double))
    {
        ut.passes("correct peak allocation");
    }
    else
    {
        ut.failure("NOT correct peak allocation");
    }
#else
    ut.passes("memory diagnostics not checked for this build");
#endif

    delete[] array;
    delete[] array2;

#if DRACO_DIAGNOSTICS & 2
    if (total_allocation()==0)
    {
        ut.passes("correct total allocation");
    }
    else
    {
        ut.failure("NOT correct total allocation");
    }
    if (peak_allocation()>=50*sizeof(double))
    {
        ut.passes("correct peak allocation");
    }
    else
    {
        ut.failure("NOT correct peak allocation");
    }
#endif

}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut( argc, argv, rtt_dsxx::release );
    try
    {
        // >>> UNIT TESTS
        tst_memory(ut);
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstProcmon.cc
//---------------------------------------------------------------------------//
