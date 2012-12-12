//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstC4_Req.cc
 * \author Kelly Thompson
 * \date   Tue Nov  1 15:49:44 2005
 * \brief  Unit test for C4_Req class.
 * \note   Copyright (C) 2006-2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "ds++/Assert.hh"
#include "../global.hh"
#include "../SpinLock.hh"
// #include "ds++/Release.hh"
#include "c4_test.hh"
#include "C4_Functions.hh"

#ifdef C4_MPI
#include "../MPI_Traits.hh"
#include <mpi.h>
#endif

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstCopyConstructor()
{
    using rtt_c4::C4_Req;
    
    C4_Req requestA;
    C4_Req requestB( requestA );

    // The behavior of the copy constructor is not obvious.  If requestA has
    // not been used (inuse() returns 0) then requestA != requestB.
    
    if( requestA.inuse() == 0 && requestA == requestB  )
    {
        FAILMSG("requestA.inuse()==0, so requestA cannot == requestB.");
    }

    if( requestA.inuse() == 0 && requestA != requestB  )
    {
        PASSMSG("requestA.inuse()==0 and requestA != requestB.");
    }

    if( requestA.inuse() == 1 && requestA == requestB  )
    {
        PASSMSG("requestA.inuse()=1 and requestA == requestB.");
    }

    if( requestA.inuse() == 1 && requestA != requestB  )
    {
        FAILMSG("requestA.inuse()=1, so requestA must == requestB.");
    }
    
    return;
}

//---------------------------------------------------------------------------//

void tstTraits()
{
    using std::cout;
    using std::endl;
    using rtt_c4::C4_Traits;

    {
        rtt_c4::HSyncSpinLock headsyncspinlock;
        
        if( C4_Traits<unsigned char>::tag  != 432 ) ITFAILS;
        if( C4_Traits<short>::tag          != 433 ) ITFAILS;
        if( C4_Traits<unsigned short>::tag != 434 ) ITFAILS;
        if( C4_Traits<unsigned int>::tag   != 436 ) ITFAILS;
        if( C4_Traits<unsigned long>::tag  != 438 ) ITFAILS;
        if( C4_Traits<long double>::tag    != 441 ) ITFAILS;
    }
#ifdef C4_MPI
    {
        using rtt_c4::MPI_Traits;
        rtt_c4::TSyncSpinLock tailsyncspinlock;
        
        if( MPI_Traits<unsigned char>::element_type()  != MPI_UNSIGNED_CHAR ) ITFAILS;
        if( MPI_Traits<short>::element_type()          != MPI_SHORT ) ITFAILS;
        if( MPI_Traits<unsigned short>::element_type() != MPI_UNSIGNED_SHORT ) ITFAILS;
        if( MPI_Traits<unsigned int>::element_type()   != MPI_UNSIGNED ) ITFAILS;
        if( MPI_Traits<unsigned long>::element_type()  != MPI_UNSIGNED_LONG ) ITFAILS;
        if( MPI_Traits<long double>::element_type()    != MPI_LONG_DOUBLE ) ITFAILS;
    }
#endif
    
    return;
}

//---------------------------------------------------------------------------------------//
void tstWait()
{
    using namespace rtt_c4;
    
    if (rtt_c4::node()>0)
    {
        cout << "sending from processor " << processor_name() << ':' << endl;
        int buffer[1];
        buffer[0] = node();
        C4_Req outgoing = send_async(buffer, 1U, 0);
        unsigned result = wait_any(1U, &outgoing);
        if (result!=0)
        {
            ITFAILS;
        }
    }
    else
    {
        cout << "receiving to processor " << processor_name() << ':' << endl;
        Check(rtt_c4::nodes()<5);
        C4_Req requests[4];
        bool done[4];
        for (int p=1; p<nodes(); ++p)
        {
            int buffer[4][1];
            requests[p] = receive_async(buffer[p], 1U, p);
            done[p] = false;
        }
        for (int c=1; c<nodes(); ++c)
        {
            unsigned result = wait_any(nodes(), requests);
            if (done[result])
            {
                ITFAILS;
            }
            done[result] = true;
        }
        for (int p=1; p<nodes(); ++p)
        {
            if (!done[p])
            {
                ITFAILS;
            }
        }
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;

    // Use deprecated form of  rtt_c4::initialize(argc, argv) to ensure that
    // it gets tested:
    C4::Init(argc,argv);

    try
    {
        // >>> UNIT TESTS
        tstCopyConstructor();
        tstTraits();
        tstWait();
    }
    catch (std::exception &err)
    {
        cout << "ERROR: While testing tstC4_Req, " 
                  << err.what() << endl;
        rtt_c4::abort();
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstC4_Req, " 
                  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << endl;
        rtt_c4::abort();
        return 1;
    }

    {
        rtt_c4::HTSyncSpinLock slock;

        // status of test
        cout <<   "\n*********************************************\n";
        if (rtt_c4_test::passed) 
            cout << "**** tstC4_Req Test: PASSED on " << rtt_c4::node();
        cout <<   "\n*********************************************\n" << endl;
    }
    
    // Use depreicated form of rtt_c4::global_barrier() to ensure that it gets
    // tested:
    C4::gsync();
    cout << "Done testing tstC4_Req on " << rtt_c4::node() << endl;
    // Use depreicated form of rtt_c4::finalize() to ensure that it gets
    // tested.
    C4::Finalize();
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstC4_Req.cc
//---------------------------------------------------------------------------//
