//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstScalar.cc
 * \author Kelly Thompson
 * \date   Tue Nov  1 13:24:19 2005
 * \brief  Test functions provided by C4_Serial.cc/.hh
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "ds++/Assert.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "ds++/Release.hh"
#include "../C4_Tags.hh"
#include "c4_test.hh"

using namespace std;

//---------------------------------------------------------------------------//
// TESTS:
//
// These tests are designed to test functionality provided in
// C4_Serial.cc/.hh. 
//---------------------------------------------------------------------------//

void tstScalar()
{

// Skip the tests if code not configured with the option --with-c4=scalar.
    
#ifndef C4_SCALAR

    if( rtt_c4::isScalar() )
        FAILMSG("Incorrectly identified process as scalar.")
    else
        PASSMSG("Correctly identified process as parallel.")

#else

    // Check the isScalar function.
    if( rtt_c4::isScalar() )
        PASSMSG("Correctly identified process as scalar.")
    else
        FAILMSG("Incorrectly identified process as parallel.")
    
// For --with-c4=scalar, probe(int,int,int) always returns false.

    int int3(99);
    bool const probeResult( rtt_c4::probe(0,0,int3) );
    if( probeResult )
    {
        FAILMSG( "For --with-c4=scalar, probe(int,int,int) returned true." );
    }
    else
    {
        PASSMSG( "For --with-c4=scalar, probe(int,int,int) returned false." );
    }

    // Test broadcast function

    int result(0);
    result = rtt_c4::broadcast( &int3, 1, 0 );
    if( result == rtt_c4::C4_SUCCESS )
    {
        PASSMSG("For --with-c4=scalar, broadcast() returned C4_SUCCCESS.");
    }
    else
    {
        PASSMSG("For --with-c4=scalar, broadcast() did not return C4_SUCCCESS.");
    }

    // Test gather and scatter
    
    result = rtt_c4::gather(&int3, &int3, 1);
    if( result == rtt_c4::C4_SUCCESS )
    {
        PASSMSG("For --with-c4=scalar, gather() returned C4_SUCCCESS.");
    }
    else
    {
        PASSMSG("For --with-c4=scalar, gather() did not return C4_SUCCCESS.");
    }
    
    result = rtt_c4::scatter(&int3, &int3, 1);
    if( result == rtt_c4::C4_SUCCESS )
    {
        PASSMSG("For --with-c4=scalar, scatter() returned C4_SUCCCESS.");
    }
    else
    {
        PASSMSG("For --with-c4=scalar, scatter() did not return C4_SUCCCESS.");
    }

    // Test send function
    result = rtt_c4::send( &int3, 1, 0, 0 );
    if( result == rtt_c4::C4_SUCCESS )
    {
        PASSMSG("For --with-c4=scalar, send() returned C4_SUCCCESS.");
    }
    else
    {
        PASSMSG("For --with-c4=scalar, send() did not return C4_SUCCCESS.");
    }
    
    // Test receive function
    result = rtt_c4::receive( &int3, 1, 0, 0 );
    if( result == rtt_c4::C4_SUCCESS )
    {
        PASSMSG("For --with-c4=scalar, receive() returned C4_SUCCCESS.");
    }
    else
    {
        PASSMSG("For --with-c4=scalar, receive() did not return C4_SUCCCESS.");
    }

    // Test send_async<T>(const T*,int,int,int) and
    // receive_async<T>(T*,int,int,int) functions:
    rtt_c4::C4_Req request( rtt_c4::send_async( &int3, 1, 0, 0 ) );
    rtt_c4::receive_async( request, &int3, 1, 0, 0 );

    request = rtt_c4::receive_async( &int3, 1, 0, 0 );
    rtt_c4::send_async( request, &int3, 1, 0, 0 );
    
    if( request.inuse() == 0 && request.count() == 0 )
    {
        PASSMSG("For --with-c4=scalar, successful test of send_async(const T*,int,int,int) and receive_async( C4_Req&,T*,int,int,int).");
    }
    else
    {
        FAILMSG("For --with-c4=scalar, unsuccessful test of send_async(const T*,int,int,int) and receive_async( C4_Req&,T*,int,int,int).");
    }

#endif

    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;
    std::string unitTestName( "tstScalar" );

    rtt_c4::initialize(argc, argv);

    try
    {
        // >>> UNIT TESTS
        tstScalar();
    }
    catch (std::exception &err)
    {
        cout << "ERROR: While testing " << unitTestName << ", "
             << err.what() << endl;
        rtt_c4::abort();
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << unitTestName << ", " 
             << "An unknown exception was thrown on processor "
             << rtt_c4::node() << endl;
        rtt_c4::abort();
        return 1;
    }

    {
        rtt_c4::HTSyncSpinLock slock;
        
        // status of test
        cout << "\n*********************************************\n";
        if (rtt_c4_test::passed) 
            cout << "**** " << unitTestName
                 << " Test: PASSED on " << rtt_c4::node();
        cout << "\n*********************************************\n"
             << endl;
    }
    
    rtt_c4::global_barrier();
    cout << "Done testing " << unitTestName << " on "
         << rtt_c4::node() << endl; 
    rtt_c4::finalize();
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstScalar.cc.cc
//---------------------------------------------------------------------------//
