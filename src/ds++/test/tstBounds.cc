//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstBounds.cc.cc
 * \author Kelly Thompson
 * \date   Tue Aug  2 11:00:19 2005
 * \brief  Test the Bounds Class (Bounds.hh)
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "../Assert.hh"
#include "../Release.hh"
#include "../Bounds.hh"
#include "ds_test.hh"

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (string(argv[arg]) == "--version")
        {
            cout << argv[0] << ": version " 
                      << rtt_dsxx::release() 
                      << endl;
            return 0;
        }

    try
    {
        // >>> UNIT TESTS
        int const bmin( -3 );
        int const bmax(  7 );
        rtt_dsxx::Bounds b( bmin, bmax );

        if( b.min() == bmin )
            PASSMSG("Found correct min value.");
        else
            FAILMSG("Found incorrect min value.");

        if( b.max() == bmax )
            PASSMSG("Found correct max value.");
        else
            FAILMSG("Found incorrect max value.");
        
        if( b.len() == bmax - bmin + 1 )
            PASSMSG("Found correct length value.");
        else
            FAILMSG("Found incorrect length value.");
        
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstBounds.cc, " 
             << err.what()
             << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstBounds.cc, " 
             << "An unknown exception was thrown."
             << endl;
        return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" 
              << endl;
    if (rtt_ds_test::passed) 
    {
        cout << "**** tstBounds.cc Test: PASSED" 
                  << endl;
    }
    cout <<     "*********************************************" 
              << endl;
    cout << endl;
    
    cout << "Done testing tstBounds.cc." << endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstBounds.cc.cc
//---------------------------------------------------------------------------//
