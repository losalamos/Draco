//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/test/tstTiming.cc
 * \author Thomas M. Evans
 * \date   Mon Dec 12 15:32:10 2005
 * \brief  
 * \note   Copyright (C) 2004-2010 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//


#include "diagnostics_test.hh"
#include "ds++/Release.hh"
#include "../Timing.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;
using rtt_dsxx::soft_equiv;
typedef rtt_diagnostics::Timing_Diagnostics D;

int nodes = 0;
int node  = 0;

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//

void do_A()
{
    TIMER(A_timer);
    TIMER_START(A_timer);

    // do a mat-vec multiply

    int S = 9000;
    
    vector<double> b(S, 0.0);
    vector<double> x(S, 0.0);

    double A = 1.0;
    double B = 1.1;
    double C = 1.01;

    x[0] = 0.2;
    for (int i = 1; i < S; i++)
        x[i] = i + 1.1 + x[i-1];

    for (int i = 1; i < (S - 1); ++i)
    {
        b[i] = x[i-1] * A + x[i] * B + x[i+1] * C;
    }

    b[0]   = B * x[0] + C * x[1];
    b[S-1] = A * x[S-2] + B * x[S-1]; 

    TIMER_STOP(A_timer);
    TIMER_RECORD("A_iteration", A_timer);
}

//---------------------------------------------------------------------------//

void do_B()
{
    TIMER(B_timer);
    TIMER_START(B_timer);

    // do a mat-vec multiply

    int S = 6000;
    
    vector<double> b(S, 0.0);
    vector<double> x(S, 0.0);

    double A = 1.0;
    double B = 1.1;
    double C = 1.01;

    x[0] = 0.2;
    for (int i = 1; i < S; i++)
        x[i] = i + 1.1 + x[i-1];

    for (int i = 1; i < (S - 1); ++i)
    {
        b[i] = x[i-1] * A + x[i] * B + x[i+1] * C;
    }

    b[0]   = B * x[0] + C * x[1];
    b[S-1] = A * x[S-2] + B * x[S-1]; 

    TIMER_STOP(B_timer);
    TIMER_RECORD("B_iteration", B_timer);
}

//---------------------------------------------------------------------------//

void do_C()
{
    TIMER(C_timer);
    TIMER_START(C_timer);

    // do a mat-vec multiply

    int S = 3000;
    
    vector<double> b(S, 0.0);
    vector<double> x(S, 0.0);

    double A = 1.0;
    double B = 1.1;
    double C = 1.01;

    x[0] = 0.2;
    for (int i = 1; i < S; i++)
        x[i] = i + 1.1 + x[i-1];

    for (int i = 1; i < (S - 1); ++i)
    {
        b[i] = x[i-1] * A + x[i] * B + x[i+1] * C;
    }

    b[0]   = B * x[0] + C * x[1];
    b[S-1] = A * x[S-2] + B * x[S-1]; 

    TIMER_STOP(C_timer);
    TIMER_RECORD("C_iteration", C_timer);
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void timing_active()
{
#ifdef DRACO_TIMING_ON
    cout << ">>> Testing timing macros with value " << DRACO_TIMING << endl;
#else
    cout << ">>> Timing macros inactive" << endl;
#endif
}

//---------------------------------------------------------------------------//

void test_timing()
{
    // add to some timers
    D::update_timer("A", 1.2);
    D::update_timer("B", 1.1);
    D::update_timer("B", 2.3);

    if (!soft_equiv(D::timer_value("A"), 1.2)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 3.4)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;

    D::reset_timer("B");
    D::update_timer("A", 1.3);
    if (!soft_equiv(D::timer_value("A"), 2.5)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;

    vector<string> timers = D::timer_keys();
    if (timers.size() != 3) ITFAILS;
    if (timers[0] != "A") ITFAILS;
    if (timers[1] != "B") ITFAILS;
    if (timers[2] != "C") ITFAILS;

    D::delete_timer("B");
    timers = D::timer_keys();
    if (timers.size() != 2) ITFAILS;
    if (timers[0] != "A") ITFAILS;
    if (timers[1] != "C") ITFAILS;

    // calling timer_value on B will get it back
    if (!soft_equiv(D::timer_value("A"), 2.5)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;
    timers = D::timer_keys();
    if (timers.size() != 3)   ITFAILS;
    if (D::num_timers() != 3) ITFAILS;

    // delete all timers
    D::delete_timers();
    if (D::num_timers() != 0) ITFAILS;
    timers = D::timer_keys();
    if (timers.size() != 0)   ITFAILS;

    D::update_timer("B", 12.4);
    D::update_timer("C", 1.3);
    if (!soft_equiv(D::timer_value("A"), 0.0))  ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 12.4)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 1.3))  ITFAILS;

    // reset all timers
    D::reset_timers();
    if (!soft_equiv(D::timer_value("A"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;
    
    if (rtt_diagnostics_test::passed)
    {
        PASSMSG("Diagnostics timer lists ok.");
    }
}

//---------------------------------------------------------------------------//

void test_macros()
{
    // delete all existing timers
    D::delete_timers();

    // make timers and do results
    TIMER(outer_timer);
    TIMER_START(outer_timer);

    do_A();
    do_B();
    do_C();

    TIMER_STOP(outer_timer);
    TIMER_RECORD("Outer", outer_timer);

    // if the timers are off we get no timing data
    vector<string> keys = D::timer_keys();
    if (DRACO_TIMING == 0)
    {
        if (keys.size() != 0)     ITFAILS;
        if (D::num_timers() != 0) ITFAILS;
    }
    else
    {
        if (keys.size() != 4) ITFAILS;
        cout << setw(15) << "Routine" << setw(15) << "Fraction" << endl;
        cout << "------------------------------" << endl;

        // get the keys and print a table
        double total = D::timer_value("Outer");
        if (total <= 0.0) ITFAILS;

        cout.precision(4);
        cout.setf(ios::fixed, ios::floatfield);
        
        for (int i = 0, N = keys.size(); i < N; ++i)
        {
            double fraction = D::timer_value(keys[i]) / total;
            cout << setw(15) << keys[i] << setw(15) << fraction << endl;
        }

        cout << "The total time was " << total << endl;
        cout << endl;
    }

    TIMER_REPORT(outer_timer, cout, "Total time");

    if (rtt_diagnostics_test::passed)
    {
        PASSMSG("Timer macros ok.");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    nodes = rtt_c4::nodes();
    node  = rtt_c4::node();

    if (node > 0)
    {
        rtt_c4::finalize();
        return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
        {
            if (rtt_c4::node() == 0)
                cout << argv[0] << ": version " 
                     << rtt_dsxx::release() 
                     << endl;
            rtt_c4::finalize();
            return 0;
        }

    try
    {
        // >>> UNIT TESTS
        timing_active();
        test_timing();
        test_macros();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstTiming, " 
                  << err.what()
                  << std::endl;
        rtt_c4::abort();
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstTiming, " 
                  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
        rtt_c4::abort();
        return 1;
    }

    {
        // status of test
        std::cout << std::endl;
        std::cout <<     "*********************************************" 
                  << std::endl;
        if (rtt_diagnostics_test::passed) 
        {
            std::cout << "**** tstTiming Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
        }
        std::cout <<     "*********************************************" 
                  << std::endl;
        std::cout << std::endl;
    }

    std::cout << "Done testing tstTiming on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstTiming.cc
//---------------------------------------------------------------------------//
