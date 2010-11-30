//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/test/tstDiagnostics.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  9 16:16:27 2005
 * \brief  Diagnostics test.
 * \note   Copyright (C) 2004-2010 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
//! \version $Id$
//---------------------------------------------------------------------------//

#include "diagnostics_test.hh"
#include "../Release.hh"
#include "../Diagnostics.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using namespace rtt_diagnostics;
using rtt_dsxx::soft_equiv;

int nodes = 0;
int node  = 0;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_ints()
{
    // add an integer quantity
    Diagnostics::integers["A"];
    if (Diagnostics::integers["A"] != 0) ITFAILS;

    // now update some values
    Diagnostics::integers["A"] = 21;
    Diagnostics::integers["A"]++;

    if (Diagnostics::integers["A"] != 22) ITFAILS;

    // add another
    Diagnostics::integers["B"] = 51;
    
    if (rtt_diagnostics_test::passed)
    {
        PASSMSG("Diagnostics integers ok.");
    }
}

//---------------------------------------------------------------------------//

void test_floats()
{
    if (Diagnostics::integers["A"] != 22) ITFAILS;
    if (Diagnostics::integers["B"] != 51) ITFAILS;

    // make a fraction entry
    Diagnostics::doubles["A_of_B"] = Diagnostics::integers["A"] /
                                     static_cast<double>(
                                         Diagnostics::integers["B"]);
    
    // check it
    if (!soft_equiv(Diagnostics::doubles["A_of_B"], 22.0/51.0)) ITFAILS;

    // erase and check
    if (Diagnostics::doubles.erase("A_of_B") != 1) ITFAILS;
    if (Diagnostics::doubles.count("A_of_B") != 0) ITFAILS;

    Diagnostics::integers.erase("A");
    
    if (rtt_diagnostics_test::passed)
    {
        PASSMSG("Diagnostics doubles ok.");
    }
}

//---------------------------------------------------------------------------//

void test_vectors()
{
    Diagnostics::vec_integers["A"];
    Diagnostics::vec_integers["B"];

    if (!Diagnostics::vec_integers["A"].empty()) ITFAILS;
    if (!Diagnostics::vec_integers["B"].empty()) ITFAILS;
    
    Diagnostics::vec_integers["A"].resize(10);
    if (Diagnostics::vec_integers["A"].size() != 10) ITFAILS;
    
    Diagnostics::vec_doubles["B"].resize(2);
    Diagnostics::vec_doubles["B"][0] = 1.1;
    Diagnostics::vec_doubles["B"][1] = 2.4;

    if (Diagnostics::vec_doubles["B"].size() != 2) ITFAILS;

    vector<double> ref(2, 1.1);
    ref[1] = 2.4;
    if (!soft_equiv(Diagnostics::vec_doubles["B"].begin(),
                    Diagnostics::vec_doubles["B"].end(),
                    ref.begin(), ref.end()))     ITFAILS;

    Diagnostics::vec_doubles["B"].clear();
    if (!Diagnostics::vec_integers["B"].empty()) ITFAILS;
    
    if (Diagnostics::integers["A"] != 0)  ITFAILS;
    if (Diagnostics::integers["B"] != 51) ITFAILS;
    
    if (rtt_diagnostics_test::passed)
    {
        PASSMSG("Diagnostics vectors ok.");
    }
}

//---------------------------------------------------------------------------//

void test_macro()
{
    cout << endl;
    
    int level[4];
    level[0] = 1;
    level[1] = 0;
    level[2] = 0;
    level[3] = 0;
    
#ifdef DRACO_DIAGNOSTICS_LEVEL_1
    cout << ">>> Testing Level 1 Block diagnostics." << endl;
    level[1] = 1;
    level[0] = 0;
#endif

#ifdef DRACO_DIAGNOSTICS_LEVEL_2
    cout << ">>> Testing Level 2 Block diagnostics." << endl;
    level[2] = 1;
    level[0] = 0;
#endif
    
#ifdef DRACO_DIAGNOSTICS_LEVEL_3
    cout << ">>> Testing Level 3 Block diagnostics." << endl;
    level[3] = 1;
    level[0] = 0;
#endif

    cout << endl;

    DIAGNOSTICS_ONE(integers["L1"] = 1);
    DIAGNOSTICS_TWO(integers["L2"] = 1);
    DIAGNOSTICS_THREE(integers["L3"] = 1);

    if (level[0] == 1)
    {
        if (Diagnostics::integers.count("L1") != 0) ITFAILS;
        if (Diagnostics::integers.count("L2") != 0) ITFAILS;
        if (Diagnostics::integers.count("L3") != 0) ITFAILS;
    }

    if (level[1] == 1)
    {
        if (Diagnostics::integers["L1"] != 1) ITFAILS;
    }

    if (level[2] == 1)
    {
        if (Diagnostics::integers["L2"] != 1) ITFAILS;
    }

    if (level[3] == 1)
    {
        if (Diagnostics::integers["L3"] != 1) ITFAILS;
    }
    
    if (rtt_diagnostics_test::passed)
    {
        PASSMSG("Diagnostics macro ok.");
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

    // print the version tag
    if (rtt_c4::node() == 0)
        cout << argv[0] << ": version " 
             << rtt_diagnostics::release() 
             << endl;
    
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
        {
            // Version tag has already printed.  If --version is found on the
            // command line then finalize and exit.
            rtt_c4::finalize();
            return 0;
        }

    try
    {
        // >>> UNIT TESTS

        test_ints();
        test_floats();
        test_vectors();

        test_macro();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstDiagnostics, " 
                  << err.what()
                  << std::endl;
        rtt_c4::abort();
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstDiagnostics, " 
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
            std::cout << "**** tstDiagnostics Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
        }
        std::cout <<     "*********************************************" 
                  << std::endl;
        std::cout << std::endl;
    }

    std::cout << "Done testing tstDiagnostics on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstDiagnostics.cc
//---------------------------------------------------------------------------//
