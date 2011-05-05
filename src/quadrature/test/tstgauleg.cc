//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tstgauleg.cc
 * \author Kent Budge
 * \date   Tue Sep 27 12:49:39 2005
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"
#include "quadrature_test.hh"
#include "square.hh"
#include "cube.hh"
#include "../gauleg.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstgauleg()
{
    vector<double> x, w;
    gauleg(0.0, 1.0, x, w, 3);
    double sum = 0.0;
    for (unsigned i=0; i<3; ++i)
    {
        sum += square(square(x[i]))*w[i];
    }
    if (!soft_equiv(sum, 0.2))
    {
        FAILMSG("gauleg NOT accurate");
    }
    else
    {
        PASSMSG("gauleg accurate");
    }        
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
        {
            std::cout << argv[0] << ": version " 
                      << rtt_dsxx::release() 
                      << std::endl;
            return 0;
        }

    try
    {
        tstgauleg();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstgauleg, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstgauleg, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
              << std::endl;
    if (rtt_quadrature_test::passed) 
    {
        std::cout << "**** tstgauleg Test: PASSED" 
                  << std::endl;
    }
    std::cout <<     "*********************************************" 
              << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstgauleg." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstgauleg.cc
//---------------------------------------------------------------------------//
