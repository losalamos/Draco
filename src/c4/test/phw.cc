//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/phw.cc
 * \author Kelly Thompson
 * \date   Tue Jun  6 15:03:08 2006
 * \brief  Parallel application used by the unit test for tstApplicationUnitTest.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include "../global.hh"

using namespace rtt_c4;

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    initialize(argc,argv);
    if( argc > 1 )
    {
        std::cout << "Hello, world!" << std::endl;
        finalize();
        return 0;
    }
    else
    {
        std::cout << "Goodbye cruel world..." << std::endl;
        finalize();
        return 1;
    }
    finalize();
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of phw.cc
//---------------------------------------------------------------------------//
