//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   min/test/tstmnbrak.cc
 * \author Kent Budge
 * \date   Tue Aug 26 13:12:30 2008
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "../mnbrak.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_min;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

double func(double const x)
{
    return x*x;
}

//---------------------------------------------------------------------------//
void tstmnbrak(UnitTest &ut)
{
    double ax, bx, cx, fa, fb, fc;

    ax = 1.0;
    bx = 2.0;
    
    mnbrak(ax,
            bx,
            cx,
            fa,
            fb,
            fc,
           func);

    if (min(ax, min(bx, cx)) < 0.0 && max(ax, max(bx, cx)) > 0.0)
    {
        ut.passes("minimum bracketed");
    }
    else
    {
        ut.failure("minimum NOT bracketed");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstmnbrak(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstmnbrak, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstmnbrak, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstmnbrak.cc
//---------------------------------------------------------------------------//
