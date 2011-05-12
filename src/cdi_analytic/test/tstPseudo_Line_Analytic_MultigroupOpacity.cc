//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstPseudo_Line_Analytic_MultigroupOpacity.cc
 * \author Kent G. Budge
 * \date   Tue Apr  5 09:01:03 2011
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "../Pseudo_Line_Analytic_MultigroupOpacity.hh"
#include "parser/Constant_Expression.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_cdi_analytic;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

Pseudo_Line_Analytic_MultigroupOpacity::Averaging const NONE =
    Pseudo_Line_Analytic_MultigroupOpacity::NONE;

Pseudo_Line_Analytic_MultigroupOpacity::Averaging const ROSSELAND =
    Pseudo_Line_Analytic_MultigroupOpacity::ROSSELAND;

Pseudo_Line_Analytic_MultigroupOpacity::Averaging const PLANCK =
    Pseudo_Line_Analytic_MultigroupOpacity::PLANCK;

void tstPseudo_Line_Analytic_MultigroupOpacity(UnitTest &ut)
{
    unsigned const NG = 2000; // 300 for full resolution
    unsigned const number_of_lines = 200;
    unsigned const number_of_edges = 10;
    unsigned seed = 1;

    SP<Expression const> const continuum(new Constant_Expression(1,1.0e-2));;
    double const peak = 1e1;
    double const width = 0.002; // keV
    double const edge_ratio = 10.0; 
    double const emax = 10.0; // keV
    double const emin = 0.0;
    
    vector<double> group_bounds(NG+1);
    for (unsigned i=0; i<=NG; i++)
    {
        group_bounds[i] = i*(emax-emin)/NG + emin;
    }

    {
        Pseudo_Line_Analytic_MultigroupOpacity model(group_bounds,
                                                     rtt_cdi::ABSORPTION,
                                                     continuum,
                                                     number_of_lines,
                                                     peak,
                                                     width,
                                                     number_of_edges,
                                                     edge_ratio,
                                                     1.0,
                                                     0.0,
                                                     emin,
                                                     emax,
                                                     NONE,
                                                     seed);
        
        ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");
        
        vector<double> sigma = model.getOpacity(1.0, 1.0);
        
        ofstream out("pseudo_none.dat");
        for (unsigned g=0; g<NG; ++g)
        {
            out << ((g+0.5)*(emax-emin)/NG + emin) << ' ' << sigma[g] << endl;
        }
        
        ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
    }

    {
        Pseudo_Line_Analytic_MultigroupOpacity model(group_bounds,
                                                     rtt_cdi::ABSORPTION,
                                                     continuum,
                                                     number_of_lines,
                                                     peak,
                                                     width,
                                                     number_of_edges,
                                                     edge_ratio,
                                                     1.0,
                                                     0.0,
                                                     emin,
                                                     emax,
                                                     ROSSELAND,
                                                     seed);
        
        ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");
        
        vector<double> sigma = model.getOpacity(1.0, 1.0);
        
        ofstream out("pseudo_rosseland.dat");
        for (unsigned g=0; g<NG; ++g)
        {
            out << ((g+0.5)*(emax-emin)/NG + emin) << ' ' << sigma[g] << endl;
        }
        
        ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
    }

    {
        Pseudo_Line_Analytic_MultigroupOpacity model(group_bounds,
                                                     rtt_cdi::ABSORPTION,
                                                     continuum,
                                                     number_of_lines,
                                                     peak,
                                                     width,
                                                     number_of_edges,
                                                     edge_ratio,
                                                     1.0,
                                                     0.0,
                                                     emin,
                                                     emax,
                                                     PLANCK,
                                                     seed);
        
        ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");
        
        vector<double> sigma = model.getOpacity(1.0, 1.0);
        
        ofstream out("pseudo_planck.dat");
        for (unsigned g=0; g<NG; ++g)
        {
            out << ((g+0.5)*(emax-emin)/NG + emin) << ' ' << sigma[g] << endl;
        }
        
        ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstPseudo_Line_Analytic_MultigroupOpacity(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstPseudo_Line_Analytic_MultigroupOpacity, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstPseudo_Line_Analytic_MultigroupOpacity, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstPseudo_Line_Analytic_MultigroupOpacity.cc
//---------------------------------------------------------------------------//
