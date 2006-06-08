//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/test/tstYlm.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Soft_Equivalence.hh"
#include "ds++/ScalarUnitTest.hh"
#include "units/PhysicalConstants.hh"

#include "../Ylm.hh"
#include "../Release.hh"

using namespace std;
using namespace rtt_sf;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

void tstYlm( rtt_dsxx::UnitTest & ut )
{
    using rtt_units::PI;
    using rtt_dsxx::soft_equiv;

    double const theta = 0.4, phi = 0.21;
    if (soft_equiv(Ylm(2,-2,theta,phi),
                   -sqrt(15/(16*PI))*sin(theta)*sin(theta)*sin(2*phi)))
    {
        ut.passes("Ylm function returned correct double");
    }
    else
    {
        ut.failure("Ylm function did NOT return correct double.");
    }
    if (soft_equiv(Ylm(2,-1,theta,phi),
                   -sqrt(15/(4*PI))*sin(theta)*cos(theta)*cos(phi)))
    {
        ut.passes("Ylm function returned correct double");
    }
    else
    {
        ut.failure("Ylm function did NOT return correct double.");
    }
    if (soft_equiv(Ylm(2,0,theta,phi),
                   sqrt(5/(16*PI))*(3*cos(theta)*cos(theta)-1)))
    {
        ut.passes("Ylm function returned correct double");
    }
    else
    {
        ut.failure("Ylm function did NOT return correct double.");
    }
    if (soft_equiv(Ylm(2,1,theta,phi),
                   -sqrt(15/(4*PI))*sin(theta)*cos(theta)*sin(phi)))
    {
        ut.passes("Ylm function returned correct double");
    }
    else
    {
        ut.failure("Ylm function did NOT return correct double.");
    }
    if (soft_equiv(Ylm(2,2,theta,phi),
                   -sqrt(15/(16*PI))*sin(theta)*sin(theta)*cos(2*phi)))
    {
        ut.passes("Ylm function returned correct double");
    }
    else
    {
        ut.failure("Ylm function did NOT return correct double.");
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        rtt_dsxx::ScalarUnitTest ut( argc, argv, release );
        tstYlm( ut );
        ut.status();
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown on processor " << endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of testYlm.cc
//---------------------------------------------------------------------------//
