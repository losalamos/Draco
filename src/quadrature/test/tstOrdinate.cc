//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tstOrdinate.cc
 * \author Kelly Thompson
 * \date   Tue June 20 14:25 2006
 * \brief  Unit test for Ordinate class.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include "special_functions/Ylm.hh"
#include "mesh_element/Geometry.hh"
#include "../QuadCreator.hh"
#include "../Quadrature.hh"
#include "../Q1DGaussLeg.hh"
#include "../QuadServices.hh"
#include "ds++/Release.hh"
#include "../Ordinate.hh"

using namespace rtt_quadrature;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_ordinate_ctor( rtt_dsxx::UnitTest & ut )
{
    using rtt_units::PI;
    
    { // Test default constructor.
        Ordinate Omega;
        if( soft_equiv( Omega.mu(), 0.0 ) &&
            soft_equiv( Omega.eta(), 0.0 ) &&
            soft_equiv( Omega.xi(), 0.0 ) &&
            soft_equiv( Omega.wt(), 0.0 ) )
            ut.passes("Default constructor initializes data to zero.");
        else
            ut.failure("Default constructor did not initialize data to zero.");
    }
    { // Test 1D constructor
        double const sqrtThird(std::sqrt(1.0/3.0)), wt(PI/2);
        Ordinate Omega( sqrtThird, wt );
        if( soft_equiv( Omega.mu(), sqrtThird ) &&
            soft_equiv( Omega.wt(), wt ) )
            ut.passes("Constructor initializes data to correct values.");
        else
            ut.failure("Constructor did not initialize data to correct values.");
    }
    { // Test 3D constructor
        double const sqrtThird(std::sqrt(1.0/3.0)), wt(PI/2);
        Ordinate Omega( sqrtThird, sqrtThird, sqrtThird, wt );
        if( soft_equiv( Omega.mu(), sqrtThird ) &&
            soft_equiv( Omega.eta(), sqrtThird ) &&
            soft_equiv( Omega.xi(), sqrtThird ) &&
            soft_equiv( Omega.wt(), wt ) )
            ut.passes("Constructor initializes data to correct values.");
        else
            ut.failure("Constructor did not initialize data to correct values.");
    }
    { // Test 3D constructor
        double const sqrtThird(std::sqrt(1.0/3.0)), wt(PI/2);
        Ordinate Omega( sqrtThird, sqrtThird, sqrtThird, wt );
        if( soft_equiv( Omega.mu(), sqrtThird ) &&
            soft_equiv( Omega.eta(), sqrtThird ) &&
            soft_equiv( Omega.xi(), sqrtThird ) &&
            soft_equiv( Omega.wt(), wt ) )
            ut.passes("Constructor initializes data to correct values.");
        else
            ut.failure("Constructor did not initialize data to correct values.");

        if (Omega == Ordinate( sqrtThird, sqrtThird, sqrtThird, wt ))
        {
            ut.passes("Ordinate tests equal to self");
        }
        else
        {
            ut.failure("Ordinate does NOT test equal to self");
        }
        if (Omega == Ordinate( 1.9, sqrtThird, sqrtThird, wt ))
        {
            ut.failure("Ordinate does NOT fail equality test on mod mu");
        }
        else
        {
            ut.passes("Ordinate doest not test equal on mod mu");
        }
        if (Omega == Ordinate( sqrtThird, 1.9, sqrtThird, wt ))
        {
            ut.failure("Ordinate does NOT fail equality test on mod eta");
        }
        else
        {
            ut.passes("Ordinate doest not test equal on mod eta");
        }
        if (Omega == Ordinate( sqrtThird, sqrtThird, 1.9, wt ))
        {
            ut.failure("Ordinate does NOT fail equality test on mod xi");
        }
        else
        {
            ut.passes("Ordinate doest not test equal on mod xi");
        }
        if (Omega == Ordinate( sqrtThird, sqrtThird, sqrtThird, 0.2 ))
        {
            ut.failure("Ordinate does NOT fail equality test on mod wt");
        }
        else
        {
            ut.passes("Ordinate doest not test equal on mod wt");
        }
    }
    return;
}

//---------------------------------------------------------------------------//

void test_ordinate_sort( UnitTest & ut )
{
    using namespace rtt_quadrature;
    
    // Generate a quadrature set.
    QuadCreator qc;
    SP< Quadrature const > const spQ( qc.quadCreate( QuadCreator::GaussLeg, 8 ) );
    int const numOrdinates( spQ->getNumOrdinates() );
    vector< Ordinate > ordinates;
    for( int i=0; i<numOrdinates; ++i )
        ordinates.push_back( Ordinate( spQ->getMu(i), spQ->getWt(i) ) );

    std::sort(ordinates.begin(), ordinates.end(), Ordinate::SnCompare );

    bool sorted(true);
    for( int i=0; i<numOrdinates-1; ++i )
    {
        if( ordinates[i+1].mu() < ordinates[i].mu() )
        {
            sorted = false;
            break;
        }
    }
    if( sorted )
        ut.passes("Ordinates were sorted correctly using the comparator SnCompare().");
    else
        ut.failure( "Ordinates not sorted correctly when using the comparator SnCompare().");
    
    return;
}

//---------------------------------------------------------------------------//

void test_create_ordinate_set( UnitTest & ut )
{
    using namespace rtt_quadrature;
    
    // Generate a quadrature set.
    QuadCreator qc;
    SP< Quadrature const > const spQ( qc.quadCreate( QuadCreator::GaussLeg, 8 ) );
    int const numOrdinates( spQ->getNumOrdinates() );

    // Call the function that we are testing.
    int const dim( 1 );
    OrdinateSet const ordinate_set( spQ, rtt_mesh_element::CARTESIAN, dim, false);
    vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();

    // Check the result
    bool looksGood(true);
    for( int i=0; i<numOrdinates; ++i )
    {
        if( ! soft_equiv( ordinates[i].mu(), spQ->getMu(i) ) )
        {
            looksGood = false;
            break;
        }
        if( ! soft_equiv( ordinates[i].wt(), spQ->getWt(i) ) )
        {
            looksGood = false;
            break;
        }
    }
    if( looksGood )
        ut.passes("Create_Ordinate_Set works!");
    else
        ut.passes("Create_Ordinate_Set failed for 1D S8.");

    // test accessor
    SP< Quadrature const > const spQ2( ordinate_set.getQuadrature() );
    if( spQ == spQ2 )
        ut.passes("Quadrature sets match.");
    else
        ut.failure("Quadrature sets do not match.");

    {
        // 1d set for 2d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym2D, 8 ) );

        {
            // Call the function that we are testing.
            int const dim( 1 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::AXISYMMETRIC,
                                            dim,
                                            false);
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim,
                                            false);
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 10 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim,
                                            false );
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 12 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim ,
                                            false);
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 14 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim,
                                            false); 

            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 16 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim,
                                            false );
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 18 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim,
                                            false );
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 20 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim,
                                            false );
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 22 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim ,
                                            false);
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    {
        // 3d set for 3d quadrature
        SP< Quadrature const > const
            spQ( qc.quadCreate( QuadCreator::LevelSym, 24 ) );

        {
            // Call the function that we are testing.
            int const dim( 3 );
            OrdinateSet const ordinate_set( spQ,
                                            rtt_mesh_element::CARTESIAN,
                                            dim ,
                                            false);
            
            vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();
            double sum = 0;
            for (unsigned i=0; i<ordinates.size(); ++i)
            {
                sum += ordinates[i].wt();
            }
            if (soft_equiv(sum, ordinate_set.getNorm()))
                ut.passes("norm matches.");
            else
                ut.failure("norm does NOT match");
        }
    }
    
    return;
}

//---------------------------------------------------------------------------//
void test_Y( UnitTest & ut)
{
    using namespace rtt_quadrature;
    
    // Generate a quadrature set.
    QuadCreator qc;
    int const quadOrder(2);
    SP< Quadrature const > const spQ( qc.quadCreate( QuadCreator::LevelSym2D, quadOrder ) );
    double const sumwt(     spQ->getNorm()      );

    // Call the function that we are testing.
    int const dim( 2 );
    OrdinateSet const ordinate_set( spQ, rtt_mesh_element::CARTESIAN, dim, false );
    vector<Ordinate> const &ordinates = ordinate_set.getOrdinates();

    if (ordinate_set.getNorm() != sumwt)
    {
        ut.failure("did NOT get right norm");
    }
    else
    {
        ut.passes("got right norm");
    }

    for( int ell=0; ell < 3; ++ell )
    {
        int k=-1*ell;
        for( ; k<=ell ; ++k )
        {
            unsigned m=0;
            double sph( Ordinate::Y(ell,k, ordinates[m], sumwt ) );
            double phi( QuadServices::compute_azimuthalAngle( ordinates[m].mu(), ordinates[m].eta(), ordinates[m].xi() ) );
            
            double sfYlm( rtt_sf::galerkinYlk(ell,
                                              k,
                                              ordinates[m].xi(),
                                              phi,
                                              sumwt ) );
            
            if( soft_equiv( sfYlm, sph ) )
            {
                std::ostringstream msg;
                msg << "Y(l,k,Omega,sumwt) == galerkinYlk(l,k,xi,phi,sumwt) "
                    "for l=" << ell
                    << " and k=" << k << "." << std::endl;
                ut.passes( msg.str() );
            }
            else
            {
                std::ostringstream msg;
                msg << "Y(l,k,Omega,sumwt) != galerkinYlk(l,k,xi,phi,sumwt) "
                    "for l=" << ell
                    << " and k=" << k << ".\n"
                    << "\tFound           Y(" << ell << "," << k
                    << ",Omega,sumwt) = " << sph << "\n"
                    << "\tFound galerkinYlk(" << ell << "," << k
                    << ",xi,phi,sumwt) = " << sfYlm << "\n"
                    << "\tFound   phi = " << phi << "\n"
                    << "\tFound atan2 = " << std::atan2(ordinates[m].eta(),
                                                        ordinates[m].mu()) 
                    << std::endl;
                ut.failure( msg.str() );
            }
        }
    }
    return;
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;
    
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member
        // function setTestName).
        rtt_dsxx::ScalarUnitTest ut( argc, argv, release );
        test_ordinate_ctor(ut);
        test_ordinate_sort(ut);
        test_create_ordinate_set(ut);
        test_Y(ut);
    }
    catch( rtt_dsxx::assertion &err )
    {
        cout << "ERROR: While testing " << argv[0] << ", " << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown" << endl;
        return 1;
    }

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tQuadServices.cc
//---------------------------------------------------------------------------//
