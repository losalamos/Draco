n//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tQuadrature.cc
 * \brief  quadrature package test.
 * \note   Copyright (C) 2006-2012 LANSLLC All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//


#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "../Q1DGaussLeg.hh"
#include "../Q2DSquareChebyshevLegendre.hh"
#include "../Q2DTriChebyshevLegendre.hh"
#include "../Q3DTriChebyshevLegendre.hh"
#include "../Q2DLevelSym.hh"
#include "../Q3DLevelSym.hh"
#include "../GeneralQuadrature.hh"

#include "units/PhysicalConstants.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define ITFAILS    ut.failure( __LINE__, __FILE__ )

using namespace std;
using namespace rtt_quadrature;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*!
 * \brief Tests the Quadcrator and Quadtrature constructors and access
 * routines. 
 *
 * To add a quadrature to this test the following items must be changed: add
 * new enumeration to Qid[] array.  add new mu[0] value to mu0[] array.
 * verify nquads is set to the correct number of quadrature sets being
 * tested.
 */
void quadrature_test( rtt_dsxx::UnitTest & ut )
{
    using rtt_dsxx::soft_equiv;

    // double precesion values will be tested for correctness against this
    // tolerance. 
    const double TOL = 1.0e-08; 

    // create an object that is responsible for creating quadrature objects.
    QuadCreator QuadratureCreator;
    
    // we will only look at S4 Sets in this test.
    const int sn_order = 4;

    // total number of quadrature sets to be tested.
    const int nquads = 8;

    // Declare an enumeration object that specifies the Quadrature set to be
    // tested.

    // Quadrature sets to be tested:
    //
    // #   Qid         Description
    // --- --------    ------------------
    // 0   TriCL2D     2D Triangular Chebyshev-Legendre
    // 1   SquareCL    2D Square Chebyshev-Legendre
    // 2   GaussLeg    1D Gauss-Legendre
    // 3   Lobatto     1D Lobatto
    // 4   DoubleGauss 1D DoubleGauss
    // 5   LevelSym2D  2D Level Symmetric
    // 6   LevelSym    3D Level Symmetric
    // 7   TriCL       3D Triangular Chebyshev-Legendre

    QuadCreator::Qid qid[nquads] = { QuadCreator::TriCL2D,
                                     QuadCreator::SquareCL,
                                     QuadCreator::GaussLeg,
                                     QuadCreator::Lobatto,
                                     QuadCreator::DoubleGauss,
				     QuadCreator::LevelSym2D,
				     QuadCreator::LevelSym,
                                     QuadCreator::TriCL };

    // mu0 holds mu for the first direction for each quadrature set tested.
    double mu0[nquads] = { -0.359474792477992,
                           -0.469676450658365,
                           -0.861136311594053,
                           -1.0,
                           -0.7886751346,
			   -0.350021174581541,
                           0.350021174581541,
                           -0.861136311594053 };
    
    SP< const Quadrature > spQuad;

    // loop over quadrature types to be tested.

    for ( int ix = 0; ix < nquads; ++ix )
    {
	
	// Verify that the enumeration value matches its int value.
	if ( qid[ix] != ix ) 
	{
	    ut.failure("Setting QuadCreator::Qid enumeration failed.");
	    break;
	}
	else 
	{
	    // Instantiate the quadrature object.
	    spQuad = QuadratureCreator.quadCreate( qid[ix], sn_order ); 

	    // print the name of the quadrature set that we are testing.
	    string qname = spQuad->name();
	    cout << "Testing the "  << qname
		 << " quadrature set."
                 << "\n   Sn Order            = " << spQuad->getSnOrder()
                 << "\n   Number of Ordinates = " << spQuad->getNumOrdinates()
                 << "\n   Dimensionality      = " << spQuad->dimensionality()
                 << "\n"
                 << "\n   Parser Name         = " << spQuad->parse_name()
                 << "\n   Class               = " << spQuad->getClass()
                 << endl;

	    // Extra tests for improving coverage analysis
	    if( qid[ix] == QuadCreator::GaussLeg)
	    {
		double const expected_sumwt( 2.0 );
		if( soft_equiv( spQuad->iDomega(), expected_sumwt ) )
		{
		    ut.passes("Sumwt for GaussLeg quad is 2.0, as expected.");
		}
		else
		{
		    ostringstream msg;
		    msg << "Unexpected value returned from spQuad->iDomega()."
			<< endl
			<< "Expected iDomega() to return " << expected_sumwt
			<< ", but instead it returned" << spQuad->iDomega()
			<< "." << endl;
		    ut.failure( msg.str() );
		}
	    }
	    else if( qid[ix] == QuadCreator::LevelSym )
	    {
		vector<double> const veta( spQuad->getEta() );
		vector<double> const vxi( spQuad->getXi() );

		// ensure that the length of the discrete ordinate vector is 1.
		for( size_t m=0; m<spQuad->getNumOrdinates(); ++m )
		{
		    double const mu( spQuad->getMu(m) );
		    double const eta( spQuad->getEta(m) );
		    double const xi( spQuad->getXi(m) );
		    double const len( mu*mu + eta*eta + xi*xi );
		    if( soft_equiv( veta[m], eta ) )
			ut.passes( "vector and single eta accessor agree." );
		    else
			ut.failure( "vector and single eta accessor disagree." );
		    if( soft_equiv( vxi[m], xi ) )
			ut.passes( "vector and single xi accessor agree." );
		    else
			ut.failure( "vector and single xi accessor disagree." );
		    if( soft_equiv( len, 1.0 ) )
		    {
			ostringstream msg;
			msg << "Length of direction ordinate " << m
			    << " has correct value of 1.0" << endl;
			ut.passes(msg.str());		       
		    }
		    else
		    {
			ostringstream msg;
			msg << "Length of direction ordinate " << m
			    << " does not have the correct value of 1.0" 
			    << endl
			    << "It was found to have length " << len
			    << ", which is incorrect." << endl;
			ut.failure(msg.str());
		    }
		}
		
	    }
            else if( qid[ix] == QuadCreator::Lobatto )
            {
                size_t const expected_dim(1);
                if( spQuad->dimensionality() == expected_dim )
                {
                    ut.passes("Dimensionality of Lobatto quadrature set is 1.");
                }
                else
                {
                    ut.failure("Dimensionality of Lobatto quadrature set is incorrect.");
                }
            }

	    // If the object was constructed sucessfully then we continue
	    // with the tests.
	    if ( ! spQuad )
		ut.failure("QuadCreator failed to create a new quadrature set.");
	    else 
	    {
		// get the mu vector
		vector<double> mu = spQuad->getMu();

		// get the omega vector for direction m=1.
		vector<double> omega_1 = spQuad->getOmega(1);

		// get the total omega object.
		vector< vector<double> > omega = spQuad->getOmega();

		// compare values.
		if ( mu.size() != spQuad->getNumOrdinates() )
		    ut.failure("The direction vector has the wrong length.");
                else if ( fabs(mu[0]-mu0[ix]) > TOL ) 
                    ut.failure("mu[0] has the wrong value.");
                else if ( fabs(mu[1]-omega_1[0]) > TOL )
                    ut.failure("mu[1] != omega_1[0].");
                else if ( fabs(omega[1][0]-omega_1[0]) > TOL )
                    ut.failure("omega[1][0] != omega_1[0].");
		else 
		{
		    spQuad->display();
		    cout << endl << endl; // end of this quadrature type
		}
	    }
	}
    }
    return;
} // end of quadrature_test

//---------------------------------------------------------------------------//

void Q2DLevelSym_tests( rtt_dsxx::UnitTest & ut )
{
    using rtt_quadrature::Q2DLevelSym;
    using std::ostringstream;
    using std::endl;

    int    const sn_order( 4 );
    double const sumwt( 1.0 );
    Q2DLevelSym const quad( sn_order, sumwt );

    size_t const expected_nlevels((sn_order+2)*sn_order/8);
    if( quad.getLevels() == expected_nlevels)
    {
	ut.passes("Found expected number of levels in quadrature set.");
    }
    else
    {
	ostringstream msg;
	msg << "Found the wrong number of quadrature levels." << endl
	    << "quad.getLevels() returned " << quad.getLevels()
	    << ", but we expected to find " << expected_nlevels << "."
	    << endl;
	ut.failure( msg.str() );
    }
    return;
} // end of Q2DLevelSym_tests()
//---------------------------------------------------------------------------//

void Q3DLevelSym_tests( rtt_dsxx::UnitTest & ut )
{
    using rtt_dsxx::soft_equiv;
    using rtt_quadrature::Q3DLevelSym;
    using std::ostringstream;
    using std::endl;

    { // Test low order quadrature sest.
        
        int    const sn_order( 4 );
        double const assigned_sumwt( 1.0 );
        Q3DLevelSym const quad( sn_order, assigned_sumwt );
        
        size_t const expected_nlevels((sn_order+2)*sn_order/8);
        if( quad.getLevels() == expected_nlevels)
        {
            ut.passes("Found expected number of levels in quadrature set.");
        }
        else
        {
            ostringstream msg;
            msg << "Found the wrong number of quadrature levels." << endl
                << "quad.getLevels() returned " << quad.getLevels()
                << ", but we expected to find " << expected_nlevels << "."
                << endl;
            ut.failure( msg.str() );
        }
        
        double const sumwt( 1.0 );
        if( soft_equiv( sumwt, assigned_sumwt ) )
        {
            ut.passes("Stored sumwt matches assigned value.");
        }
        else
        {
            ostringstream msg;
            msg << "Stored sumwt does not match assigned value as retrieved by iDomega()."
                << endl
                << "quad.iDomega() returned " << quad.iDomega()
                << ", but we expected to find " << assigned_sumwt << "."
                << endl;
            ut.failure( msg.str() );
        }
        
        // Test renormalize member function
        {
            // must create a non-const quadrature
            Q3DLevelSym myQuad( sn_order, assigned_sumwt );
            double const wt1a( myQuad.getWt(1) );
            double const fourpi( 4.0*rtt_units::PI );
            myQuad.renormalize( fourpi );
            
            if( ! soft_equiv( myQuad.getNorm(), fourpi     ) )         ITFAILS;
            if( ! soft_equiv( myQuad.getWt(1), fourpi*wt1a ) )         ITFAILS;
        }
    }
    
    { // Test a high order quadrature set
        int    const sn_order( 22 );
        double const assigned_sumwt( 1.0 );
        Q3DLevelSym const quad( sn_order, assigned_sumwt );
        
        size_t const expected_nlevels((sn_order+2)*sn_order/8);
        if( quad.getLevels() == expected_nlevels)
            ut.passes("Found expected number of levels in quadrature set.");
        else
        {
            ostringstream msg;
            msg << "Found the wrong number of quadrature levels." << endl
                << "quad.getLevels() returned " << quad.getLevels()
                << ", but we expected to find " << expected_nlevels << "."
                << endl;
            ut.failure( msg.str() );
        }
        
        double const sumwt( 1.0 );
        if( soft_equiv( sumwt, assigned_sumwt ) )
            ut.passes("Stored sumwt matches assigned value.");
        else
        {
            ostringstream msg;
            msg << "Stored sumwt does not match assigned value as retrieved"
                << " by iDomega()." << endl
                << "quad.iDomega() returned " << quad.iDomega()
                << ", but we expected to find " << assigned_sumwt << "."
                << endl;
            ut.failure( msg.str() );
        }
    }
    return;
} // end of Q3DLevelSym_tests()

//---------------------------------------------------------------------------//
void Q2DSCL_test( rtt_dsxx::UnitTest & ut )
{
    using namespace rtt_dsxx;
    using namespace rtt_quadrature;
    using namespace std;

    cout << "\nLooking at Q2DSquareChebyshevLegendre\n" << endl;
    
    int    const sn_order( 8 );
    double const assigned_sumwt( 1.0 );
    Q2DSquareChebyshevLegendre const quad( sn_order, assigned_sumwt );

    size_t const expected_nlevels(sn_order);
    if( quad.getLevels() == expected_nlevels)
    {
	ut.passes("Found expected number of levels in quadrature set.");
    }
    else
    {
	ostringstream msg;
	msg << "Found the wrong number of quadrature levels." << endl
	    << "quad.getLevels() returned " << quad.getLevels()
	    << ", but we expected to find " << expected_nlevels << "."
	    << endl;
	ut.failure( msg.str() );
    }
    
    double const sumwt( 1.0 );
    if( soft_equiv( sumwt, assigned_sumwt ) )
    {
	ut.passes("Stored sumwt matches assigned value.");
    }
    else
    {
	ostringstream msg;
	msg << "Stored sumwt does not match assigned value as retrieved by iDomega()."
	    << endl
	    << "quad.iDomega() returned " << quad.iDomega()
	    << ", but we expected to find " << assigned_sumwt << "."
	    << endl;
	ut.failure( msg.str() );
    }

    return;
} // end of Q2DSCL_test()

//---------------------------------------------------------------------------//
void Q2DTCL_test( rtt_dsxx::UnitTest & ut )
{
    using namespace rtt_dsxx;
    using namespace rtt_quadrature;
    using namespace std;

    cout << "\nLooking at Q2DTriChebyshevLegendre\n" << endl;
    
    int    const sn_order( 8 );
    double const assigned_sumwt( 2*rtt_units::PI );
    Q2DTriChebyshevLegendre const quad( sn_order, assigned_sumwt );

    size_t const expected_nlevels(sn_order);
    if( quad.getLevels() == expected_nlevels)
    {
	ut.passes("Found expected number of levels in quadrature set.");
    }
    else
    {
	ostringstream msg;
	msg << "Found the wrong number of quadrature levels." << endl
	    << "quad.getLevels() returned " << quad.getLevels()
	    << ", but we expected to find " << expected_nlevels << "."
	    << endl;
	ut.failure( msg.str() );
    }
    
    double const sumwt(quad.iDomega());
    if( soft_equiv( sumwt, assigned_sumwt ) )
    {
	ut.passes("Stored sumwt matches assigned value.");
    }
    else
    {
	ostringstream msg;
	msg << "Stored sumwt does not match assigned value as retrieved by iDomega()."
	    << endl
	    << "quad.iDomega() returned " << quad.iDomega()
	    << ", but we expected to find " << assigned_sumwt << "."
	    << endl;
	ut.failure( msg.str() );
    }

    return;
} // end of Q2DTCL_test()

//---------------------------------------------------------------------------//
void Q3DTCL_test( rtt_dsxx::UnitTest & ut )
{
    using namespace rtt_dsxx;
    using namespace rtt_quadrature;
    using namespace std;

    cout << "\nLooking at Q3DTriChebyshevLegendre\n" << endl;
    
    int    const sn_order( 8 );
    double const assigned_sumwt( 2*rtt_units::PI );
    Q3DTriChebyshevLegendre const quad( sn_order, assigned_sumwt );

    size_t const expected_nlevels((sn_order+2)*sn_order/8);
    if( quad.getLevels() == expected_nlevels)
    {
	ut.passes("Found expected number of levels in quadrature set.");
    }
    else
    {
	ostringstream msg;
	msg << "Found the wrong number of quadrature levels." << endl
	    << "quad.getLevels() returned " << quad.getLevels()
	    << ", but we expected to find " << expected_nlevels << "."
	    << endl;
	ut.failure( msg.str() );
    }
    
    double const sumwt(quad.iDomega());
    if( soft_equiv( sumwt, assigned_sumwt ) )
    {
	ut.passes("Stored sumwt matches assigned value.");
    }
    else
    {
	ostringstream msg;
	msg << "Stored sumwt does not match assigned value as retrieved by iDomega()."
	    << endl
	    << "quad.iDomega() returned " << quad.iDomega()
	    << ", but we expected to find " << assigned_sumwt << "."
	    << endl;
	ut.failure( msg.str() );
    }

    return;
} // end of Q3DTCL_test()

//---------------------------------------------------------------------------//

void tst_general_quadrature( UnitTest & ut )
{
    cout << "\nRunning unit test function: tst_general_quadrature...\n" << endl;
    
    int const snOrder(4);
    double const norm(1.0);
    Q3DLevelSym const refQuad( snOrder, norm );

    // Create a general quadrature
    GeneralQuadrature const quad( snOrder,
                                  norm,
                                  refQuad.getMu(),
                                  refQuad.getEta(),
                                  refQuad.getXi(),
                                  refQuad.getWt(),
                                  refQuad.getLevels(),
                                  refQuad.dimensionality(),
                                  string("GeneralQuadrature"),
                                  refQuad.getClass()); 

    quad.display();
    
    if( quad.getNumOrdinates() != refQuad.getNumOrdinates() )          ITFAILS;
    if( quad.name() != string("GeneralQuadrature") )                   ITFAILS;
    if( quad.dimensionality() != refQuad.dimensionality() )            ITFAILS;
    if( quad.getSnOrder() != refQuad.getSnOrder() )                    ITFAILS;
    if( quad.getLevels() != refQuad.getLevels() )                      ITFAILS;
    if( quad.getClass() != refQuad.getClass() )                        ITFAILS;
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
	quadrature_test(ut);
	Q2DLevelSym_tests(ut);
	Q3DLevelSym_tests(ut);
        Q2DSCL_test(ut);
        Q2DTCL_test(ut);
        Q3DTCL_test(ut);
        tst_general_quadrature(ut);
    }
    catch( rtt_dsxx::assertion &err )
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
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
             << "An unknown exception was thrown" << endl;
        return 1;
    }

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tQuadrature.cc
//---------------------------------------------------------------------------//
