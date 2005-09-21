//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tQuadrature.cc
 * \author Kelly Thompson
 * \date   Tue Mar 26 12:36:41 2002
 * \brief  quadrature package test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//


#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "../Q2DLevelSym.hh"
#include "../Q3DLevelSym.hh"
#include "../Release.hh"

#include "quadrature_test.hh"

using namespace std;

using rtt_quadrature::QuadCreator;
using rtt_quadrature::Quadrature;
using rtt_dsxx::SP;

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
void quadrature_test()
{
    using rtt_dsxx::soft_equiv;

    // double precesion values will be tested for correctness against this
    // tolerance. 
    const double TOL = 1.0e-10; 

    // create an object that is responsible for creating quadrature objects.
    QuadCreator QuadratureCreator;
    
    // we will only look at S4 Sets in this test.
    const int sn_order = 4;

    // total number of quadrature sets to be tested.
    const int nquads = 5;

    // Declare an enumeration object that specifies the Quadrature set to be
    // tested.

    // Quadrature sets to be tested:
    //
    // #   Qid         Description
    // -   --------    ------------------
    // 0   GaussLeg    1D Gauss Legendre
    // 1   Lobatto     1D Labotto
    // 2   DoubleGauss 1D DoubleGauss
    // 3   LevelSym2D  2D Level Symmetric
    // 4   LevelSym    3D Level Symmetric

    QuadCreator::Qid qid[nquads] = { QuadCreator::GaussLeg,
                                     QuadCreator::Lobatto,
                                     QuadCreator::DoubleGauss,
				     QuadCreator::LevelSym2D,
				     QuadCreator::LevelSym };

    // mu0 holds mu for the first direction for each quadrature set tested.
    double mu0[nquads] = { 0.8611363116, 1.0, 0.7886751346,
			  -0.350021174581541, -0.350021174581541 };
    
    SP< const Quadrature > spQuad;

    // loop over quadrature types to be tested.

    for ( int ix = 0; ix < nquads; ++ix ) {
	
	// Verify that the enumeration value matches its int value.
	if ( qid[ix] != ix ) 
	{
	    FAILMSG("Setting QuadCreator::Qid enumeration failed.");
	    break;
	}
	else 
	{
	    // Instantiate the quadrature object.
	    spQuad = QuadratureCreator.quadCreate( qid[ix], sn_order ); 

	    // print the name of the quadrature set that we are testing.
	    string qname = spQuad->name();
	    cout << "\nTesting the "  << qname
		 << " quadrature set." << endl
                 << "   Sn Order         = " << spQuad->getSnOrder() << endl
                 << "   Number of Angles = " << spQuad->getNumAngles() << endl
                 << "   Dimensionality   = " << spQuad->dimensionality()
                 << endl << endl; 

	    // Extra tests for improving coverage analysis
	    if( qid[ix] == QuadCreator::GaussLeg)
	    {
		double const expected_sumwt( 2.0 );
		if( soft_equiv( spQuad->iDomega(), expected_sumwt ) )
		{
		    PASSMSG("Sumwt for GaussLeg quad is 2.0, as expected.");
		}
		else
		{
		    ostringstream msg;
		    msg << "Unexpected value returned from spQuad->iDomega()."
			<< endl
			<< "Expected iDomega() to return " << expected_sumwt
			<< ", but instead it returned" << spQuad->iDomega()
			<< "." << endl;
		    FAILMSG( msg.str() );
		}
	    }
	    else if( qid[ix] == QuadCreator::LevelSym )
	    {
		vector<double> const veta( spQuad->getEta() );
		vector<double> const vxi( spQuad->getXi() );

		// ensure that the length of the discrete ordinate vector is 1.
		for( size_t m=0; m<spQuad->getNumAngles(); ++m )
		{
		    double const mu( spQuad->getMu(m) );
		    double const eta( spQuad->getEta(m) );
		    double const xi( spQuad->getXi(m) );
		    double const len( mu*mu + eta*eta + xi*xi );
		    if( soft_equiv( veta[m], eta ) )
			PASSMSG( "vector and single eta accessor agree." )
		    else
			FAILMSG( "vector and single eta accessor disagree." )
		    if( soft_equiv( vxi[m], xi ) )
			PASSMSG( "vector and single xi accessor agree." )
		    else
			FAILMSG( "vector and single xi accessor disagree." )
		    if( soft_equiv( len, 1.0 ) )
		    {
			ostringstream msg;
			msg << "Length of direction ordinate " << m
			    << " has correct value of 1.0" << endl;
			PASSMSG(msg.str());		       
		    }
		    else
		    {
			ostringstream msg;
			msg << "Length of direction ordinate " << m
			    << " does not have the correct value of 1.0" 
			    << endl
			    << "It was found to have length " << len
			    << ", which is incorrect." << endl;
			FAILMSG(msg.str());
		    }
		}
		
	    }
            else if( qid[ix] == QuadCreator::Lobatto )
            {
                int const expected_dim(1);
                cout << expected_dim << endl;
                cout <<  spQuad->dimensionality() << endl;
                if( spQuad->dimensionality() == expected_dim )
                {
                    PASSMSG("Dimensionality of Lobatto quadrature set is 1.");
                }
                else
                {
                    FAILMSG("Dimensionality of Lobatto quadrature set is incorrect.");
                }
            }

	    // If the object was constructed sucessfully then we continue
	    // with the tests.
	    if ( ! spQuad )
		FAILMSG("QuadCreator failed to create a new quadrature set.")
	    else 
	    {
		// get the mu vector
		vector<double> mu = spQuad->getMu();
		// get the omega vector for direction m=1.
		vector<double> omega_1 = spQuad->getOmega(1);
		// get the total omega object.
		vector< vector<double> > omega = spQuad->getOmega();
		// compare values.
		if ( mu.size() != spQuad->getNumAngles() )
		    FAILMSG("The direction vector has the wrong length.")
		else if ( fabs( mu[0] + mu0[ix] ) >= TOL ) 
		    FAILMSG("mu[0] has the wrong value.")
		else if ( fabs( mu[1] - omega_1[0] ) >= TOL )
		    FAILMSG("mu[1] != omega_1[0].")
		else if ( fabs( omega[1][0] - omega_1[0] ) >= TOL )
		    FAILMSG("omega[1][0] != omega_1[0].")
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

void Q2DLevelSym_tests()
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
	PASSMSG("Found expected number of levels in quadrature set.");
    }
    else
    {
	ostringstream msg;
	msg << "Found the wrong number of quadrature levels." << endl
	    << "quad.getLevels() returned " << quad.getLevels()
	    << ", but we expected to find " << expected_nlevels << "."
	    << endl;
	FAILMSG( msg.str() );
    }
    return;
} // end of Q2DLevelSym_tests()
//---------------------------------------------------------------------------//

void Q3DLevelSym_tests()
{
    using rtt_dsxx::soft_equiv;
    using rtt_quadrature::Q3DLevelSym;
    using std::ostringstream;
    using std::endl;

    int    const sn_order( 4 );
    double const assigned_sumwt( 1.0 );
    Q3DLevelSym const quad( sn_order, assigned_sumwt );

    size_t const expected_nlevels((sn_order+2)*sn_order/8);
    if( quad.getLevels() == expected_nlevels)
    {
	PASSMSG("Found expected number of levels in quadrature set.");
    }
    else
    {
	ostringstream msg;
	msg << "Found the wrong number of quadrature levels." << endl
	    << "quad.getLevels() returned " << quad.getLevels()
	    << ", but we expected to find " << expected_nlevels << "."
	    << endl;
	FAILMSG( msg.str() );
    }
    
    double const sumwt( 1.0 );
    if( soft_equiv( sumwt, assigned_sumwt ) )
    {
	PASSMSG("Stored sumwt matches assigned value.");
    }
    else
    {
	ostringstream msg;
	msg << "Stored sumwt does not match assigned value as retrieved by iDomega()."
	    << endl
	    << "quad.iDomega() returned " << quad.iDomega()
	    << ", but we expected to find " << assigned_sumwt << "."
	    << endl;
	FAILMSG( msg.str() );
    }

    return;
} // end of Q3DLevelSym_tests()

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_quadrature::release() 
		 << endl;
	    return 0;
	}

    cout << "This is the quadrature package." << endl
	 << "Version " << rtt_quadrature::release() << endl;

    try
    {
	// >>> UNIT TESTS
	quadrature_test();
	Q2DLevelSym_tests();
	Q3DLevelSym_tests();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tQuadrature, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_quadrature_test::passed) 
    {
        cout << "**** tQuadrature Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tQuadrature." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tQuadrature.cc
//---------------------------------------------------------------------------//
