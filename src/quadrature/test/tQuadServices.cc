//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tQuadServices.cc
 * \author Kelly Thompson
 * \date   Mon Nov 8 10:48 2004
 * \brief  Quadrature Services.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
//#include <cmath>
#include <sstream>
#include <string>

#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include "quadrature_test.hh"
#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "../QuadServices.hh"
#include "../Release.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_quad_services_with_1D_S2_quad()
{   
    using rtt_quadrature::QuadCreator;
    using rtt_quadrature::Quadrature;
    using rtt_quadrature::QuadServices;
    using rtt_dsxx::SP;
    using rtt_dsxx::soft_equiv;
    
    using std::cout;
    using std::endl;
    using std::string;
    using std::vector;
    using std::ostringstream;

    //----------------------------------------
    // Setup Quadrature set
    
    // create an object that is responsible for creating quadrature objects.
    // QuadCreator QuadratureCreator;
    
    // we will only look at S2 Sets in this test.
    const size_t sn_ord_ref( 2                   );
    const string qname_ref ( "1D Gauss Legendre" );
    const size_t n_ang_ref ( 2                   );
    
    // Banner
    cout << "\nTesting the "  << qname_ref << "S2 quadrature set." << endl;
    
    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::GaussLeg, sn_ord_ref ); 
    
    // print the name of the quadrature set that we are testing.
    const string qname   ( spQuad->name()         );
    const size_t sn_order( spQuad->getSnOrder()   );
    const size_t numAngles( spQuad->getNumAngles() );
    
    // check basic quadrature setup.
    if( sn_order != sn_ord_ref ) 
    {
	FAILMSG("Found incorrect Sn Order.");
    }
    else 
    {
	PASSMSG("Found correct Sn Order.");
    }
    if( numAngles != n_ang_ref  )
    {
	FAILMSG("Found incorrect number of angles.");
    }
    else 
    {
	PASSMSG("Found correct number of angles.");
    }
    if( qname != qname_ref  )
    {
	cout << qname << endl;
	FAILMSG("Found incorrect name of quadrature set.");
    }
    else 
    {
	PASSMSG("Found correct name of quadrature set.");
    }
    
    // Print a table
    spQuad->display();

    //----------------------------------------
    // Setup QuadServices object
    
    QuadServices qs( spQuad );
    
    vector<double> const M( qs.getM() );
    unsigned const numMoments( qs.getNumMoments() );

    std::vector< unsigned > dims;
    dims.push_back( numAngles );
    dims.push_back( numMoments );
    
    qs.print_matrix( "Mmatrix", M, dims );

    //----------------------------------------
    //           2l+1
    // M_{n,m} = ----- * C_n * Y_n( Omega_m )
    //           sumwt
    //
    // k is always 0 for 1D ==> c_n == 1 for 1D
    //
    // M(0,0) = ( 1/2 ) * ( 1 ) * P(0,0)(mu_m)*cos(0)
    //        = 1/2 * P(0,0) = 1/2
    // M(0,m) = 1/2
    //----------------------------------------

    { // scope for testing n=0
	unsigned n(0);
	for( size_t m=0; m<numAngles; ++m )
	{
	    if( soft_equiv( M[ n + m*numMoments ], 0.5 ) )
	    {
		ostringstream msg;
		msg << "M[" << n << "," << m << "] has the expected value." << endl;
		    PASSMSG( msg.str() );
	    }
	    else
	    {		
		ostringstream msg;
		msg << "M[" << n << "," << m 
		    << "] does not have the expected value." << endl
		    << "\tFound M[" << n << "," << m << "] = " 
		    << M[ n + m*numMoments ] << ", but was expecting 0.5" << endl; 
		FAILMSG( msg.str() );		
	    }
	}
    } // end scope for testing n=0

    // When n=1 ==> (l,k) = (1,0)
    // M(1,m) = (3/2) * (1) *P(1,0)(mu_m)*cos(0)
    // M(1,m) = (3/2) * mu_m

    { // scope for testing n=1
	unsigned n(1);
	std::vector< double > const mu( spQuad->getMu() );
	for( size_t m=0; m<numAngles; ++m )
	{
	    if( soft_equiv( M[ n + m*numMoments ], 3*mu[m]/2 ) )
	    {
		ostringstream msg;
		msg << "M[" << n << "," << m << "] has the expected value." << endl;
		    PASSMSG( msg.str() );
	    }
	    else
	    {		
		ostringstream msg;
		msg << "M[" << n << "," << m 
		    << "] does not have the expected value." << endl
		    << "\tFound M[" << n << "," << m << "] = " 
		    << M[ n + m*numMoments ] << ", but was expecting " 
		    << 3*mu[m]/2 << "." << endl; 
		FAILMSG( msg.str() );		
	    }
	}
    } // end scope for testing n=1

    //-----------------------------------//

    vector<double> const D( qs.getD() );
    qs.print_matrix( "Dmatrix", D, dims );

    // The first row of D should contain the quadrature weights.
    {
	unsigned n(0);
	std::vector< double > const wt( spQuad->getWt() );
	for( size_t m=0; m<numAngles; ++m )
	{
	    if( soft_equiv( D[ m + n*numAngles ], wt[m] ) )
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numAngles ] 
		    << " matched the expected value." << endl;
		PASSMSG( msg.str() );
	    }
	    else
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numAngles ] 
		    << " did not match the expected value of " 
		    << wt[m] << "." << endl;
		FAILMSG( msg.str() );
	    }
	}
    }

    return;
}

void test_quad_services_with_1D_S8_quad()
{   
    using rtt_quadrature::QuadCreator;
    using rtt_quadrature::Quadrature;
    using rtt_quadrature::QuadServices;
    using rtt_dsxx::SP;
    using rtt_dsxx::soft_equiv;
    
    using std::cout;
    using std::endl;
    using std::string;
    using std::vector;
    using std::ostringstream;

    //----------------------------------------
    // Setup Quadrature set
    
    // create an object that is responsible for creating quadrature objects.
    // QuadCreator QuadratureCreator;
    
    // we will only look at S2 Sets in this test.
    const size_t sn_ord_ref( 8                   );
    const string qname_ref ( "1D Gauss Legendre" );
    const size_t n_ang_ref ( 8                   );
    
    // Banner
    cout << "\nTesting the "  << qname_ref << " S8 quadrature set." << endl;
    
    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::GaussLeg, sn_ord_ref ); 
    
    // print the name of the quadrature set that we are testing.
    const string qname   ( spQuad->name()         );
    const size_t sn_order( spQuad->getSnOrder()   );
    const size_t numAngles( spQuad->getNumAngles() );
    
    // check basic quadrature setup.
    if( sn_order != sn_ord_ref ) 
    {
	FAILMSG("Found incorrect Sn Order.");
    }
    else 
    {
	PASSMSG("Found correct Sn Order.");
    }
    if( numAngles != n_ang_ref  )
    {
	FAILMSG("Found incorrect number of angles.");
    }
    else 
    {
	PASSMSG("Found correct number of angles.");
    }
    if( qname != qname_ref  )
    {
	cout << qname << endl;
	FAILMSG("Found incorrect name of quadrature set.");
    }
    else 
    {
	PASSMSG("Found correct name of quadrature set.");
    }
    
    // Print a table
    spQuad->display();

    //----------------------------------------
    // Setup QuadServices object
    
    QuadServices qs( spQuad );
    
    vector<double> const M( qs.getM() );
    unsigned const numMoments( qs.getNumMoments() );
    
    //           2l+1
    // M_{n,m} = ----- * C_n * Y_n( Omega_m )
    //           sumwt
    //
    // k is always 0 for 1D ==> c_n == 1 for 1D
    //
    // M(0,0) = ( 1/2 ) * ( 1 ) * P(0,0)(mu_m)*cos(0)
    //        = 1/2 * P(0,0) = 1/2
    // M(0,m) = 1/2

    { // scope for testing n=0
	unsigned n(0);
	for( size_t m=0; m<numAngles; ++m )
	{
	    if( soft_equiv( M[ n + m*numMoments ], 0.5 ) )
	    {
		ostringstream msg;
		msg << "M[" << n << "," << m << "] has the expected value." << endl;
		    PASSMSG( msg.str() );
	    }
	    else
	    {		
		ostringstream msg;
		msg << "M[" << n << "," << m 
		    << "] does not have the expected value." << endl
		    << "\tFound M[" << n << "," << m << "] = " 
		    << M[ n + m*numMoments ] << ", but was expecting 0.5" << endl; 
		FAILMSG( msg.str() );		
	    }
	}
    } // end scope for testing n=0

    // When n=1 ==> (l,k) = (1,0)
    // M(1,m) = (3/2) * (1) *P(1,0)(mu_m)*cos(0)
    // M(1,m) = (3/2) * mu_m

    { // scope for testing n=1
	unsigned n(1);
	std::vector< double > const mu( spQuad->getMu() );
	for( size_t m=0; m<numAngles; ++m )
	{
	    if( soft_equiv( M[ n + m*numMoments ], 3*mu[m]/2 ) )
	    {
		ostringstream msg;
		msg << "M[" << n << "," << m << "] has the expected value." << endl;
		    PASSMSG( msg.str() );
	    }
	    else
	    {		
		ostringstream msg;
		msg << "M[" << n << "," << m 
		    << "] does not have the expected value." << endl
		    << "\tFound M[" << n << "," << m << "] = " 
		    << M[ n + m*numMoments ] << ", but was expecting " 
		    << 3*mu[m]/2 << "." << endl; 
		FAILMSG( msg.str() );		
	    }
	}
    } // end scope for testing n=1

    return;
}

// //---------------------------------------------------------------------------//
// /*! 
//  * \brief test_legendre_poly
//  * 
//  * Test a function that evaluates legendre polynomials.
//  */
// void test_legendre_poly()
// {
//     using rtt_dsxx::soft_equiv;
//     using rtt_dsxx::SP;
//     using rtt_quadrature::QuadCreator;
//     using rtt_quadrature::Quadrature;
//     using rtt_quadrature::QuadServices;

//     // Create a quadrature set from a temporary instance of a
//     // QuadratureCreator factory object.
//     SP< const Quadrature > spQuad;
//     spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, 2 ); 

//     // Create a QuadServices object so we can test its member functions.
//     QuadServices qs( spQuad );

//     // P_{0,0}(x) == 1.0
//     double reference( 1.0 );
//     unsigned k(0), ell(0);
//     double value( qs.legendre_polynomial( k, ell, 1.0 ) );

//     if( soft_equiv(reference,value) )
//     {
// 	PASSMSG("Correct computation of P00(x).");
//     }
//     else
//     {
// 	FAILMSG("Failed to compute P00(x).");
//     }
    
//     return;
// }


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
     using rtt_quadrature::QuadCreator;
     using rtt_quadrature::Quadrature;
     using rtt_quadrature::QuadServices;
     using rtt_dsxx::SP;
 //     using rtt_dsxx::soft_equiv;
    
     using std::cout;
     using std::endl;
     using std::string;
     using std::vector;

     //----------------------------------------
     // Setup Quadrature set

     // create an object that is responsible for creating quadrature objects.
     // QuadCreator QuadratureCreator;
    
     // we will only look at S2 Sets in this test.
     const size_t sn_ord_ref( 2                    );
     const string qname_ref ( "3D Level Symmetric" );
     const size_t n_ang_ref ( 8                    );

     // Banner
     cout << "Testing the "  << qname_ref << " quadrature set." << endl;

     // Create a quadrature set from a temporary instance of a
     // QuadratureCreator factory object.
     SP< const Quadrature > spQuad;
     spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, sn_ord_ref ); 
    
     // print the name of the quadrature set that we are testing.
     const string qname   ( spQuad->name()         );
     const size_t sn_order( spQuad->getSnOrder()   );
     const size_t numAngles   ( spQuad->getNumAngles() );

     // check basic quadrature setup.
     if( sn_order != sn_ord_ref ) ITFAILS;
     if( numAngles    != n_ang_ref  ) ITFAILS;
     if( qname    != qname_ref  ) ITFAILS;

     // Print a table
     spQuad->display();

     //----------------------------------------
     // Setup QuadServices object

     QuadServices qs( spQuad );

     vector<double> const M( qs.getM() );
     unsigned const numMoments( M.size()/numAngles );

     return;
} // end of quadrature_test

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::string;
    using std::cout;
    using std::endl;

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_quadrature::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	test_quad_services_with_1D_S2_quad();
	test_quad_services_with_1D_S8_quad();
	// quadrature_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tQuadServices, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_quadrature_test::passed) 
    {
        cout << "**** tQuadServices Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tQuadServices." << endl;

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tQuadServices.cc
//---------------------------------------------------------------------------//
