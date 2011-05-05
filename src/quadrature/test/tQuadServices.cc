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
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>

#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/ScalarUnitTest.hh"
#include "special_functions/Factorial.hh"
#include "special_functions/KroneckerDelta.hh"
#include "units/PhysicalConstants.hh"

#include "quadrature_test.hh"
#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "../QuadServices.hh"
#include "ds++/Release.hh"

using namespace rtt_quadrature;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// Legendre Polynomials
double P( int const ell, double const x )
{
    Require( ell >= 0 );
    Require( ell < 7 );
    Require( std::fabs(x) <= 1.0 );

    if( ell == 0 ) return 1.0;
    if( ell == 1 ) return x;
    if( ell == 2 ) return (3.0*x*x-1.0)/2.0;
    if( ell == 3 ) return (5.0*x*x*x - 3.0*x )/2.0;
    if( ell == 4 ) return (35.0*x*x*x*x - 30.0*x*x + 3)/8.0;
    if( ell == 5 ) return (63.0*x*x*x*x*x - 70.0*x*x*x + 15.0*x)/8.0;

    Ensure( ell == 6 );
    return (231.0*x*x*x*x*x*x - 315.0*x*x*x*x +105.0*x*x - 5.0)/16.0;
}

// Associated Legendre Polynomials
double P( unsigned const ell, unsigned const k, double const x )
{
    Require( k <= ell );
    Require( std::fabs(x) <= 1.0 );
    Require( ell < 4 );

    if( ell == 0 ) return 1.0;
    if( ell == 1 && k == 0 ) return x;
    if( ell == 1 && k == 1 ) return -1.0 * std::sqrt(1.0-x*x);
    if( ell == 2 && k == 0 ) return (3.0*x*x-1.0)/2.0;
    if( ell == 2 && k == 1 ) return -3.0*x*std::sqrt(1.0-x*x);
    if( ell == 2 && k == 2 ) return 3.0*(1.0-x*x);
    if( ell == 3 && k == 0 ) return x/2.0 * (5.0*x*x-3.0);
    if( ell == 3 && k == 1 ) return 1.5*(1.0-5.0*x*x)*std::sqrt(1.0-x*x);
    if( ell == 3 && k == 2 ) return 15.0*x*(1.0-x*x);
    if( ell == 3 && k == 3 ) return -15.0*std::pow(1.0-x*x,1.5);

    return -9999999999.99;
}

//---------------------------------------------------------------------------//

double getclk( unsigned const ell, int const k )
{
    using std::sqrt;
    using std::abs;
    using rtt_sf::factorial;
    using rtt_sf::kronecker_delta;
    return sqrt( (2.0 - kronecker_delta(k,0) ) 
		 * factorial(ell-abs( k)) / (1.0*factorial(ell+abs( k)) ));
}

//---------------------------------------------------------------------------//

void test_quad_services_with_1D_S2_quad( rtt_dsxx::UnitTest & ut )
{   
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
    size_t const sn_ord_ref( 2                   );
    string const qname_ref ( "1D Gauss Legendre" );
    size_t const n_ang_ref ( 2                   );
    
    // Banner
    cout << "\nTesting the "  << qname_ref << "S"
	 << sn_ord_ref << " quadrature set." << endl << endl;
    
    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::GaussLeg, sn_ord_ref ); 
    
    // print the name of the quadrature set that we are testing.
    string const qname   (  spQuad->name()         );
    size_t const snOrder(   spQuad->getSnOrder()   );
    size_t const numOrdinates( spQuad->getNumOrdinates() );
    double const sumwt(     spQuad->getNorm() );
    
    // check basic quadrature setup.
    if( snOrder != sn_ord_ref ) 
    {
	ut.failure("Found incorrect Sn Order.");
    }
    else 
    {
	ut.passes("Found correct Sn Order.");
    }
    if( numOrdinates != n_ang_ref  )
    {
	ut.failure("Found incorrect number of ordinates.");
    }
    else 
    {
	ut.passes("Found correct number of ordinates.");
    }
    if( qname != qname_ref  )
    {
	cout << qname << endl;
	ut.failure("Found incorrect name of quadrature set.");
    }
    else 
    {
	ut.passes("Found correct name of quadrature set.");
    }
    
    // Print a table
    spQuad->display();

    //----------------------------------------
    // Setup QuadServices object

    unsigned const expansionOrder( 1 ); // 2 moments.
    QuadServices qs( spQuad, SN, expansionOrder ); // 
    
    vector<double> const M( qs.getM() );
    unsigned const numMoments( qs.getNumMoments() );

    std::vector< unsigned > dims;
    dims.push_back( numMoments );
    dims.push_back( numOrdinates );
    
    qs.print_matrix( "Mmatrix", M, dims );

    //----------------------------------------
    // For 1D Quadrature we have the following:
    //
    // k == 0, n == l, C_lk == 1, Omega_m == mu_m
    //
    //           2*n+1
    // M_{n,m} = ----- * 1 * Y_n( mu_m )
    //           sumwt
    //
    // Y_n( mu_m ) = P(l=0,k=0)(mu_m)*cos(k*theta)
    //             = P(n,mu_m)
    //----------------------------------------

    std::vector< double > const mu( spQuad->getMu() );
    double const clk(1.0);

    for( size_t n=0; n<numMoments; ++n )
    { 
	double const c( (2.0*n+1.0)/sumwt );

	for( size_t m=0; m<numOrdinates; ++m )
	{
	    if( soft_equiv( M[ n + m*numMoments ], c*clk*P(n,mu[m]) ) )
	    {
		ostringstream msg;
		msg << "M[" << n << "," << m << "] has the expected value." << endl;
		    ut.passes( msg.str() );
	    }
	    else
	    {		
		ostringstream msg;
		msg << "M[" << n << "," << m 
		    << "] does not have the expected value." << endl
		    << "\tFound M[" << n << "," << m << "] = " 
		    << M[ n + m*numMoments ] << ", but was expecting " 
		    << c*clk*P(n,mu[m]) << endl; 
		ut.failure( msg.str() );		
	    }
	}
    } 

    //-----------------------------------//

    vector<double> const D( qs.getD() );
    qs.print_matrix( "Dmatrix", D, dims );

    // The first row of D should contain the quadrature weights.
    {
	unsigned n(0);
	std::vector< double > const wt( spQuad->getWt() );
	for( size_t m=0; m<numOrdinates; ++m )
	{
	    if( soft_equiv( D[ m + n*numOrdinates ], wt[m] ) )
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " matched the expected value." << endl;
		ut.passes( msg.str() );
	    }
	    else
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " did not match the expected value of " 
		    << wt[m] << "." << endl;
		ut.failure( msg.str() );
	    }
	}
    }

    // Ensure D = M^{-1}
    // ------------------------------------------------------------
    {
	if( qs.D_equals_M_inverse() )
	{
	    ut.passes("Found D = inverse(M) for 1D S2.");
	}
	else
	{
	    ut.failure("Oh no! D != inverse(M) for 1D S2.");
	}
    }
    
    // Test applyD function
    // ------------------------------------------------------------
    {
	vector<double> const angularFlux( numOrdinates, 7.0 );
	vector<double> const fluxMoments( qs.applyD( angularFlux ) );

	if( soft_equiv( fluxMoments[0], 14.0 ) &&
	    soft_equiv( fluxMoments[1], 0.0 ) )
	{
	    ut.passes("applyD() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyD() failed to work as expected." << endl
		<< "Expected phi = { 14.0, 0.0} but found phi = { "
		<< fluxMoments[0] << ", " << fluxMoments[1] << " }." << endl;
	    ut.failure(msg.str());
	}
    }

    // Test applyM function
    // ------------------------------------------------------------
    {
	double fm[2] = { 7.0, 0.0 };
	vector<double> const fluxMoments( fm, fm+2 );
	vector<double> const angularFlux( qs.applyM( fluxMoments ) );
	
	if( soft_equiv( angularFlux[0], 3.5 ) &&
	    soft_equiv( angularFlux[1], 3.5 ) )
	{
	    ut.passes("applyM() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyM() failed to work as expected." << endl
		<< "Expected psi = { 3.5, 3.5 } but found psi = { "
		<< angularFlux[0] << ", " << angularFlux[1] << " }." << endl;
	    ut.failure(msg.str());
	}
    }	    


    return;
}

//---------------------------------------------------------------------------//

void test_quad_services_with_1D_S8_quad( rtt_dsxx::UnitTest & ut )
{   
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
    cout << "\nTesting the "  << qname_ref << " S"
	 << sn_ord_ref << " quadrature set." << endl << endl;
    
    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::GaussLeg, sn_ord_ref ); 
    
    // print the name of the quadrature set that we are testing.
    const string qname   (  spQuad->name()         );
    const size_t snOrder(   spQuad->getSnOrder()   );
    const size_t numOrdinates( spQuad->getNumOrdinates() );
    
    // check basic quadrature setup.
    if( snOrder != sn_ord_ref ) 
    {
	ut.failure("Found incorrect Sn Order.");
    }
    else 
    {
	ut.passes("Found correct Sn Order.");
    }
    if( numOrdinates != n_ang_ref  )
    {
	ut.failure("Found incorrect number of ordinates.");
    }
    else 
    {
	ut.passes("Found correct number of ordinates.");
    }
    if( qname != qname_ref  )
    {
	cout << qname << endl;
	ut.failure("Found incorrect name of quadrature set.");
    }
    else 
    {
	ut.passes("Found correct name of quadrature set.");
    }
    
    // Print a table
    spQuad->display();

    //----------------------------------------
    // Setup QuadServices object
    
    QuadServices qs( spQuad, SN, 7 ); // 8 moments
    
    vector<double> const M( qs.getM() );
    unsigned const numMoments( qs.getNumMoments() );

    std::vector< unsigned > dims;
    dims.push_back( numMoments );
    dims.push_back( numOrdinates );
    
    qs.print_matrix( "Mmatrix", M, dims );

    //----------------------------------------
    // For 1D Quadrature we have the following:
    //
    // k == 0, n == l, C_lk == 1, Omega_m == mu_m
    //
    //           2*n+1
    // M_{n,m} = ----- * 1 * Y_n( mu_m )
    //           sumwt
    //
    // Y_n( mu_m ) = P(l=0,k=0)(mu_m)*cos(k*theta)
    //             = P(n,mu_m)
    //----------------------------------------

    std::vector< double > const mu( spQuad->getMu() );
    double const clk(1.0);

    for( size_t n=0; n<numMoments && n<6; ++n )
    { 
	double const c( (2.0*n+1.0)/2.0 );

	for( size_t m=0; m<numOrdinates; ++m )
	{
	    if( soft_equiv( M[ n + m*numMoments ], c*clk*P(n,mu[m]) ) )
	    {
		ostringstream msg;
		msg << "M[" << n << "," << m << "] has the expected value." << endl;
		    ut.passes( msg.str() );
	    }
	    else
	    {		
		ostringstream msg;
		msg << "M[" << n << "," << m 
		    << "] does not have the expected value." << endl
		    << "\tFound M[" << n << "," << m << "] = " 
		    << M[ n + m*numMoments ] << ", but was expecting " 
		    << c*clk*P(n,mu[m]) << endl; 
		ut.failure( msg.str() );		
	    }
	}
    } 

    //-----------------------------------//

    vector<double> const D( qs.getD() );
    qs.print_matrix( "Dmatrix", D, dims );
    
    // The first row of D should contain the quadrature weights.
    {
	unsigned n(0);
	std::vector< double > const wt( spQuad->getWt() );
	for( size_t m=0; m<numOrdinates; ++m )
	{
	    if( soft_equiv( D[ m + n*numOrdinates ], wt[m] ) )
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " matched the expected value." << endl;
		ut.passes( msg.str() );
	    }
	    else
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " did not match the expected value of " 
		    << wt[m] << "." << endl;
		ut.failure( msg.str() );
	    }
	}
    }

    return;
}

//---------------------------------------------------------------------------//

void test_quad_services_with_3D_S2_quad( rtt_dsxx::UnitTest & ut )
{   
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
    
    // we will only look at S2 Set in this test.
    const size_t sn_ord_ref( 2                    );
    const string qname_ref ( "3D Level Symmetric" );
    const size_t n_ang_ref ( 8                   );
    
    // Banner
    cout << "\nTesting the "  << qname_ref << " S"
	 << sn_ord_ref << " quadrature set." << endl << endl;
    
    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, sn_ord_ref ); 
    
    // print the name of the quadrature set that we are testing.
    const string qname   (  spQuad->name()         );
    const size_t snOrder(   spQuad->getSnOrder()   );
    const size_t numOrdinates( spQuad->getNumOrdinates() );
    const double sumwt(     spQuad->getNorm() );

    // check basic quadrature setup.
    if( snOrder != sn_ord_ref ) 
    {
	ut.failure("Found incorrect Sn Order.");
    }
    else 
    {
	ut.passes("Found correct Sn Order.");
    }
    if( numOrdinates != n_ang_ref  )
    {
	ut.failure("Found incorrect number of ordinates.");
    }
    else 
    {
	ut.passes("Found correct number of ordinates.");
    }
    if( qname != qname_ref  )
    {
	cout << qname << endl;
	ut.failure("Found incorrect name of quadrature set.");
    }
    else 
    {
	ut.passes("Found correct name of quadrature set.");
    }
    
    // Print a table
    spQuad->display();

    //----------------------------------------
    // Setup QuadServices object

    QuadServices qs( spQuad );
    
    vector<double> const M( qs.getM() );
    unsigned const numMoments( qs.getNumMoments() );

    {
        std::vector< unsigned > dims;
        dims.push_back( numMoments );
        dims.push_back( numOrdinates );
        qs.print_matrix( "Mmatrix", M, dims );
    }

    //----------------------------------------
    // For 3D Quadrature we have the following:
    //
    // n maps to the index pair (l,k) via qs.n2kl
    // 
    //                       2*l+1
    // M_{n,m} = M_{l,k,m} = ----- * c_{l,k} * Y_{l,k}( mu_m )
    //                       sumwt
    //
    // Y_n( mu_m ) = P(l=0,k=0)(mu_m)*cos(k*theta)
    //             = P(n,mu_m)
    //----------------------------------------

    std::vector< double > const mu( spQuad->getMu() );
    std::vector< double > const eta( spQuad->getEta() );
    std::vector< double > const xi( spQuad->getXi() );

    for( size_t n=0; n<numMoments; ++n )
    { 
	unsigned const ell( qs.lkPair( n ).first  );
	int      const k(   qs.lkPair( n ).second );
	double   const c(   ( 2.0*ell+1.0 ) / sumwt );
	
        if( ell < 4 )
        for( size_t m=0; m<numOrdinates; ++m )
        {
            double expVal = c*getclk(ell,k)*P(ell,std::abs( k),xi[m]);
            double phi    = QuadServices::compute_azimuthalAngle( mu[m], eta[m], xi[m] );
            if( k<0 )
                expVal *= std::sin(-1*k*phi);
            else
                expVal *= std::cos(k*phi);
            
            if( soft_equiv( M[ n + m*numMoments ], expVal ) )
            {
                ostringstream msg;
                msg << "M[" << n << "," << m 
                    << "] has the expected value." << endl;
                ut.passes( msg.str() );
            }
            else
            {		
                ostringstream msg;
                msg << "M[" << n << "," << m 
                    << "] does not have the expected value." << endl
                    << "\tFound M[" << n << "," << m << "] = " 
                    << M[ n + m*numMoments ] << ", but was expecting " 
                    << expVal << "\n"
                    << "\t(l,k) = " << ell << ", " << k << "\n"
                    << endl; 
                ut.failure( msg.str() );		
            }
	}
    } 

    //-----------------------------------//

    vector<double> const D( qs.getD() );
    {
        std::vector< unsigned > dims;
        dims.push_back( numOrdinates );
        dims.push_back( numMoments );
        qs.print_matrix( "Dmatrix", D, dims );
    }
    
    // The first row of D should contain the quadrature weights.
    {
	unsigned n(0);
	std::vector< double > const wt( spQuad->getWt() );
	for( size_t m=0; m<numOrdinates; ++m )
	{
	    if( soft_equiv( D[ m + n*numOrdinates ], wt[m] ) )
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " matched the expected value." << endl;
		ut.passes( msg.str() );
	    }
	    else
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " did not match the expected value of " 
		    << wt[m] << "." << endl;
		ut.failure( msg.str() );
	    }
	}
    }

    // Ensure D = M^{-1}
    // ------------------------------------------------------------
    {
	if( qs.D_equals_M_inverse() )
	{
	    ut.passes("Found D = inverse(M) for 3D S2.");
	}
	else
	{
	    ut.failure("Oh no! D != inverse(M) for 3D S2.");
	}
    }
    
    // Test applyD function
    // ------------------------------------------------------------
    {
        // Isotropic angular flux -> only 1 non-zero moment.
        double magnitude(7.0);
	vector<double> const angularFlux( numOrdinates, magnitude );
	vector<double> const fluxMoments( qs.applyD( angularFlux ) );
        vector<double> expectedPhi( numMoments, 0.0 );
        expectedPhi[0]=magnitude*sumwt;

	if( soft_equiv( fluxMoments.begin(), fluxMoments.end(),
                        expectedPhi.begin(), expectedPhi.end() ) )
	{
	    ut.passes("applyD() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyD() failed to work as expected." << endl
		<< "\tExpected phi = { " << expectedPhi[0] << ", 0.0, ... 0.0 } "
                << "but found \n\tphi = {";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << "\n\t" << fluxMoments[i];
            msg << " }." << endl;
	    ut.failure(msg.str());
	}
    }

    // Test applyM function
    // ------------------------------------------------------------
    {
        // moments that are all zero except first entry are equal to an
        // isotropic angular flux.x
        double magnitude(7.0);
	vector<double> fluxMoments( numMoments, 0.0 );
        fluxMoments[0]=magnitude;
	vector<double> const angularFlux( qs.applyM( fluxMoments ) );
        vector<double> expectedPsi( numOrdinates, magnitude/sumwt );
	
	if( soft_equiv( angularFlux.begin(), angularFlux.end(),
                        expectedPsi.begin(), expectedPsi.end() ) )
	{
	    ut.passes("applyM() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyM() failed to work as expected." << endl
		<< "Expected psi = { ";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << expectedPsi[i] << "\n";
            msg << " }, but found psi = { ";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << angularFlux[i] << "\n"; 
            msg << " }." << endl;
	    ut.failure(msg.str());
	}
    }	    

    // Test applyM and applyD for anisotropic angular flux.
    // ------------------------------------------------------------
    {
        QuadServices qsm( spQuad, GALERKIN );
        double magnitude(7.0);
        for( size_t i=0; i< numOrdinates; ++i )
        {
            vector<double> psi(numOrdinates,0.0);
            psi[i]=magnitude;
            vector<double> phi( qsm.applyD( psi ) );
            vector<double> psi2( qsm.applyM( phi ) );
            if( soft_equiv( psi.begin(),  psi.end(),
                            psi2.begin(), psi2.end() ) )
            {
                ostringstream msg;
                msg << "Recovered psi = M D psi for case i=" << i << endl;
                ut.passes(msg.str());
            }
            else
            {
                ostringstream msg;
                msg << "Failed to recover psi = M D psi for case i=" << i << endl;
                ut.failure(msg.str());
            }
        }
    }
    return;
}

//---------------------------------------------------------------------------//

void test_quad_services_with_3D_S4_quad( rtt_dsxx::UnitTest & ut )
{   
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
    
    // we will only look at S2 Set in this test.
    const size_t sn_ord_ref( 4                    );
    const string qname_ref ( "3D Level Symmetric" );
    const size_t n_ang_ref ( 24                   );
    
    // Banner
    cout << "\nTesting the "  << qname_ref << " S"
	 << sn_ord_ref << " quadrature set." << endl << endl;
    
    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, sn_ord_ref ); 
    
    // print the name of the quadrature set that we are testing.
    const string qname   (  spQuad->name()         );
    const size_t snOrder(   spQuad->getSnOrder()   );
    const size_t numOrdinates( spQuad->getNumOrdinates() );
    const double sumwt(     spQuad->getNorm() );

    // check basic quadrature setup.
    if( snOrder != sn_ord_ref ) 
    {
	ut.failure("Found incorrect Sn Order.");
    }
    else 
    {
	ut.passes("Found correct Sn Order.");
    }
    if( numOrdinates != n_ang_ref  )
    {
	ut.failure("Found incorrect number of ordinates.");
    }
    else 
    {
	ut.passes("Found correct number of ordinates.");
    }
    if( qname != qname_ref  )
    {
	cout << qname << endl;
	ut.failure("Found incorrect name of quadrature set.");
    }
    else 
    {
	ut.passes("Found correct name of quadrature set.");
    }
    
    // Print a table
    spQuad->display();

    //----------------------------------------
    // Setup QuadServices object
    
    QuadServices qs( spQuad );
    
    vector<double> const M( qs.getM() );
    unsigned const numMoments( qs.getNumMoments() );

    {
        std::vector< unsigned > dims;
        dims.push_back( numMoments );
        dims.push_back( numOrdinates );
        qs.print_matrix( "Mmatrix", M, dims );
    }

    //----------------------------------------
    // For 3D Quadrature we have the following:
    //
    // n maps to the index pair (l,k) via qs.n2kl
    // 
    //                       2*l+1
    // M_{n,m} = M_{l,k,m} = ----- * c_{l,k} * Y_{l,k}( mu_m )
    //                       sumwt
    //
    // Y_n( mu_m ) = P(l=0,k=0)(mu_m)*cos(k*theta)
    //             = P(n,mu_m)
    //----------------------------------------

    std::vector< double > const mu( spQuad->getMu() );
    std::vector< double > const eta( spQuad->getEta() );
    std::vector< double > const xi( spQuad->getXi() );

    for( size_t n=0; n<numMoments; ++n )
    { 
	unsigned const ell( qs.lkPair( n ).first  );
	int      const k(   qs.lkPair( n ).second );
	double   const c(   ( 2.0*ell+1.0 ) / sumwt );
	
        if( ell < 4 )
        for( size_t m=0; m<numOrdinates; ++m )
        {
            double expVal = c*getclk(ell,k)*P(ell,std::abs( k),xi[m]);
            double phi    = QuadServices::compute_azimuthalAngle( mu[m], eta[m], xi[m] );
            if( k<0 )
            {
                expVal *= std::sin(-1.0*k*phi) ;
            }
            else
                expVal *= std::cos(k*phi) ;
            
            if( soft_equiv( M[ n + m*numMoments ], expVal ) )
            {
                ostringstream msg;
                msg << "M[" << n << "," << m 
                    << "] has the expected value." << endl;
                ut.passes( msg.str() );
            }
            else
            {		
                ostringstream msg;
                msg << "M[" << n << "," << m 
                    << "] does not have the expected value." << endl
                    << "\tFound M[" << n << "," << m << "] = " 
                    << M[ n + m*numMoments ] << ", but was expecting " 
                    << expVal << "\n"
                    << "\t(l,k) = " << ell << ", " << k << "\n"
                    << "\tOmega = " << mu[m] << ", " << eta[m] << ", "
                    << xi[m] << "\n" << endl; 
                ut.failure( msg.str() );		
            }
	}
    } 

    //-----------------------------------//

    vector<double> const D( qs.getD() );
    {
        std::vector< unsigned > dims;
        dims.push_back( numOrdinates );
        dims.push_back( numMoments );
        qs.print_matrix( "Dmatrix", D, dims );
    }
    
    // The first row of D should contain the quadrature weights.
    {
	unsigned n(0);
	std::vector< double > const wt( spQuad->getWt() );
	for( size_t m=0; m<numOrdinates; ++m )
	{
	    if( soft_equiv( D[ m + n*numOrdinates ], wt[m] ) )
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " matched the expected value." << endl;
		ut.passes( msg.str() );
	    }
	    else
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " did not match the expected value of " 
		    << wt[m] << "." << endl;
		ut.failure( msg.str() );
	    }
	}
    }

    // Ensure D = M^{-1}
    // ------------------------------------------------------------
    {
	if( qs.D_equals_M_inverse() )
	{
	    ut.passes("Found D = inverse(M) for 3D S4.");
	}
	else
	{
	    ut.failure("Oh no! D != inverse(M) for 3D S4.");
	}
    }

    // Test applyD function
    // ------------------------------------------------------------
    {
        // Isotropic angular flux -> only 1 non-zero moment.
        double magnitude(7.0);
	vector<double> const angularFlux( numOrdinates, magnitude );
	vector<double> const fluxMoments( qs.applyD( angularFlux ) );
        vector<double> expectedPhi( numMoments, 0.0 );
        expectedPhi[0]=magnitude*sumwt;

	if( soft_equiv( fluxMoments.begin(), fluxMoments.end(),
                        expectedPhi.begin(), expectedPhi.end() ) )
	{
	    ut.passes("applyD() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyD() failed to work as expected." << endl
		<< "Expected phi = { " << expectedPhi[0] << ", 0.0, ... 0.0 } "
                << "but found phi = { \n";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << fluxMoments[i] << "\n";
            msg << " }." << endl;
	    ut.failure(msg.str());
	}
    }

    // Test applyM function
    // ------------------------------------------------------------
    {
        // moments that are all zero except first entry are equal to an
        // isotropic angular flux.x
        double magnitude(7.0);
	vector<double> fluxMoments( numMoments, 0.0 );
        fluxMoments[0]=magnitude;
	vector<double> const angularFlux( qs.applyM( fluxMoments ) );
        vector<double> expectedPsi( numOrdinates, magnitude/sumwt );
	
	if( soft_equiv( angularFlux.begin(), angularFlux.end(),
                        expectedPsi.begin(), expectedPsi.end() ) )
	{
	    ut.passes("applyM() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyM() failed to work as expected." << endl
		<< "Expected psi = { ";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << expectedPsi[i] << "\n";
            msg << " }, but found psi = { ";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << angularFlux[i] << "\n"; 
            msg << " }." << endl;
	    ut.failure(msg.str());
	}
    }	    

    // Test applyM and applyD for anisotropic angular flux.
    // ------------------------------------------------------------
    {
        QuadServices qsm( spQuad, GALERKIN );
        double magnitude(7.0);
        for( size_t i=0; i< numOrdinates; ++i )
        {
            vector<double> psi(numOrdinates,0.0);
            psi[i]=magnitude;
            vector<double> phi( qsm.applyD( psi ) );
            vector<double> psi2( qsm.applyM( phi ) );
            if( soft_equiv( psi.begin(),  psi.end(),
                            psi2.begin(), psi2.end() ) )
            {
                ostringstream msg;
                msg << "Recovered psi = M D psi for case i=" << i << endl;
                ut.passes(msg.str());
            }
            else
            {
                ostringstream msg;
                msg << "Failed to recover psi = M D psi for case i=" << i << endl;
                ut.failure(msg.str());
            }
        }
    }
    
    return;
}

//---------------------------------------------------------------------------//

void test_quad_services_with_2D_S6_quad( rtt_dsxx::UnitTest & ut )
{   
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
    const size_t sn_ord_ref( 6                    );
    const string qname_ref ( "2D Level Symmetric" );
    const size_t n_ang_ref ( 24                    );
    
    // Banner
    cout << "\nTesting the "  << qname_ref << " S"
	 << sn_ord_ref << " quadrature set." << endl << endl;
    
    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym2D, sn_ord_ref ); 
    
    // print the name of the quadrature set that we are testing.
    const string qname   (  spQuad->name()         );
    const size_t snOrder(   spQuad->getSnOrder()   );
    const size_t numOrdinates( spQuad->getNumOrdinates() );
    const double sumwt(     spQuad->getNorm() );

    // check basic quadrature setup.
    if( snOrder != sn_ord_ref ) 
    {
	ut.failure("Found incorrect Sn Order.");
    }
    else 
    {
	ut.passes("Found correct Sn Order.");
    }
    if( numOrdinates != n_ang_ref  )
    {
	ut.failure("Found incorrect number of ordinates.");
    }
    else 
    {
	ut.passes("Found correct number of ordinates.");
    }
    if( qname != qname_ref  )
    {
	cout << qname << endl;
	ut.failure("Found incorrect name of quadrature set.");
    }
    else 
    {
	ut.passes("Found correct name of quadrature set.");
    }
    
    // Print a table
    spQuad->display();

    //----------------------------------------
    // Setup QuadServices object
       
    QuadServices qs( spQuad );
    
    vector<double> const M( qs.getM() );
    unsigned const numMoments( qs.getNumMoments() );

    {
        std::vector< unsigned > dims;
        dims.push_back( numMoments );
        dims.push_back( numOrdinates );
        qs.print_matrix( "Mmatrix", M, dims );
    }

    //----------------------------------------
    // For 3D Quadrature we have the following:
    //
    // n maps to the index pair (l,k) via qs.n2kl
    // 
    //                       2*l+1
    // M_{n,m} = M_{l,k,m} = ----- * c_{l,k} * Y_{l,k}( mu_m )
    //                       sumwt
    //
    // Y_n( mu_m ) = P(l=0,k=0)(mu_m)*cos(k*theta)
    //             = P(n,mu_m)
    //----------------------------------------

    std::vector< double > const mu( spQuad->getMu() );

    for( size_t n=0; n<numMoments; ++n )
    { 
	unsigned const ell( qs.lkPair( n ).first  );
	int      const k(   qs.lkPair( n ).second );
	double   const c(   ( 2.0*ell+1.0 ) / sumwt );
        
	if( k == 0 ) 
	{
	    for( size_t m=0; m<numOrdinates; ++m )
	    {
                double expectedValue = c*getclk(ell,k)*P( ell, mu[m] );
		if( soft_equiv( M[ n + m*numMoments ], expectedValue ) )
		{
		    ostringstream msg;
		    msg << "M[" << n << "," << m 
			<< "] has the expected value." << endl;
		    ut.passes( msg.str() );
		}
		else
		{		
		    ostringstream msg;
		    msg << "M[" << n << "," << m 
			<< "] does not have the expected value." << endl
			<< "\tFound M[" << n << "," << m << "] = " 
			<< M[ n + m*numMoments ] << ", but was expecting " 
			<< c*getclk(ell,k)*P(ell,mu[m]) << endl; 
		    ut.failure( msg.str() );		
		}
	    }
	}
    } 

    //-----------------------------------//

    vector<double> const D( qs.getD() );

    {
        std::vector< unsigned > dims;
        dims.push_back( numOrdinates );
        dims.push_back( numMoments );
        qs.print_matrix( "Dmatrix", D, dims );
    }
    
    // The first row of D should contain the quadrature weights.
    {
	unsigned n(0);
	std::vector< double > const wt( spQuad->getWt() );
	for( size_t m=0; m<numOrdinates; ++m )
	{
	    if( soft_equiv( D[ m + n*numOrdinates ], wt[m] ) )
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " matched the expected value." << endl;
		ut.passes( msg.str() );
	    }
	    else
	    {
		ostringstream msg;
		msg << "D[" << m << "," << n << "] = " 
		    << D[ m + n*numOrdinates ] 
		    << " did not match the expected value of " 
		    << wt[m] << "." << endl;
		ut.failure( msg.str() );
	    }
	}
    }

    // Ensure D = M^{-1}
    // ------------------------------------------------------------
    {
	if( qs.D_equals_M_inverse() )
	{
	    ut.passes("Found D = inverse(M) for 3D S4.");
	}
	else
	{
	    ut.failure("Oh no! D != inverse(M) for 3D S4.");
	}
    }

    // Test applyD function
    // ------------------------------------------------------------
    {
        // Isotropic angular flux -> only 1 non-zero moment.
        double magnitude(7.0);
	vector<double> const angularFlux( numOrdinates, magnitude );
	vector<double> const fluxMoments( qs.applyD( angularFlux ) );
        vector<double> expectedPhi( numMoments, 0.0 );
        expectedPhi[0]=magnitude*sumwt;

	if( soft_equiv( fluxMoments.begin(), fluxMoments.end(),
                        expectedPhi.begin(), expectedPhi.end() ) )
	{
	    ut.passes("applyD() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyD() failed to work as expected." << endl
		<< "Expected phi = { " << expectedPhi[0] << ", 0.0, ... 0.0 } "
                << "but found phi = { \n";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << fluxMoments[i] << "\n";
            msg << " }." << endl;
	    ut.failure(msg.str());
	}
    }

    // Test applyM function
    // ------------------------------------------------------------
    {
        // moments that are all zero except first entry are equal to an
        // isotropic angular flux.x
        double magnitude(7.0);
	vector<double> fluxMoments( numMoments, 0.0 );
        fluxMoments[0]=magnitude;
	vector<double> const angularFlux( qs.applyM( fluxMoments ) );
        vector<double> expectedPsi( numOrdinates, magnitude/sumwt );
	
	if( soft_equiv( angularFlux.begin(), angularFlux.end(),
                        expectedPsi.begin(), expectedPsi.end() ) )
	{
	    ut.passes("applyM() appears to work.");
	}
	else
	{
	    ostringstream msg;
	    msg << "applyM() failed to work as expected." << endl
		<< "Expected psi = { ";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << expectedPsi[i] << "\n";
            msg << " }, but found psi = { ";
            for( size_t i=0; i< numOrdinates; ++i)
                msg << angularFlux[i] << "\n"; 
            msg << " }." << endl;
	    ut.failure(msg.str());
	}
    }	    

    // Test applyM and applyD for anisotropic angular flux.
    // ------------------------------------------------------------
    {
        QuadServices qsm( spQuad, GALERKIN );
        double magnitude(7.0);
        for( size_t i=0; i< numOrdinates; ++i )
        {
            vector<double> psi(numOrdinates,0.0);
            psi[i]=magnitude;
            vector<double> phi( qsm.applyD( psi ) );
            vector<double> psi2( qsm.applyM( phi ) );
            if( soft_equiv( psi.begin(),  psi.end(),
                            psi2.begin(), psi2.end() ) )
            {
                ostringstream msg;
                msg << "Recovered psi = M D psi for case i=" << i << endl;
                ut.passes(msg.str());
            }
            else
            {
                ostringstream msg;
                msg << "Failed to recover psi = M D psi for case i=" << i << endl;
                ut.failure(msg.str());
            }
        }
    }
    
    return;
}

//---------------------------------------------------------------------------//

void test_quad_services_alt_constructor( rtt_dsxx::UnitTest & ut )
{
    using rtt_dsxx::SP;
    using rtt_dsxx::soft_equiv;

    using std::endl;
    using std::vector;
    using std::ostringstream;

    typedef std::pair< unsigned, int > lk_index;

    //----------------------------------------
    // Setup Quadrature set
    
    // we will only look at S2 Sets in this test.
    size_t const snOrder( 2 );

    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, snOrder ); 
    
    // Create a vector that designates the (l,k) moments that will be used
    unsigned const numMoments( spQuad->getNumOrdinates() );
    unsigned n(0);
    vector< lk_index > lkMoments;
    
    // Copy algorithm from compute_n2lk_3D()
    // -------------------------------------
    // Choose: l= 0, ..., N-1, k = -l, ..., l
    for( unsigned ell=0; ell< snOrder; ++ell )
	for( int k(-1*static_cast<int>(ell));
             std::abs( k ) <= static_cast<int>(ell); ++k, ++n )
	    lkMoments.push_back( lk_index(ell,k) );

    // Add ell=N and k<0
    {
	unsigned ell( snOrder );
	for( int k(-1*static_cast<int>(ell)); k<0; ++k, ++n )
	    lkMoments.push_back( lk_index(ell,k) );
    }

    // Add ell=N, k>0, k odd
    {
	unsigned ell( snOrder );
	for( int k=1; k<=static_cast<int>(ell); k+=2, ++n )
	    lkMoments.push_back( lk_index(ell,k) );
    }

    // Add ell=N+1 and k<0, k even
    {
	unsigned ell( snOrder+1 );
	for( int k(-1*static_cast<int>(ell)+1); k<0; k+=2, ++n )
	    lkMoments.push_back( lk_index(ell,k) );
    }

    //----------------------------------------
    // Setup QuadServices object using alternate constructor.
    
    QuadServices qsStd( spQuad, GALERKIN );
    QuadServices qsAlt( spQuad, lkMoments, GALERKIN );

    for( unsigned n=0; n<numMoments; ++n )
    {
	lk_index stdIndexValues( qsStd.lkPair(n) );
	lk_index altIndexValues( qsAlt.lkPair(n) );
	if( stdIndexValues.first == altIndexValues.first &&
	    stdIndexValues.second == altIndexValues.second )
	{
	    ostringstream msg;
	    msg << "Alternate Constructor -- lk_index has expected value for moment "
		<< n << "." << endl;
	    ut.passes(msg.str());
	}
	else
	{
	    ostringstream msg;
	    msg << "Alternate Constructor -- "
		<< "lk_index does not have the expected value for moment "
		<< n << "." << endl
		<< "Found lk_index = (" << altIndexValues.first << ", "
		<< altIndexValues.second << ") but expected (" << stdIndexValues.first
		<< ", " << stdIndexValues.second << ")." << endl;
	    ut.passes(msg.str());
	}
    }
    return;
}

//---------------------------------------------------------------------------//

void test_quad_services_SVD( rtt_dsxx::UnitTest & ut )
{
    using rtt_dsxx::SP;
    using rtt_dsxx::soft_equiv;

    using std::endl;
    using std::vector;
    using std::ostringstream;

    typedef std::pair< unsigned, int > lk_index;

    //----------------------------------------
    // Setup Quadrature set
    
    // we will only look at S2 Sets in this test.
    size_t const snOrder( 2 );

    // Create a quadrature set from a temporary instance of a
    // QuadratureCreator factory object.
    SP< const Quadrature > spQuad;
    spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, snOrder ); 
    
    // Create a vector that designates the (l,k) moments that will be used
    unsigned n(0);
    vector< lk_index > lkMoments;
    
    // Copy algorithm from compute_n2lk_3D()
    // -------------------------------------
    // Choose: l= 0, ..., N-1, k = -l, ..., l
    for( unsigned ell=0; ell< snOrder; ++ell )
	for( int k(-1*static_cast<int>(ell));
             std::abs( k ) <= static_cast<int>(ell); ++k, ++n )
	    lkMoments.push_back( lk_index(ell,k) );

    // Add ell=N and k<0
    {
	unsigned ell( snOrder );
	for( int k(-1*static_cast<int>(ell)); k<0; ++k, ++n )
	    lkMoments.push_back( lk_index(ell,k) );
    }

    // Add ell=N, k>0, k odd
    {
	unsigned ell( snOrder );
	for( int k=1; k<=static_cast<int>(ell); k+=2, ++n )
	    lkMoments.push_back( lk_index(ell,k) );
    }

    // Add ell=N+1 and k<0, k even
    {
	unsigned ell( snOrder+1 );
	for( int k(-1*static_cast<int>(ell)+1); k<0; k+=2, ++n )
	    lkMoments.push_back( lk_index(ell,k) );
    }

    //----------------------------------------
    // Setup QuadServices object using alternate constructor.
    
    QuadServices qsStd( spQuad, SVD );
    QuadServices qsAlt( spQuad, lkMoments, SVD );

    unsigned const numMoments( qsStd.getNumMoments() );

    for( unsigned n=0; n<numMoments; ++n )
    {
	lk_index stdIndexValues( qsStd.lkPair(n) );
	lk_index altIndexValues( qsAlt.lkPair(n) );
	if( stdIndexValues.first == altIndexValues.first &&
	    stdIndexValues.second == altIndexValues.second )
	{
	    ostringstream msg;
	    msg << "Alternate Constructor -- lk_index has expected value for moment "
		<< n << "." << endl;
	    ut.passes(msg.str());
	}
	else
	{
	    ostringstream msg;
	    msg << "Alternate Constructor -- "
		<< "lk_index does not have the expected value for moment "
		<< n << "." << endl
		<< "Found lk_index = (" << altIndexValues.first << ", "
		<< altIndexValues.second << ") but expected (" << stdIndexValues.first
		<< ", " << stdIndexValues.second << ")." << endl;
	    ut.passes(msg.str());
	}
    }
    return;
}

//---------------------------------------------------------------------------//

void test_dnz( UnitTest & ut )
{
    vector<double> vec(9,1.0);
    if( QuadServices::diagonal_not_zero(vec,3,3) )
        ut.passes("Found non-zero diagonal.");
    else
        ut.failure("Incorrectly reported zero on the diagonal.");

    vec[0]=0.0; // Put a zero on the diagonal.
    if( QuadServices::diagonal_not_zero(vec,3,3) )
        ut.failure("Incorrectly reported no zeroes on the diagonal.");
    else
        ut.passes("Correctly found a zero on the diagonal.");
    
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

        test_dnz(ut);
  	test_quad_services_with_1D_S2_quad(ut);
   	test_quad_services_with_1D_S8_quad(ut);
 	test_quad_services_with_3D_S2_quad(ut);
  	test_quad_services_with_3D_S4_quad(ut);
  	test_quad_services_with_2D_S6_quad(ut);
   	test_quad_services_alt_constructor(ut);
   	test_quad_services_SVD(ut);
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
