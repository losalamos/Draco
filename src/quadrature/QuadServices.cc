//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices.cc
 * \author Kelly Thompson
 * \date   Mon Nov  8 11:17:12 2004
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>
#include <cmath>
#include <sstream>

// Vendor software
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>

// Draco software
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

#include "QuadServices.hh"

#include <iostream>

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor builds square D and M operators using Morel's
 * Galerkin-Sn heuristic. 
 * \param spQuad_ a smart pointer to a Quadrature object.
 * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
 */
QuadServices::QuadServices( rtt_dsxx::SP< const Quadrature > const spQuad_ )
    : spQuad(     spQuad_ ), 
      numMoments( spQuad->getNumAngles() ),
      n2lk(       compute_n2lk() ),
      Mmatrix(    computeM() ),
      Dmatrix(    computeD() )	
{ 
    using rtt_dsxx::soft_equiv;
    using rtt_units::PI;

    Ensure( factorial(0) == 1 );
    Ensure( factorial(1) == 1 );
    Ensure( factorial(2) == 2 );
    Ensure( factorial(3) == 6 );
    Ensure( factorial(4) == 24 );

    Ensure( kronecker_delta(0,0) == 1 );
    Ensure( kronecker_delta(1,0) == 0 );
    Ensure( kronecker_delta(1,1) == 1 );
    Ensure( kronecker_delta(0,1) == 0 );
    Ensure( kronecker_delta(-1,0) == 0 );

    Ensure( soft_equiv( compute_clk( 0, 0 ), 1.0 ));
    Ensure( soft_equiv( compute_clk( 1, -1 ), 1.0 ));
    Ensure( soft_equiv( compute_clk( 1, 0 ), 1.0 ));
    Ensure( soft_equiv( compute_clk( 1, 1 ), 1.0 ));

    double const mu1( std::sqrt(2.0)/2.0 );
    double const mu2( std::sqrt(3.0)/3.0 );
    Ensure( soft_equiv( compute_azimuthalAngle( 1.0, 0.0, 0.0 ), 0.0 ) );
    // Ensure( compute_azimuthalAngle( mu1, mu1, 0.0 ) == 0.0 );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, mu2, mu2 ), PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, -1.0*mu2, mu2 ), 3.0*PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, -1.0*mu2, -1.0*mu2 ), 5.0*PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, mu2, -1.0*mu2 ), 7.0*PI/4.0 )  );

    Check( soft_equiv(gsl_sf_legendre_Plm( 0, 0, 0.5 ), 1.0 ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 1, 0, 0.5 ), 0.5 ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 1, 1, mu2 ), -1.0*std::sqrt(1.0-mu2*mu2) ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 2, 2, mu2 ), 3.0*(1.0-mu2*mu2) ));

    Check( n2lk[0].first == 0 );
    Check( n2lk[1].first == 1 );
    Check( n2lk[2].first == 1 );
    Check( n2lk[3].first == 1 );
    Check( n2lk[4].first == 2 );
    Check( n2lk[5].first == 2 );
    Check( n2lk[6].first == 2 );
    Check( n2lk[7].first == 3 );

    Check( n2lk[0].second == 0 );
    Check( n2lk[1].second == -1 );
    Check( n2lk[2].second == 0 );
    Check( n2lk[3].second == 1 );
    Check( n2lk[4].second == -2 );
    Check( n2lk[5].second == -1 );
    Check( n2lk[6].second == 1 );
    Check( n2lk[7].second == -2 );

//---------------------------------------------------------------------------//

    std::vector< unsigned > dims;
    dims.push_back( spQuad->getNumAngles() );
    dims.push_back( numMoments );
    print_matrix( "Mmatrix", Mmatrix, dims );

    Ensure( D_equals_M_inverse() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor that allows the user to pick the (k,l) moments to use.
 * \param spQuad_     a smart pointer to a Quadrature object.
 * \param lkMoments_  vector of tuples that maps from index n to (k,l).
 * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
 */
QuadServices::QuadServices( 
    rtt_dsxx::SP< const Quadrature > const   spQuad_,
    std::vector< lk_index >          const & lkMoments_ )
    : spQuad(     spQuad_ ), 
      numMoments( lkMoments_.size() ),
      n2lk(       lkMoments_ ),
      Mmatrix(    computeM() ),
      Dmatrix(    computeD() )	
{ 
    Ensure( D_equals_M_inverse() );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Apply the action of \f$ M \f$ to the moment-based solution vector
 * \f$ \mathbf\Phi \f$.  That is, \f$ \mathbf\Psi = \mathbf{M}\mathbf\Phi \f$
 * 
 * \param phi Moment based solution vector, \f$ \mathbf\Phi \f$.
 * \return The discrete angular flux, \f$ \mathbf\Psi \f$.
 */
std::vector< double > QuadServices::applyM( std::vector< double > const & phi ) const
{
    Require( phi.size() == numMoments );

    size_t const numAngles( spQuad->getNumAngles() );
    std::vector< double > psi( numAngles, 0.0 );

    for( size_t m=0; m<numAngles; ++m )
	for( size_t n=0; n<numMoments; ++n )
	    psi[ m ] += Mmatrix[ n + m*numMoments ] * phi[n];
    
    return psi;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Apply the action of \f$ D \f$ to the discrete-angle based solution vector
 * \f$ \mathbf\Psi \f$. That is, \f$ \mathbf\Phi = \mathbf{D}\mathbf\Psi \f$.
 * 
 * \param psi Discrete angle-based solution vector, \f$ \mathbf\Psi \f$.
 * \return The moment-based solution vector, \f$ \mathbf\Phi \f$.
 */
std::vector< double > QuadServices::applyD( std::vector< double > const & psi ) const
{
    size_t const numAngles( spQuad->getNumAngles() );
    Require( psi.size() == numAngles );

    std::vector< double > phi( numMoments, 0.0 );

    for( size_t m=0; m<numAngles; ++m )
	for( size_t n=0; n<numMoments; ++n )
	    phi[ n ] += Dmatrix[ m + n*numAngles ] * psi[m];
    
    return phi;
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the discrete-to-moment matrix. 
 *
 * Computes \f$ \mathbf{D} \equiv \mathbf{M}^{-1} \f$.  This private function
 * is called by the constuctor. 
 */
std::vector< double > QuadServices::computeD() const
{
    int n( numMoments );
    int m( spQuad->getNumAngles() );

    // create a copy of Mmatrix to use a temp space.
    std::vector< double > M( Mmatrix );
    std::vector< double > Dmatrix( n*m );

    // Create GSL matrix views of our M and D matrices.
    // LU will get a copy of M.  This matrix will be decomposed into LU. 
    gsl_matrix_view LU = gsl_matrix_view_array( &M[0],       m, n );
    gsl_matrix_view D  = gsl_matrix_view_array( &Dmatrix[0], m, n );

    // Create some local space for the permutation matrix.
    gsl_permutation *p = gsl_permutation_alloc( m );

    // Store information aobut sign changes in this variable.
    int signum(0);

    // Factorize the square matrix M into the LU decomposition PM = LU.  On
    // output the diagonal and upper triangular part of the input matrix M
    // contain the matrix U.  The lower triangular part of the input matrix
    // (excluding the diagonal) contains L. The diagonal elements of L are
    // unity, and are not stored.
    //
    // The permutation matrix P is encoded in the permutation p.  The j-th
    // column of the matrix P is given by the k-th column of the identity,
    // where k=p[j] thej-th element of the permutation vector.  The sign of
    // the permutation is given by signum.  It has the value \f$ (-1)^n \f$,
    // where n is the number of interchanges in the permutation.
    //
    // The algorithm used in the decomposition is Gaussian Elimination with
    // partial pivoting (Golub & Van Loan, Matrix Computations, Algorithm
    // 3.4.1).
    
    int result = gsl_linalg_LU_decomp( &LU.matrix, p, &signum );
    Check( result == 0 );

    // Compute the inverse of the matrix LU from its LU decomposition (LU,p),
    // storing the results in the matrix Dmatrix.  The inverse is computed by
    // solving the system (LU) x = b for each column of the identity matrix.

    result = gsl_linalg_LU_invert( &LU.matrix, p, &D.matrix );
    Check( result == 0 );

    return Dmatrix;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Create the M array (moment-to-discrete matrix).
 * \return The moment-to-discrete matrix.
 *
 * This private member function is called by the constructor. 
 * 
 * The moment-to-discrete matrix will be num_moments by num_angles in size.
 * If the default constructor is used num_moments == num_angles. 
 */
std::vector< double > QuadServices::computeM() const
{
    unsigned const numAngles( spQuad->getNumAngles() );
    unsigned const dim(       spQuad->dimensionality() );
    double sumwt( 0.0 );
    for( size_t m=0; m<numAngles; ++m )
	sumwt += spQuad->getWt(m);

    // resize the M matrix.
    std::vector< double > Mmatrix( numAngles*numMoments, -9999.0 );

    for( unsigned n=0; n<numMoments; ++n )
    {
	unsigned const ell ( n2lk[n].first  );
	int      const k   ( n2lk[n].second ); 
	double   const norm( (2*ell+1)/sumwt );
	double   const clk ( compute_clk(ell,k) );
	
	// Loop over all angles that use these values.
	for( unsigned m=0; m<numAngles; ++m )
	    Mmatrix[ n + m*numMoments ] 
		= norm * clk * spherical_harmonic(m,ell,k);
    }
    return Mmatrix;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the \f$ (\ell,k) \f$ spherical harmonic evaluated at \f$
 * \Omega_m \f$.
 * \param m The index for the current discrete angle.
 * \param ell The \f$ \ell \f$ index for the current spherical harmonic function.
 * \param k The k index for the current spherical harmonic function. 
 * \return The \f$ (\ell,k) \f$ sphercial harmonic evaluated at \f$ \Omega_m. \f$
 *
 * \sa <a
 * href="http://mathworld.wolfram.com/SphericalHarmonic.html">Mathworld's
 * entry for Spherical Harmonic</a>.
 */
double QuadServices::spherical_harmonic( unsigned const m, 
					 unsigned const ell,
					 int      const k   ) const
{
    Require( std::abs(k) <= ell );
    Require( m           <  spQuad->getNumAngles() );

    unsigned const dim( spQuad->dimensionality() );
    double   const mu ( spQuad->getMu(m) );
    double   const eta( dim>1 ? spQuad->getEta(m) : 0.0 );
    double   const xi ( dim>2 ? spQuad->getXi( m) : 0.0 );

    // Compute the azimuthal angle.
    double const azimuthalAngle( compute_azimuthalAngle( mu, eta, xi ) );
    
    double sphHarm(0.0);
    if( k < 0 )
	sphHarm = gsl_sf_legendre_Plm( ell, std::abs(k), mu )
	    * sin( std::abs(k) * azimuthalAngle );
    else
	sphHarm = gsl_sf_legendre_Plm( ell, k, mu ) 
	    * cos( ell * azimuthalAngle );
    
    return sphHarm;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the Azimuthal angle for the current quadrature direction.
 */
double QuadServices::compute_azimuthalAngle( double const mu,
					     double const eta,
					     double const xi ) const
{
    using rtt_units::PI;
    using rtt_dsxx::soft_equiv;

    Require( std::abs(mu)  <= 1.0 );
    Require( std::abs(eta) <= 1.0 );
    Require( std::abs(xi)  <= 1.0 );

    // For 1D sets, we define this angle to be zero.
    if( soft_equiv( eta, 0.0 ) ) return 0.0;

    // For 2D sets, reconstruct xi from known information: 
    // xi*xi = 1.0 - eta*eta - mu*mu
    // Always use positive value for xi.
    double local_xi( xi );
    if( soft_equiv( local_xi,  0.0 ) )
	local_xi = std::sqrt( 1.0 - mu*mu - eta*eta );

    double azimuthalAngle(999.0);

    if( local_xi > 0.0 )
    {
	if( eta > 0.0 )
	    azimuthalAngle = atan(xi/eta);
	else
	    azimuthalAngle = PI - atan(xi/std::abs(eta));
    } 
    else 
    {
	if( eta > 0 )
	    azimuthalAngle = 2*PI - atan(std::abs(xi)/eta);
	else
	    azimuthalAngle = PI + atan(xi/eta);
    }

    // ensure that theta is in the range 0...2*PI.
    Ensure( azimuthalAngle >= 0 );
    Ensure( azimuthalAngle <= 2*PI );
    
    return azimuthalAngle;
}


//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the c(l,k) spherical harmonics coefficient.
 */
double QuadServices::compute_clk( unsigned const ell, int const k ) const
{
    return std::sqrt( ( 2 - kronecker_delta(k,0) ) 
		      * factorial( ell - std::abs(k) )
		      / ( 1.0 * factorial( ell + std::abs(k) ) ) );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Multiply M and D and compare the result to the identity matrix.
 * \return true if M = D^(-1), otherwise false.
 */
bool QuadServices::D_equals_M_inverse() const
{
    using rtt_dsxx::soft_equiv;
    using std::cout;
    using std::endl;

    int n( numMoments );
    int m( spQuad->getNumAngles() );
    int nm( std::min(n,m) );

    // create non-const versions of M and D.
    std::vector< double > Marray( Mmatrix );
    std::vector< double > Darray( Dmatrix );
    std::vector< double > Iarray( nm*nm, -999.0 );
    gsl_matrix_view M = gsl_matrix_view_array( &Marray[0], m, n );
    gsl_matrix_view D = gsl_matrix_view_array( &Darray[0], n, m );
    gsl_matrix_view I = gsl_matrix_view_array( &Iarray[0], nm, nm );
    
    // Compute the matrix-matrix product and sum:
    //
    // I = alpha * op1(M) * op2(D) + beta*I
    //
    // where op1 is one of:
    //    CblasNoTrans    <-->    Use M as provided.
    //    CblasTrans      <-->    Transpose M before multiplication.
    //    CblasConjTRans  <-->    Hermitian transpose M before mult.
    CBLAS_TRANSPOSE_t op( CblasNoTrans );
    double alpha(1.0);
    double beta( 0.0);
    
    gsl_blas_dgemm( op, op, alpha, &M.matrix, &D.matrix, beta, &I.matrix );

    for( unsigned i=0; i<nm; ++i )
	// cout << Iarray[ i ] << endl;
	if( ! soft_equiv( Iarray[ i + i*nm ], 1.0 ) ) return false;

    return true;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 *
 * This function computes the mapping as specified by Morel in "A Hybrid
 * Collocation-Galerkin-Sn Method for Solving the Boltzmann Transport
 * Equation." 
 */
std::vector< QuadServices::lk_index > QuadServices::compute_n2lk() const
{
    unsigned const dim( spQuad->dimensionality() );

    if ( dim == 3 ) return compute_n2lk_3D();
    if ( dim == 2 ) return compute_n2lk_2D();
    Check( dim == 1 );
    return compute_n2lk_1D();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::compute_n2lk_3D() const
{
    unsigned const numAngles( spQuad->getNumAngles() );
    unsigned const snOrder(   spQuad->getSnOrder()   );
    unsigned n(0);
    
    std::vector< lk_index > result;

    // Choose: l= 0, ..., N-1, k = -l, ..., l
    for( unsigned ell=0; ell< snOrder; ++ell )
	for( int k(-1*ell); std::fabs(k) <= ell; ++k, ++n )
	    result.push_back( lk_index(ell,k) );

    // Add ell=N and k<0
    {
	unsigned ell( snOrder );
	for( int k(-1*ell); k<0; ++k, ++n )
	    result.push_back( lk_index(ell,k) );
    }

    // Add ell=N, k>0, k odd
    {
	unsigned ell( snOrder );
	for( int k=1; k<=ell; k+=2, ++n )
	    result.push_back( lk_index(ell,k) );
    }

    // Add ell=N+1 and k<0, k even
    {
	unsigned ell( snOrder+1 );
	for( int k(-1*ell+1); k<0; k+=2, ++n )
	    result.push_back( lk_index(ell,k) );
    }

    Ensure( n == numMoments );
    Ensure( result.size() == numMoments );
    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::compute_n2lk_2D() const
{
    unsigned const numAngles( spQuad->getNumAngles() );
    unsigned const snOrder(   spQuad->getSnOrder()   );
    unsigned n(0);
    
    std::vector< lk_index > result;
    
    // Choose: l= 0, ..., N-1, k = 0, ..., l
    for( unsigned ell=0; ell< snOrder; ++ell )
	for( int k=0; k<= ell; ++k, ++n )
	    result.push_back( lk_index(ell,k) );

    // Add ell=N and k>0, k odd
    {
	unsigned ell( snOrder );
	for( int k=1; k<=ell; k+=2, ++n )
	    result.push_back( lk_index(ell,k) );
    }

    Ensure( n == numMoments );
    Ensure( result.size() == numMoments );
    return result;
}


//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::compute_n2lk_1D() const
{
    unsigned const numAngles( spQuad->getNumAngles() );
    unsigned const snOrder(   spQuad->getSnOrder()   );
    unsigned n(0);
    
    std::vector< lk_index > result;
    
    // Choose: l= 0, ..., N-1, k = 0
    int k(0); // k is always zero for 1D.
    for( unsigned ell=0; ell<snOrder; ++ell, ++n )
	result.push_back( lk_index(ell,k) );

    Ensure( n == numMoments );
    Ensure( result.size() == numMoments );
    return result;
}


} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of QuadServices.cc
//---------------------------------------------------------------------------//

