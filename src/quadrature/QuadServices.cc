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
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include "QuadServices.hh"

#include <iostream>

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor assumes numMoments == 1.
 * \param spQuad_     a smart pointer to a Quadrature object.
 * \param numMoments_ the number of moments used to represent the moment
 *                    space. 
 */
QuadServices::QuadServices( rtt_dsxx::SP< const Quadrature > spQuad_)
    : spQuad( spQuad_ ), 
      numMoments( spQuad->getNumAngles() ),
      Mmatrix( computeM() ),
      Dmatrix( computeD() )	
    { /* empty */ }

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the discrete-to-moment matrix. D = inverse(M).
 */
std::vector< double > QuadServices::computeD() const
{
    std::vector< double > Dmatrix( Mmatrix );

    return Dmatrix;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Create the M array (moment-to-discrete matrix).
 * 
 * \param name description
 * \return The moment-to-discrete matrix, cYnm(num_moments*num_angles)
 *
 * The M array is actually (2\ell+1)/sum_wts * cYnm(n,m).
 *
 * \f[
 * cY_{(\ell,k),m} = \left[ \frac{(\ell-\absk)!}{(\ell+\absk)!} (2 -
 * \delta_{k,0}) \right] Y_{\ell,k}(\Omega_m),
 * \f]
 *
 * where \f$ Y_{\ell,k}(\Omega_m) \f$ is the real part of the
 * \f$(\ell,k)^{th}\f$ spherical harmonic, evaluated at the 
 * \f$ m^{th} \f$ direction, \f$ \delta_{k,0} = 1 \f$ if \f$ k=0, \f$
 * otherwise it is 0, and \f$ (j)! \f$ represents a factorial. 
 *
 * The indexing of the cYnm array requires some additional description.  The
 * m index represents the M discrete directions defined by the quadrature
 * set.  The n index represents the index tuple \f$ (\ell,k) \f$ so that we
 * have the following ordering:
 *
 *      n    \ell    l
 *     -------------
 *      1    0    0
 *      2    1   -1
 *      3    1    0
 *      4    1    1
 *      5    2   -2  etc.
 */
std::vector< double > QuadServices::computeM() const
{
    unsigned const numAngles( spQuad->getNumAngles() );
    unsigned const dim(       spQuad->dimensionality() );

    std::vector< double > Mmatrix( numAngles*numMoments, -9999.0 );

    // loop over the (l,k) spherical harmonics indices. n = foo(l,k)
    // Use Morel's scheme for picking the (l,k) pairs. Ref: Jim Morel, "A
    // Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann
    // Transport Equation," NS&E 101, (1989).

    if( dim == 1 )
    {
	computeM_1D( Mmatrix );
	return Mmatrix;
    }

    // p84)  2-D XY
    else if( dim == 2 )
    {
	computeM_2D( Mmatrix );
	return Mmatrix;
    }

    // p85)  3-D XYZ
    else if( dim == 3 )
    {
	computeM_3D( Mmatrix );
	return Mmatrix;
    }

    else
    {
	std::ostringstream msg;
	msg << "FATAL ERROR in QuadServices::computeM().\n" 
	    << "\tThis class does not know how to build the "
	    << "moment-to-discrete matrix for quadrature set \"\n"
	    << spQuad->name() << "\"." << std::endl;
	    Insist( false, msg.str() );
    }

    // should never get here
    return Mmatrix;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the M matrix for 3D quadrature sets
 * \param Mmatrix the Moment-to-Discrete matrix
 * --------------------------------------------------------------
 * Use all moments for ell=0,...,N-1, k=-ell,...,ell
 * Plus                ell=N,         k<0
 * Plus                ell=N,         k>0, k odd
 * Plus                ell=N+1,       k<0, k even
 * ---------------------------------------------------------------
 */
void QuadServices::computeM_3D( std::vector< double > & Mmatrix ) const
{
    unsigned const snOrder(   spQuad->getSnOrder()   );
    unsigned n(0);

    // ell=0,...,(N-1) and k=-ell,...,ell
    for( unsigned ell=0; ell<snOrder; ++ell)
	for( int k(-1*ell); std::fabs(k) <= ell; ++k, ++n )
	    cYnm( Mmatrix, n, k, ell );
    { 
    // ell=N and k<0
	unsigned ell( snOrder );
	for( int k(-1*ell); k<0; ++k, ++n )
	    cYnm( Mmatrix, n, k, ell );

    // ell=N, k>0, k odd
	for( int k=1; k<=ell; k+=2, ++n )
	    cYnm( Mmatrix, n, k, ell );    
    }
    { 
    // ell=N+1 and k<0, k even
	unsigned ell( snOrder+1 );
	for( int k(-1*ell+1); k<0; k+=2, ++n )
	    cYnm( Mmatrix, n, k, ell );    
    }
    Ensure( n == numMoments );
    return;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the M matrix for 2D quadrature sets
 * \param Mmatrix the Moment-to-Discrete matrix
 * --------------------------------------------------------------
 * Use all moments for ell=0,...,N-1, k=0,...,ell
 * Plus                ell=N,         k=0,...,ell, k odd.
 * ---------------------------------------------------------------
 */
void QuadServices::computeM_2D( std::vector< double > & Mmatrix ) const
{
    unsigned const snOrder( spQuad->getSnOrder() );
    unsigned n(0);

    // ell=0,...,(N-1) and k=0,...,ell
    for( unsigned ell=0; ell<snOrder; ++ell)
	for( int k(0); k <= ell; ++k, ++n )
	    cYnm( Mmatrix, n, k, ell );
    { 
    // ell=N and k>0 and k odd
	unsigned ell( snOrder );
	for( int k=1; k<=ell; k+=2, ++n )
	    cYnm( Mmatrix, n, k, ell );
    }
    Ensure( n == numMoments );
    return;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the M matrix for 1D quadrature sets
 * \param Mmatrix the Moment-to-Discrete matrix
 * --------------------------------------------------------------
 * Use all moments for ell=0,...,N-1, k=0
 * ---------------------------------------------------------------
 */
void QuadServices::computeM_1D( std::vector< double > & Mmatrix ) const
{
    unsigned const snOrder( spQuad->getSnOrder() );
    unsigned n(0);

    // ell=0,...,(N-1) and k=0
    int k(0); // k is always zero for 1D.
    for( unsigned ell=0; ell<snOrder; ++ell, ++n)
	cYnm( Mmatrix, n, k, ell );

    Ensure( n == numMoments );
    return;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the (ell,k) spherical harmonic evaluated at Omega_m.
 * 
 * \param ell The ell index for the current spherical harmonic function.
 * \param k The k index for the current spherical harmonic function. 
 * \param m The index for the current discrete angle.
 * \return The (ell,k) sphercial harmonic evaluated at Omeaga_m. 
 */
double QuadServices::spherical_harmonic( unsigned const ell, 
					 int      const k,
					 unsigned const m ) const
{
    using rtt_units::PI;
    using rtt_dsxx::soft_equiv;

    Require( abs(k) <= ell );
    Require( ell < spQuad->getNumAngles() );
    Require( m < spQuad->getNumAngles() );

    unsigned const dim( spQuad->dimensionality() );
    double   const mu ( spQuad->getMu(m) );
    double   const eta( dim>1 ? spQuad->getEta(m) : 0.0 );
    double   const xi ( dim>2 ? spQuad->getXi(m)  : 0.0 );

    // Compute the azimuthal angle.
    double theta;
    if( soft_equiv( eta, 0.0 ) )
	theta = PI/2;
    else 
    {
	if( eta >= 0 )
	    theta = atan(xi/eta);
	else
	    theta = atan(xi/eta) + PI;
    }

    // ensure that theta is in the range 0...2*PI.
    if( theta < 0 )
	theta += 2*PI;
    Check( theta >= 0 );
    Check( theta <= 2*PI );
    
    return legendre_polynomial(ell,std::abs(k),mu) * cos( std::abs(k)*theta );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the (ell,k) associated Legendre polynomial evaluated at
 * mu_m. 
 * \param ell
 * \param k
 * \param mu_m
 * \return The (ell,k) associated Legendre polynomial evaluated at mu_m.
 * 
 * The associated Legendre polynomials, P_{ell,k}(x), are solutions to the
 * associated Legendre differential equation, where ell is a positive integer
 * and k = 0, ..., ell.  They can be given in terms of the unassociated
 * polynomials by:
 *
 * \f[
 * \begin{align}
 * P^k_ell(x) & = (-1)^k(1-x^2)^{k/2}\frac{d^k}{dx^k}P_ell(x) \\
 *             & = \frac{(-1)^k}{2^ellk
 *             !}(1-x^2)^{k/2}\frac{d^{ell+k}}{dx^{ell+k}}(x^2-1)^ell.
 * \end{align}
 * \f]
 *
 * where \f$ P_ell(x) \f$ are the unassociated Legendre polynomials.
 *
 * Associated polynomials are sometimes called Ferrers' functions (Sansone
 * 1991, p. 246). If m=0, they reduce to the unassociated polynomials.  The
 * associated Legendre functions are part of the <a
 * href="http://mathworld.wolfram.com/SphericalHarmonic.html">spherical
 * harmonics</a>, which are the solution of <a
 * href="http://mathworld.wolfram.com/LaplacesEquation.html">Laplace's
 * equation</a> in spherical coordinates.  They are orthogonal over W(x)=1
 * with the weighting funciton 1.
 * 
 * The associated Legendre Polynomials obey the following recurrence
 * relation:
 * 
 * \f[
 * (ell-k)P^k_ell(x) = x(2ell-1)P^k_{ell-1}(x) - (k+k-1)P^k_{ell-2}(x).
 * \f]
 *
 * An identity relating associated polynomials with negative k to the
 * corresponding functions with positive k is:
 *
 * \f[
 * P^{-k}_ell(x) = (-1)^k\frac{(ell-k)!}{(ell+k)!P^k_ell(x)}.
 * \f]
 *
 * Also,
 *
 * \f[
 * P^ell_ell(x) = (-1)^ell (2ell-1)!! (1-x^2)^{ell/2},
 * \f]
 *
 * and 
 *
 * \f[
 * P^ell_{ell+1}(x) = x (2ell+1)P^ell_ell(x).
 * \f]
 *
 * Including the factor of (-1)^k, the first few associated Legendre
 * polynomials are:
 * 
 * P^0_0(x) = 1
 * P^0_1(x) = x
 * P^1_1(x) = 1(1-x^2)^0.5
 * P^0_2(x) = 0.5 ( 3x^2 - 1 )
 * P^1_2(x) = -3x (1-x^2)^0.5
 * P^2_2(x) = 3(1-x^2)
 * etc.
 *
 * \sa http://mathworld.wolfram.com/LegendrePolynomial.html
 * 
 * \bug replace with gnu scientific library function gsl_sf_legendre ?
 */
double QuadServices::legendre_polynomial( unsigned const ell, 
					  unsigned const k,
					  double   const x ) const
{
    Require( k <= ell );

    // P_{0,0}(x) is defined to be 1.0.
    double pll(1.0);  // l==0

    if( k > 0 )
    {
	double somx2 = sqrt((1-x)*(1+x));
	double fact = 1.0;
	for( unsigned i=0; i<k; ++i )
	{
	    pll  *= -1.0*fact*somx2;
	    fact += 2.0;
	}
    }

    if( ell==k ) 
	return pll;
    else
    {
	// Compute P{k, k+1}
	double pllp1 = x * (2*k+1) * pll;
	
	if( ell == (k+1) )
	    return pllp1;
	else
	{
	    double pkk;
	    for( unsigned kk=k+2; kk < ell; ++kk )
	    {
		pkk = (x*(2*kk-1)*pllp1-(kk+k-1)*pll)/(kk-k);
		pll = pllp1;
		pllp1 = pkk;
	    }
	    return pkk;
	}
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the gamma(k,l) coefficient
 * 
 * \param name description
 * \return description
 */
double QuadServices::gamma( int k, unsigned ell ) const
{
    return std::sqrt( 2 - kronecker_delta(k,0) 
		      * factorial( ell - std::abs(k) )
		      / factorial( ell + std::abs(k) ) );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the beta(l) coefficient
 * 
 * \param name description
 * \return description
 */
double QuadServices::beta( unsigned ell ) const
{
    return ( 2*ell+1 ) / spQuad->getNorm();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the (n,m)-th entry of the Mmatrix.
 * 
 * \param Mmatrix the Moment-to-discrete matrix
 * \param n What column in the Mmatrix are we working on?
 * \param k The index n corresponds to the tuple (k,ell).
 * \param ell The index n corresponds to the tuple (k,ell).
 */
void QuadServices::cYnm( std::vector< double > & Mmatrix,
			   unsigned n, int k, unsigned ell ) const
{
    Check(n<numMoments);

    unsigned const numAngles( spQuad->getNumAngles() );
    double   const b( beta( ell ) );
    double   const g( gamma( k, ell ) );

    // Loop over all angles that use these values.
    for( unsigned m=0; m<numAngles; ++m )
    {
	Mmatrix[ m + n*numAngles ] 
	    = b * g * spherical_harmonic(ell,k,m); 
    }
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of QuadServices.cc
//---------------------------------------------------------------------------//
