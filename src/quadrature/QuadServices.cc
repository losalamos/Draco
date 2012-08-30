//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices.cc
 * \author Kelly Thompson
 * \date   Mon Nov  8 11:17:12 2004
 * \brief  Provide Moment-to-Discrete and Discrete-to-Moment operators.
 * \note   © Copyright 2006 LANSLLC All rights reserved. 
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
#include "ds++/Array.hh"
#include "special_functions/Factorial.hh"
#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"

#include "QuadServices.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor builds square D and M operators using Morel's
 * Galerkin-Sn heuristic. 
 * \param spQuad_ a smart pointer to a Quadrature object.
 * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
 */
QuadServices::QuadServices( vector<Ordinate> const &ordinates,
                            double const norm,
                            unsigned const dimension,
                            unsigned const expansionOrder)
    : ordinates_(ordinates),
      norm_(norm),
      dimension_(dimension),
      expansionOrder_(expansionOrder)
{ 
    using rtt_dsxx::soft_equiv;
    using rtt_units::PI;

    Remember( double const mu2( std::sqrt(3.0)/3.0 ); );

//    std::cout << " + + " << compute_azimuthalAngle(  mu2,  mu2, 0.0 ) << std::endl;
//    std::cout << " + - " << compute_azimuthalAngle(  mu2, -mu2, 0.0 ) << std::endl;
//    std::cout << " - + " << compute_azimuthalAngle( -mu2,  mu2, 0.0 ) << std::endl;
//    std::cout << " - - " << compute_azimuthalAngle( -mu2, -mu2, 0.0 ) << std::endl;

    Ensure( soft_equiv( compute_azimuthalAngle( 1.0, 0.0, 0.0 ), PI ) );
    Ensure( soft_equiv( compute_azimuthalAngle(  mu2,  mu2, 0.0 ), PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle( -mu2,  mu2, 0.0 ), 3.0*PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle(  mu2, -mu2, 0.0 ), 7.0*PI/4.0 )  );

    Check( soft_equiv(gsl_sf_legendre_Plm( 0, 0, 0.5 ), 1.0 ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 1, 0, 0.5 ), 0.5 ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 1, 1, mu2 ), -1.0*std::sqrt(1.0-mu2*mu2) ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 2, 2, mu2 ), 3.0*(1.0-mu2*mu2) ));
}

//---------------------------------------------------------------------------//
/*! \brief Return the moment-to-discrete operator.
 */
//---------------------------------------------------------------------------//

std::vector< double > QuadServices::getM() const
{
    std::vector< double > M(getM_());
    std::vector< double > Result(M);

    return Result;
}

//---------------------------------------------------------------------------//
/*! \brief Return the discrete-to-moment operator.
 */
//---------------------------------------------------------------------------//
std::vector< double > QuadServices::getD() const
{
    std::vector< double > D(getD_());
    std::vector< double > Result(D);

    return Result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Apply the action of \f$ M \f$ to the moment-based solution vector
 * \f$ \mathbf\Phi \f$.  That is, \f$ \mathbf\Psi = \mathbf{M}\mathbf\Phi \f$
 * 
 * \param phi Moment based solution vector, \f$ \mathbf\Phi \f$.
 * \return The discrete angular flux, \f$ \mathbf\Psi \f$.
 */
std::vector< double > QuadServices::applyM( std::vector< double > const & phi) const
{
    size_t const numOrdinates( ordinates_.size() );
    std::vector< double > psi( numOrdinates, 0.0 );

    std::vector< double > M(getM());

    unsigned const numMoments(getNumMoments());
    Require( phi.size() == numMoments );
    
    for( size_t m=0; m<numOrdinates; ++m )
        for( size_t n=0; n<numMoments; ++n )
            psi[ m ] += M[ n + m*numMoments ] * phi[n];
    
    return psi;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Apply the action of \f$ D \f$ to the discrete-ordinate based solution vector
 * \f$ \mathbf\Psi \f$. That is, \f$ \mathbf\Phi = \mathbf{D}\mathbf\Psi \f$.
 * 
 * \param psi Discrete ordinate-based solution vector, \f$ \mathbf\Psi \f$.
 * \return The moment-based solution vector, \f$ \mathbf\Phi \f$.
 */
std::vector< double > QuadServices::applyD( std::vector< double > const & psi) const
{
    size_t const numOrdinates( ordinates_.size() );
    Require( psi.size() == numOrdinates );

    std::vector< double > D(getD());
    unsigned const numMoments(getNumMoments());

    std::vector< double > phi;

    phi.resize( numMoments, 0.0 );
    
    for( size_t m=0; m<numOrdinates; ++m )
        for( size_t n=0; n<numMoments; ++n )
            phi[ n ] += D[ m + n*numOrdinates ] * psi[m];

    return phi;
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the discrete-to-moment matrix. 
 *
 * Normally, M will not be square because we only have isotropic scatting.
 * For isotropic scattering M will be (numOrdinates x 1 moment).  We will use the
 * Moore-Penrose Pseudo-Inverse Matrix, \f$ D = (M^T * M)^-1 * M^T.\f$
 */

/* static */
std::vector< double > QuadServices::computeD_SVD(std::vector<Ordinate> const &ordinates,
                                                 std::vector< lk_index > const &n2lk,
                                                 std::vector<double> const &mM,
                                                 unsigned const,
                                                 double const)
{
    int n = n2lk.size();
    int m = ordinates.size();

    // SVD:
    //! \f$ M = U S V^T \f$
    
    // create a copy of Mmatrix to use a temp space.
    std::vector< double > M( mM );
    std::vector< double > V( n*n );
    std::vector< double > D( m*n );
    std::vector< double > S( n );

    // Create GSL matrix views of our M and D matrices.
    // LU will get a copy of M.  This matrix will be decomposed into LU. 
    gsl_matrix_view gsl_M = gsl_matrix_view_array( &M[0], m, n );
    gsl_vector_view gsl_S = gsl_vector_view_array( &S[0], n );
    gsl_matrix_view gsl_V = gsl_matrix_view_array( &V[0], n, n );
    // gsl_matrix_view gsl_D = gsl_matrix_view_array( &D[0], n, m );
    
    // Singular Value Decomposition
    //
    // A general rectangular m-by-n matrix M has a singular value
    // decomposition (svd) into the product of an m-by-n orthogonal matrix U,
    // an n-by-n diagonal matrix of singular values S and the transpose of an
    // n-by-n orthogonal square matrix V,
    //
    //      M = U S V^T
    //
    // The singular values \sigma_i = S_{ii} are all non-negative and are
    // generally chosen to form a non-increasing sequence
    // \sigma_1 >= \sigma_2 >= ... >= \sigma_n >= 0.
    //
    // The singular value decomposition of a matrix has many practical
    // uses. The condition number of the matrix is given by the ratio of the
    // largest singular value to the smallest singular value. The presence of
    // a zero singular value indicates that the matrix is singular. The number
    // of non-zero singular values indicates the rank of the matrix. In
    // practice singular value decomposition of a rank-deficient matrix will
    // not produce exact zeroes for singular values, due to finite numerical
    // precision. Small singular values should be edited by choosing a
    // suitable tolerance.
    //
    // This function factorizes the m-by-n matrix M into the singular value
    // decomposition M = U S V^T for m >= n. On output the matrix M is
    // replaced by U. The diagonal elements of the singular value matrix S are
    // stored in the vector S. The singular values are non-negative and form a
    // non-increasing sequence from S_1 to S_n. The matrix V contains the
    // elements of V in untransposed form. To form the product U S V^T it is
    // necessary to take the transpose of V. M workspace of length n is
    // required in work.
    //
    // This routine uses the Golub-Reinsch SVD algorithm.
    //
    
    // returns U in storage "M"
    {
        gsl_vector    * gsl_w = gsl_vector_alloc(n);
        Remember(int result = )
            gsl_linalg_SV_decomp( &gsl_M.matrix, &gsl_V.matrix,
                                  &gsl_S.vector, gsl_w );
        gsl_vector_free(gsl_w);
        Check( result == 0 );
    }

    //
    // Now D = V S^-1 M^T
    //
    
    // Invert S
    for( int nn=0; nn<n; ++nn )
        S[nn]=1.0/S[nn];

    // Now, D = V S M^T

    // mult V(nxn) by S(diag,n)

    for( int nn=0; nn<n; ++nn )
        for( int nnn=0; nnn<n; ++nnn )
            V[nnn+n*nn] *= S[nnn];
    
    // mult VS by U^T

    for( int nn=0; nn<n; ++nn )
    {
        for( int mm=0; mm<m; ++mm )
        {
            double sum(0);
            for( int nnn=0; nnn<n; ++nnn )
            {
                sum += V[nnn+n*nn] * M[nnn+mm*n];
            }
            D[m*nn+mm] = sum;
        }
    }

    return D;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Check the diagonal of a matrix
 * \return true if all entries are non-zero.
 *
 * This private member function used by DBC routines in computeD.
 */
bool QuadServices::diagonal_not_zero( std::vector<double> const & vec,
                                      int m, int n ) 
{
    int dim( std::min(m,n) );
    for(int i=0;i<dim;++i)
        if( rtt_dsxx::soft_equiv(vec[i+i*m],0.0) )
            return false;
    return true;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the Azimuthal angle for the current quadrature direction.
 */
double QuadServices::compute_azimuthalAngle( double const mu,
					     double const eta,
					     double const Remember(xi) ) const
{
    using rtt_units::PI;
    using rtt_dsxx::soft_equiv;

    Require( std::abs(mu)  <= 1.0 );
    Require( std::abs(eta) <= 1.0 );
    Require( std::abs(xi)  <= 1.0 );

    if( soft_equiv( eta, 0.0 ) ) return PI;

//-------------------------------------------------------------

    // It is important to remember that the positive mu axis points to the
    // left and the positive eta axis points up, when the unit sphere is
    // projected on the plane of the mu- and eta-axis. In this case, phi is
    // measured from the mu-axis counterclockwise.
    //
    // This accounts for the fact that the aziumuthal angle is discretized
    // on levels of the xi-axis, making the computation of the azimuthal angle
    // here consistent with the discretization by using the eta and mu
    // ordinates to define phi.
    
    double azimuthalAngle ( std::atan2( eta, mu) );
    
    if( azimuthalAngle < 0.0 )
        azimuthalAngle += 2*PI;
    
//-------------------------------------------------------------
//     // For 2D sets, reconstruct xi from known information: 
//     // xi*xi = 1.0 - eta*eta - mu*mu
//     // Always use positive value for xi.
//     double local_xi( xi );
//     if( soft_equiv( local_xi,  0.0 ) )
// 	local_xi = std::sqrt( 1.0 - mu*mu - eta*eta );
//
//     double azimuthalAngle(999.0);
//
//     if( local_xi > 0.0 )
//     {
// 	if( eta > 0.0 )
// 	    azimuthalAngle = std::atan(xi/eta);
// 	else
// 	    azimuthalAngle = PI - std::atan(xi/std::abs(eta));
//     } 
//     else 
//     {
// 	if( eta > 0 )
// 	    azimuthalAngle = 2*PI - std::atan(std::abs(xi)/eta);
// 	else
// 	    azimuthalAngle = PI + std::atan(xi/eta);
//     }
//-------------------------------------------------------------

    // Ensure that theta is in the range 0...2*PI.
    Ensure( azimuthalAngle >= 0 );
    Ensure( azimuthalAngle <= 2*PI );
    
    return azimuthalAngle;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */

std::vector< QuadServices::lk_index >
QuadServices::compute_n2lk_1D( unsigned const L )
{
    std::vector< lk_index > result;
    
    // Choose: l= 0, ..., L-1, k = 0
    int k(0); // k is always zero for 1D.
    
    for( unsigned ell=0; ell<L; ++ell )
	result.push_back( lk_index(ell,k) );

    return result;
}
 
//---------------------------------------------------------------------------//
/*! 
 * \brief Computes a list of the number of moments at each spherical harmonics
 * index, based on the maximum expansion order found.
 */

std::vector<unsigned> QuadServices::compute_moments(unsigned const L,
                                                    std::vector< QuadServices::lk_index > const &n2lk) 
{
    Insist(!n2lk.empty(), "n2lk is unexpectedly empty.");

    std::vector<unsigned> moments(L+1, 0);
    for(unsigned n=0; n<n2lk.size(); ++n)
    {
        unsigned const l(n2lk[n].first); 
        moments[l] += 1;
    }

    //for(unsigned l=0; l<=L; ++l)
    //{
    //    std::cout << " moment[" << l << "]: " << moments[l] <<  std::endl;
    //}

    return moments;
}

unsigned QuadServices::max_available_expansion_order(std::vector< QuadServices::lk_index > const &n2lk) 
{
    Insist(!n2lk.empty(), "n2lk is unexpectedly empty.");

    unsigned L=0;
    for(unsigned n=0; n<n2lk.size(); ++n)
    {
        unsigned const l(n2lk[n].first); 
        L = std::max(L,l);
    //    std::cout << " n = " << n << " ... l=" << n2lk[n].first << ", k=" << n2lk[n].second << std::endl;
    }

    //std::cout << " maximum available expansion order is " << L << std::endl;

    return L;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute extra entries for M (for starting directions).
 * \param n Spherical Harmonic moment n -> (l,k).
 * \param Omega A discrete direction ordinate.
 * \return An appropriate value for augmenting M.
 *
 * These values are not used in determing D.  They only extent M for starting
 * directions.  For example, for S2 RZ there are 2 additional starting
 * directions.  The S2 XY M operator has 16 values (4 ordinates and 4 moments).
 * The starting directions will have the same moment values as the first 4
 * ordinates (Y00, Y10, Y11, Y21), but will be evaluated as the starting
 * direction ordinates.
 */
double QuadServices::augmentM( unsigned n, Ordinate const &Omega,  std::vector< lk_index > const &n2lk ) const
{
    // If you trigger this exception, you may have requested too many
    // moments.  Your quadrature set must have more ordinates than the number of
    // moments requested.
    Require(n<n2lk.size());

    // The n-th moment is the (l,k) pair used to evaluate Y_{l,k}.
    //return Ordinate::Y( n2lk[n].first, n2lk[n].second, Omega, spQuad->getNorm() );

    using rtt_sf::galerkinYlk;

    double mu(Omega.mu());
    double xi(Omega.xi());
    double eta(Omega.eta());
    double phi( compute_azimuthalAngle( xi, eta, mu ) );
    return galerkinYlk( n2lk[n].first, n2lk[n].second, mu, phi, norm_);
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of QuadServices.cc
//---------------------------------------------------------------------------//

