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
 *
 * For \c dimension<3, the redundant moments are eliminated from the
 * expansion.  This means that only the \f$Y_l^0\f$ moments are included in
 * 1-D while only moments that are even in the azimuthal angle are included in
 * 2-D.
 *
 * \param dimensions Dimensionality of the simulation.
 * \param expansion_order Moment expansion order.
 *
 * \pre <code> dimensions==1 || dimensions==2 || dimensions==3
 * </code>
 */
/* static */
unsigned
QuadServices::compute_number_of_moments(unsigned const dimensions,
                                        unsigned const expansion_order)
{
    Require(dimensions==1 || dimensions==2 || dimensions==3);

    switch(dimensions)
    {
        case 1:
            return expansion_order+1;
        case 2:
            return (expansion_order+1)*(expansion_order+2)/2;
        case 3:
            return (expansion_order+1)*(expansion_order+1);
        default:
            Insist(false, "bad case");
            return 0;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor builds square D and M operators using Morel's
 * Galerkin-Sn heuristic. 
 * \param spQuad_ a smart pointer to a Quadrature object.
 * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
 */
QuadServices::QuadServices( rtt_dsxx::SP< const Quadrature > const spQuad_,
                            QIM const qm,
                            unsigned const expansionOrder,
                            comparator_t comparator_ )
    : spQuad( spQuad_ ),
      qm( qm ),
      n2lk( compute_n2lk( expansionOrder ) ),
      numMoments( n2lk.size() ),
      ordinates( compute_ordinates(spQuad_,
                                   comparator_) ),
      Mmatrix( computeM() ),
      Dmatrix( computeD() )
{ 
    using rtt_dsxx::soft_equiv;
    using rtt_units::PI;

    Remember( double const mu2( std::sqrt(3.0)/3.0 ); );
    Ensure( soft_equiv( compute_azimuthalAngle( 1.0, 0.0, 0.0 ), 0.0 ) );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, mu2, mu2 ), PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, -1.0*mu2, mu2 ), 7.0*PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, -1.0*mu2, -1.0*mu2 ), 7.0*PI/4.0 )  );
    Ensure( soft_equiv( compute_azimuthalAngle( mu2, mu2, -1.0*mu2 ), PI/4.0 )  );

    Check( soft_equiv(gsl_sf_legendre_Plm( 0, 0, 0.5 ), 1.0 ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 1, 0, 0.5 ), 0.5 ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 1, 1, mu2 ), -1.0*std::sqrt(1.0-mu2*mu2) ));
    Check( soft_equiv(gsl_sf_legendre_Plm( 2, 2, mu2 ), 3.0*(1.0-mu2*mu2) ));

    if( qm == GALERKIN ) Ensure( D_equals_M_inverse() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor that allows the user to pick the (k,l) moments to use.
 * \param spQuad_     a smart pointer to a Quadrature object.
 * \param lkMoments_  vector of tuples that maps from index n to (k,l).
 * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
 */
QuadServices::QuadServices( rtt_dsxx::SP< const Quadrature > const spQuad_,
                            std::vector< lk_index > const & lkMoments_,
                            QIM const   qm,
                            comparator_t comparator_ )
    : spQuad( spQuad_ ),
      qm( qm ),
      n2lk( lkMoments_ ),
      numMoments( n2lk.size() ),
      ordinates( compute_ordinates(spQuad_,
                                   comparator_) ),
      Mmatrix( computeM() ),
      Dmatrix( computeD() )
{ 
    Ensure( D_equals_M_inverse() );
}

// //---------------------------------------------------------------------------//
// /*!
//  * \brief Constructor that allows construction from an OrdinateSet.
//  * \param os_ An OrdinateSet
//  * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
//  */
// QuadServices::QuadServices( OrdinateSet const & os )
//     : spQuad( os.getQuadrature() ),   
//       numMoments( spQuad->getNumOrdinates() ),    
//       n2lk(       compute_n2lk() ),
//       Mmatrix(    computeM() ),
//       Dmatrix(    computeD() )	
// {
//     Ensure( D_equals_M_inverse() );
// }


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

    size_t const numOrdinates( spQuad->getNumOrdinates() );
    std::vector< double > psi( numOrdinates, 0.0 );

    for( size_t m=0; m<numOrdinates; ++m )
	for( size_t n=0; n<numMoments; ++n )
	    psi[ m ] += Mmatrix[ n + m*numMoments ] * phi[n];
    
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
std::vector< double > QuadServices::applyD( std::vector< double > const & psi ) const
{
    size_t const numOrdinates( spQuad->getNumOrdinates() );
    Require( psi.size() == numOrdinates );

    std::vector< double > phi( numMoments, 0.0 );

    for( size_t m=0; m<numOrdinates; ++m )
        for( size_t n=0; n<numMoments; ++n )
            phi[ n ] += Dmatrix[ m + n*numOrdinates ] * psi[m];
    
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
 *
 * Normally, M will not be square because we only have isotropic scatting.
 * For isotropic scattering M will be (numOrdinates x 1 moment).  We will use the
 * Moore-Penrose Pseudo-Inverse Matrix, \f$ D = (M^T * M)^-1 * M^T.\f$
 */
std::vector< double > QuadServices::computeD(void) const
{
    if( qm == GALERKIN ) return computeD_morel();
    if( qm == SN )       return computeD_traditional();
    if( qm == SVD )      return computeD_svd();

    Check(spQuad->getNumOrdinates() == ordinates.size());

    // Should never get here.
    Insist( qm == GALERKIN || qm == SN || qm == SVD,
            "qm has an unknown value!");
    return std::vector<double>();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the discrete-to-moment matrix. 
 *
 * Computes \f$ \mathbf{D} \equiv \mathbf{M}^{-1} \f$.  This private function
 * is called by the constuctor.
 *
 * Under Morel's method, M is always square and we do a direct inversion to
 * obtain D, \f$ D = M^{-1} \f$.
 */
std::vector< double > QuadServices::computeD_morel(void) const
{
    int n( numMoments );
    int m( spQuad->getNumOrdinates() );

    Require( n == m );

    // create a copy of Mmatrix to use a temp space.
    std::vector< double > M( Mmatrix );
    std::vector< double > D( m*n );

    // Create GSL matrix views of our M and D matrices.
    // LU will get a copy of M.  This matrix will be decomposed into LU. 
    gsl_matrix_view gsl_M = gsl_matrix_view_array( &M[0], m, n );
    gsl_matrix_view gsl_D = gsl_matrix_view_array( &D[0], n, m );
    
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

    // Store the LU decomposition in the matrix M.
    int result = gsl_linalg_LU_decomp( &gsl_M.matrix, p, &signum );
    Check( result == 0 );
    Check( diagonal_not_zero( M, n, m ) );

    // Compute the inverse of the matrix LU from its LU decomposition (LU,p),
    // storing the results in the matrix Dmatrix.  The inverse is computed by
    // solving the system (LU) x = b for each column of the identity matrix.

    result = gsl_linalg_LU_invert( &gsl_M.matrix, p, &gsl_D.matrix );
    Check( result == 0 );

    // Free the space reserved for the permutation matrix.
    gsl_permutation_free( p );

    return D;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the discrete-to-moment matrix. 
 *
 * This function does not attempt to ensure \f$ D = M^{-1} \f$.  Instead, it
 * simply computes the D so that,
 * \f[
 *    \Phi_{(l,k)} = \int\limits_{4\pi}{d\Omega c_{(l,k)} Y_{(l,k)}(\Omega)}.
 * \f]
 */
std::vector< double > QuadServices::computeD_traditional(void) const
{
    using rtt_sf::galerkinYlk;

    unsigned const numOrdinates( spQuad->getNumOrdinates() );
    double   const sumwt(     spQuad->getNorm() );
    unsigned const dim(       spQuad->dimensionality() );
    
    std::vector< double > D( numOrdinates*numMoments );
    
    for( unsigned m=0; m<numOrdinates; ++m )
    {
        double mu( ordinates[m].mu() );
        double wt( ordinates[m].wt() );

        for( unsigned n=0; n<numMoments; ++n )
        {
            unsigned const ell ( n2lk[n].first  );
            int      const k   ( n2lk[n].second );  
            // Must mult by (sumwt/(2ell+1)) twice to back out this
            // coefficient.
//            double c( sumwt*sumwt/(2*ell+1)/(2*ell+1) );
            double c( sumwt/(2*ell+1) );

            if( dim == 1 ) // 1D mesh, 1D quadrature
            { // for 1D, mu is the polar direction and phi == 0, k==0
                D[ m + n*numOrdinates ] = c*wt*galerkinYlk( ell, std::abs(k), mu, 0.0, sumwt );
            }

            else if( dim == 2 ) // 2D mesh, 2D quadrature
            { // for 2D, mu is taken to be the polar direction.
                double eta( ordinates[m].eta() );
                double xi ( ordinates[m].xi() );
                double phi( compute_azimuthalAngle( xi, eta, mu ) );
                D[ m + n*numOrdinates ] = c*wt*galerkinYlk( ell, k, mu, phi, sumwt );
            }
            else // 3D mesh, 3D quadrature
            {
                Check( dim == 3);
                double eta( ordinates[m].eta() );
                double xi ( ordinates[m].xi() );
                double phi( compute_azimuthalAngle( mu, eta, xi ) );
                D[ m + n*numOrdinates ] = c*wt*galerkinYlk( ell, k, xi, phi, sumwt );
            }
            
        } // n: end moment loop
    } // m: end ordinate loop

    return D;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the discrete-to-moment matrix. 
 *
 * Computes \f$ \mathbf{D} \equiv \mathbf{M}^{-1} \f$.  This private function
 * is called by the constuctor.
 *
 * Normally, M will not be square because we only have isotropic scatting.
 * For isotropic scattering M will be (numOrdinates x 1 moment).  We will use the
 * Moore-Penrose Pseudo-Inverse Matrix, \f$ D = (M^T * M)^-1 * M^T.\f$
 */
std::vector< double > QuadServices::computeD_svd(void) const
{
    int n( numMoments );
    int m( spQuad->getNumOrdinates() );

    // SVD:
    //! \f$ M = U S V^T \f$
    
    // create a copy of Mmatrix to use a temp space.
    std::vector< double > M( Mmatrix );
    std::vector< double > V( n*n );
    std::vector< double > D( m*n );
    std::vector< double > S( n );

    // Create GSL matrix views of our M and D matrices.
    // LU will get a copy of M.  This matrix will be decomposed into LU. 
    gsl_matrix_view gsl_M = gsl_matrix_view_array( &M[0], m, n );
    gsl_vector_view gsl_S = gsl_vector_view_array( &S[0], n );
    gsl_matrix_view gsl_V = gsl_matrix_view_array( &V[0], n, n );
    gsl_matrix_view gsl_D = gsl_matrix_view_array( &D[0], n, m );
    
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
 * \brief Create the M array (moment-to-discrete matrix).
 * \return The moment-to-discrete matrix.
 *
 * This private member function is called by the constructor. 
 * 
 * The moment-to-discrete matrix will be num_moments by num_ordinates in size.
 */
std::vector< double > QuadServices::computeM(void) const
{
    using std::sqrt;
    using rtt_sf::galerkinYlk;

    unsigned const numOrdinates( spQuad->getNumOrdinates() );
    unsigned const dim(       spQuad->dimensionality() );
    double   const sumwt(     spQuad->getNorm() );

    Check(numOrdinates == ordinates.size());

    // resize the M matrix.
    std::vector< double > Mmatrix( numMoments*numOrdinates, -9999.0 );

    for( unsigned n=0; n<numMoments; ++n )
    {
        for( unsigned m=0; m<numOrdinates; ++m )
        {
            unsigned const ell ( n2lk[n].first  );
            int      const k   ( n2lk[n].second ); 
        
            if( dim == 1 ) // 1D mesh, 1D quadrature
            { // for 1D, mu is the polar direction and phi == 0, k==0
                double mu ( ordinates[m].mu() );
                Mmatrix[ n + m*numMoments ] = galerkinYlk( ell, k, mu, 0.0, sumwt );
            }
            else if( dim == 2 ) // 2D mesh, 2D quadrature
            {
                // for 2D, mu is taken to be the polar direction.
                // xi is always positive (a half-space).
                //! \todo this is the same computation as Ordinate.cc::Y(l,k,Ordinate,norm). Try to prevent code duplication.

                double mu ( ordinates[m].mu() );
                double eta( ordinates[m].eta() );
                double xi(  ordinates[m].xi() );
                double phi( compute_azimuthalAngle( xi, eta, mu ) );
                Mmatrix[ n + m*numMoments ] = galerkinYlk( ell, k, mu, phi, sumwt );
            }
            else // 3D mesh, 3D quadrature
            {
                Check( dim == 3);
                double mu ( ordinates[m].mu()  );
                double eta( ordinates[m].eta() );
                double xi ( ordinates[m].xi() );
                double phi( compute_azimuthalAngle( mu, eta, xi ) );
                Mmatrix[ n + m*numMoments ] = galerkinYlk( ell, k, xi, phi, sumwt );
            } 
        } // n: end moment loop
    } // m: end ordinate loop
    return Mmatrix;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the Azimuthal angle for the current quadrature direction.
 */
double QuadServices::compute_azimuthalAngle( double const mu,
					     double const eta,
					     double const Remember(xi) ) 
{
    using rtt_units::PI;
    using rtt_dsxx::soft_equiv;

    Require( std::abs(mu)  <= 1.0 );
    Require( std::abs(eta) <= 1.0 );
    Require( std::abs(xi)  <= 1.0 );

    // For 1D sets, we define this angle to be zero.
    if( soft_equiv( eta, 0.0 ) ) return 0.0;

    double azimuthalAngle ( std::atan2( eta, mu ) );
    if( azimuthalAngle < 0.0 )
        azimuthalAngle += 2.0*PI;
    
//     // For 2D sets, reconstruct xi from known information: 
//     // xi*xi = 1.0 - eta*eta - mu*mu
//     // Always use positive value for xi.
//     double local_xi( xi );
//     if( soft_equiv( local_xi,  0.0 ) )
// 	local_xi = std::sqrt( 1.0 - mu*mu - eta*eta );

//     double azimuthalAngle(999.0);

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

    // ensure that theta is in the range 0...2*PI.
    Ensure( azimuthalAngle >= 0 );
    Ensure( azimuthalAngle <= 2*PI );
    
    return azimuthalAngle;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Multiply M and D and compare the result to the identity matrix.
 * \return true if M = D^(-1), otherwise false.
 */
bool QuadServices::D_equals_M_inverse() const
{
    using rtt_dsxx::soft_equiv;

    unsigned n( numMoments );
    int m( spQuad->getNumOrdinates() );
//    int nm( std::min( n,m ) );

    // create non-const versions of M and D.
    std::vector< double > Marray( Mmatrix );
    std::vector< double > Darray( Dmatrix );
    std::vector< double > Iarray( n*n, -999.0 );
    gsl_matrix_view M = gsl_matrix_view_array( &Marray[0], m, n );
    gsl_matrix_view D = gsl_matrix_view_array( &Darray[0], n, m );
    gsl_matrix_view I = gsl_matrix_view_array( &Iarray[0], n, n );
    
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
    
    gsl_blas_dgemm( op, op, alpha, &D.matrix, &M.matrix, beta, &I.matrix );

    for( unsigned i=0; i<n; ++i )
	if( ! soft_equiv( Iarray[ i + i*n ], 1.0 ) ) return false;

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
std::vector< QuadServices::lk_index > QuadServices::
compute_n2lk( unsigned const expansionOrder ) const
{
    unsigned const dim( spQuad->dimensionality() );
    unsigned const L( expansionOrder+1 );

    if( dim == 3 )
    {
        if( qm == GALERKIN )  
            return compute_n2lk_3D_morel();
        else
            return compute_n2lk_3D_traditional(L);
    }

    Check( dim < 3 );
    if( dim == 2 )
    {
        if( qm == GALERKIN ) 
            return compute_n2lk_2D_morel();
        else
            return compute_n2lk_2D_traditional(L);
    }
    
    Check( dim == 1 );
    if( dim == 1 )
    {
        if( qm == GALERKIN )
            return compute_n2lk_1D( spQuad->getNumOrdinates() );
        else
            return compute_n2lk_1D(L);
    }

    // Should never get here.
    Insist( dim <= 3 && dim >= 1, "I only know about dim = {1,2,3}.");
    return std::vector< QuadServices::lk_index >();
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::
compute_n2lk_3D_morel( void ) const
{
    int const L( spQuad->getSnOrder() );

    // This algorithm only  works for level symmetric sets because it
    // assumes numOrdinates = (L)(L+2).
    Require( static_cast<int>(spQuad->getNumOrdinates()) == L*(L+2) );
    
    std::vector< lk_index > result;

    // Choose: l= 0, ..., L-1, k = -l, ..., l
    for( int ell=0; ell< L; ++ell )
	for( int k = -ell; k <= ell; ++k )
	    result.push_back( lk_index(ell,k) );

    // Add ell=L and k<0
    {
	unsigned ell( L );
	for( int k(-1*static_cast<int>(ell)); k<0; ++k )
	    result.push_back( lk_index(ell,k) );
    }

    // Add ell=L, k>0, k odd
    {
	int ell( L );
	for( int k=1; k<=ell; k+=2 )
	    result.push_back( lk_index(ell,k) );
    }

    // Add ell=L+1 and k<0, k even
    {
	unsigned ell( L+1 );
	for( int k(-1*static_cast<int>(ell)+1); k<0; k+=2 )
	    result.push_back( lk_index(ell,k) );
    }

    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::
compute_n2lk_2D_morel( void ) const
{
    int const L( spQuad->getSnOrder() );

    // This algorithm only  works for level symmetric sets because it
    // assumes numOrdinates = (L)(L+2)/2.
    Require( static_cast<int>(spQuad->getNumOrdinates()) == L*(L+2)/2 );

    std::vector< lk_index > result;
    
    // Choose: l= 0, ..., N-1, k = 0, ..., l
    for( int ell=0; ell<L; ++ell )
	for( int k=0; k<=ell; ++k )
	    result.push_back( lk_index(ell,k) );

    // Add ell=N and k>0, k odd
    int ell( L );
    for( int k=1; k<=ell; k+=2 )
        result.push_back( lk_index(ell,k) );
    
    return result;
}


//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::
compute_n2lk_3D_traditional( unsigned const L ) const
{
    std::vector< lk_index > result;

    // Choose: l= 0, ..., L, k = -l, ..., l
    for( int ell=0; ell<static_cast<int>(L); ++ell )
	for( int k(-1*static_cast<int>(ell)); std::abs(k) <= ell; ++k )
	    result.push_back( lk_index(ell,k) );

    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::
compute_n2lk_2D_traditional( unsigned const L ) const
{
    std::vector< lk_index > result;
    
    // Choose: l= 0, ..., N, k = 0, ..., l
    for( int ell=0; ell<static_cast<int>(L); ++ell )
	for( int k=0; k<=ell; ++k )
	    result.push_back( lk_index(ell,k) );

    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index > QuadServices::
compute_n2lk_1D( unsigned const L ) const
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
double QuadServices::augmentM( unsigned n, Ordinate const &Omega ) const
{
    // If you trigger this exception, you may have requested too many
    // moments.  Your quadrature set must have more ordinates than the number of
    // moments requested.
    Require(n<n2lk.size());
    // The n-th moment is the (l,k) pair used to evaluate Y_{l,k}.
    //return Ordinate::Y( n2lk[n].first, n2lk[n].second, Omega, spQuad->getNorm() );

    using rtt_sf::galerkinYlk;

    double mu(Omega.mu());
    double xi(Omega.eta());
    double eta(Omega.xi());
    double phi( compute_azimuthalAngle( xi, eta, mu ) );
    return galerkinYlk( n2lk[n].first, n2lk[n].second, mu, phi, spQuad->getNorm());
}

std::vector< Ordinate > 
QuadServices::compute_ordinates(  rtt_dsxx::SP< const Quadrature > const spQuad,
                                  comparator_t const comparator) const
{
    unsigned const numOrdinates( spQuad->getNumOrdinates() );
    unsigned const dim( spQuad->dimensionality() );    

    std::vector< Ordinate > Result(numOrdinates);

    if( dim == 1 ) 
    { 
        for (unsigned m=0; m<numOrdinates; m++)
        {
            double const mu = spQuad->getMu(m);
            double const wt = spQuad->getWt(m);
            Result[m] = Ordinate(mu, 0.0, 0.0, wt);
        }
    }
    else if ( dim == 2 ) 
    {
        if( spQuad->getEta().empty() )
        {
            for (unsigned m=0; m<numOrdinates; m++)
            {
                double const mu = spQuad->getMu(m);
                double const xi = spQuad->getXi(m);
                double const eta  = sqrt(1.0 - mu*mu - xi*xi);
                double const wt = spQuad->getWt(m);
                Result[m] = Ordinate(mu, eta, xi, wt);
            }
        }
        else
        {
            for (unsigned m=0; m<numOrdinates; m++)
            {
                double const mu = spQuad->getMu(m);
                double const eta = spQuad->getEta(m);
                double const xi = sqrt(1.0 - mu*mu - eta*eta);
                double const wt= spQuad->getWt(m);
                Result[m] = Ordinate(mu, eta, xi, wt);
            }
        }
    }
    else if ( dim == 3)
    {
        for (unsigned m=0; m<numOrdinates; m++)
        {
            double const mu = spQuad->getMu(m);
            double const eta = spQuad->getEta(m);
            double const xi = spQuad->getXi(m);
            double const wt = spQuad->getWt(m);
            Result[m] = Ordinate(mu, eta, xi, wt);
        }
    }
        
    std::sort( Result.begin(), Result.end(), comparator);        

    return Result;
} 

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of QuadServices.cc
//---------------------------------------------------------------------------//

