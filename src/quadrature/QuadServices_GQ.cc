//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices_GQ.cc
 * \author Kelly Thompson
 * \date   Mon Nov  8 11:17:12 2004
 * \brief  Provide Moment-to-Discrete and Discrete-to-Moment operators.
 * \note   © Copyright 2006 LANSLLC All rights reserved. 
 */
//---------------------------------------------------------------------------//
// $Id: QuadServices_GQ.cc 6499 2012-03-15 20:19:33Z kgbudge $
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

#include "QuadServices_GQ.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor builds square D and M operators using Morel's
 * Galerkin-Sn heuristic. 
 * \param spQuad_ a smart pointer to a Quadrature object.
 * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
 */
QuadServices_GQ::QuadServices_GQ( vector<Ordinate> const &ordinates,
                                  double const norm,
                                  unsigned const dimension,
                                  unsigned const expansionOrder,
                                  rtt_mesh_element::Geometry const geometry)
    : QuadServices(ordinates, norm, dimension, expansionOrder),
      geometry(geometry),
      n2lk( compute_n2lk( expansionOrder, dimension) ),
      maxExpansionOrder(  max_available_expansion_order(n2lk) ),
      moments( compute_moments(expansionOrder, n2lk) ),
      M( computeM() ),
      D( computeD() )
{ 
    std::vector< unsigned > dimsM;
    dimsM.push_back( getNumMoments() );
    dimsM.push_back( getNumOrdinates() );

    
    std::vector< double > Mm(getM());
    //std::cout << "M matrix is (" << dimsM[0] << " x " << dimsM[1] << std::endl
    //          << "M matrix size is: " << Mm.size() << " (should be: " << dimsM[0]*dimsM[1]  << ")" << std::endl;
    print_matrix( "M", Mm, dimsM );

    std::vector< unsigned > dimsD;
    dimsD.push_back( getNumOrdinates() );
    dimsD.push_back( getNumMoments() );

    std::vector< double > Dd(getD());
    //std::cout << "D matrix is (" << dimsD[0] << " x " << dimsD[1] << std::endl
    //          << "D matrix size is: " << Dd.size() << " (should be: " << dimsD[0]*dimsD[1]  << ")" << std::endl;
    print_matrix( "D", Dd, dimsD );

    //double s=0;
    //for (unsigned i=0; i<dimsD[0]; ++i)
    //    s += Dd[i];
    //std::cout << " sum over D " << s << std::endl;

    Ensure( geometry != rtt_mesh_element::CARTESIAN || D_equals_M_inverse() );
}

//! \brief Return the moment-to-discrete operator.
std::vector< double > QuadServices_GQ::getM_() const
{ 
    unsigned const numOrdinates(this->getNumOrdinates());
    unsigned const numMoments(this->getNumMoments());
    unsigned const numGQMoments(n2lk.size());

    std::vector< double > Mmatrix(numOrdinates*numMoments);
    
    for( unsigned n=0; n<numMoments; ++n )
        for( unsigned m=0; m<numOrdinates; ++m )
            Mmatrix[ n + m*numMoments ] = M[n + m*numGQMoments ];
    
    return Mmatrix; 
}

//! \brief Return the discrete-to-moment operator.
std::vector< double > QuadServices_GQ::getD_() const
{
    unsigned const numOrdinates(this->getNumOrdinates());
    unsigned const numMoments(this->getNumMoments());

    std::vector< double > Dmatrix(numOrdinates*numMoments);

    for( unsigned m=0; m<numOrdinates; ++m )
        for( unsigned n=0; n<numMoments; ++n )
            Dmatrix[ m + n*numOrdinates ] = D[ m + n*numOrdinates ];
    
    return Dmatrix;
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS SPECIFIC TO THE GALERKIN QUADRATURE METHOD
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the moment-to-discrete matrix.
 */

std::vector< double >
QuadServices_GQ::computeM()
{
    using rtt_sf::Ylm;

    std::vector< double > M;

    if (geometry == rtt_mesh_element::CARTESIAN)
    {
        M = computeM(this->getOrdinates(), n2lk, this->getDimension(), this->getNorm());
    }
    else
    {
        // First construct a vector of ordinates without starting directions

        std::vector<Ordinate> const ordinates(this->getOrdinates());
        unsigned const numOrdinates(ordinates.size());

        std::vector<Ordinate> cartesian_ordinates;
        std::vector<unsigned> indexes;
        unsigned count=0;
        for( unsigned m=0; m<numOrdinates; ++m )
        {
            if (ordinates[m].wt() != 0)
            {
                cartesian_ordinates.push_back(ordinates[m]);
                indexes.push_back(count++);
            }
            else
                indexes.push_back(0);
        }

        // And compute the operator for these ordinates only

        std::vector< double > cartesian_M = computeM(cartesian_ordinates, n2lk, this->getDimension(), this->getNorm());

        // Now augment the matrix and store it appropriately

        unsigned const numMoments = n2lk.size();
        M.resize(numMoments*ordinates.size()); 

        for( unsigned n=0; n<numMoments; ++n )
        {
            unsigned const ell ( n2lk[n].first  );
            int      const k   ( n2lk[n].second );

            for( unsigned m=0; m<numOrdinates; ++m )
            {
                if (ordinates[m].wt() != 0)
                {
                    M[ n + m*numMoments ] = cartesian_M[n + indexes[m]*numMoments ];
                }
                else
                {
                    if (this->getDimension() == 1)
                    {
                        double mu ( ordinates[m].mu() );
                        M[ n + m*numMoments ] = Ylm( ell, k, mu, 0.0, this->getNorm());
                    }
                    else
                    {
                        double mu ( ordinates[m].mu() );
                        double eta( ordinates[m].eta() );
                        double xi(  ordinates[m].xi() );
                        
                        double phi( compute_azimuthalAngle(mu, eta, xi) );
                        M[ n + m*numMoments ] = Ylm( ell, k, xi, phi, this->getNorm());
                    }
                }
            }
        }

    }
    
    return M;
}

std::vector< double >
QuadServices_GQ::computeM(std::vector<Ordinate> const &ordinates,
                          std::vector< lk_index > const &n2lk,
                          unsigned const dim,
                          double const sumwt)
{
    using rtt_sf::Ylm;

    unsigned const numOrdinates = ordinates.size();
    unsigned const numMoments = n2lk.size();

    // resize the M matrix.
    std::vector< double > Mmatrix( numMoments*numOrdinates, -9999.0 );

    for( unsigned n=0; n<numMoments; ++n )
    {
        for( unsigned m=0; m<numOrdinates; ++m )
        {
            unsigned const ell ( n2lk[n].first  );
            int      const k   ( n2lk[n].second ); 
        
            if( dim == 1 ) // 1D mesh, 1D quadrature
            { 
                double mu ( ordinates[m].mu() );
                Mmatrix[ n + m*numMoments ] = Ylm( ell, k, mu, 0.0, sumwt );
            }
            else 
            {
                // It is important to remember that the positive mu axis points to the
                // left and the positive eta axis points up, when the unit sphere is
                // projected on the plane of the mu- and eta-axis. In this case, phi is
                // measured from the mu-axis counterclockwise.
                //
                // This accounts for the fact that the aziumuthal angle is discretized
                // on levels of the xi-axis, making the computation of the azimuthal angle
                // here consistent with the discretization by using the eta and mu
                // ordinates to define phi.
                
                double mu ( ordinates[m].mu() );
                double eta( ordinates[m].eta() );
                double xi(  ordinates[m].xi() );

                double phi( compute_azimuthalAngle(mu, eta, xi) );
                Mmatrix[ n + m*numMoments ] = Ylm( ell, k, xi, phi, sumwt );
            }
        } // n: end moment loop
    } // m: end ordinate loop

    return Mmatrix;
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
std::vector< double > QuadServices_GQ::computeD()
{
    std::vector<double> D;

    if (geometry == rtt_mesh_element::CARTESIAN)
    {
        D = computeD(this->getOrdinates(), n2lk, M, this->getDimension(), this->getNorm());
    }
    else
    {
        // First construct a vector of ordinates without starting directions

        std::vector<Ordinate> const ordinates(this->getOrdinates());
        unsigned const numOrdinates(ordinates.size());

        std::vector<Ordinate> cartesian_ordinates;
        std::vector<unsigned> indexes;
        unsigned count=0;
        for( unsigned m=0; m<numOrdinates; ++m )
        {
            if (ordinates[m].wt() != 0)
            {
                cartesian_ordinates.push_back(ordinates[m]);
                indexes.push_back(count++);
            }
            else
                indexes.push_back(0);
        }

        // And compute the operators for these ordinates only

        std::vector< double > cartesian_M = computeM(cartesian_ordinates, n2lk, this->getDimension(), this->getNorm());
        std::vector< double > cartesian_D = computeD(cartesian_ordinates, n2lk, cartesian_M, this->getDimension(), this->getNorm());

        // Now augment the matrix and store it appropriately

        unsigned const numMoments = n2lk.size();
        D.resize(numMoments*ordinates.size()); 
        
        for( unsigned m=0; m<numOrdinates; ++m )
        {
            for( unsigned n=0; n<numMoments; ++n )
            {
                if (ordinates[m].wt() != 0)
                {
                    D[ m + n*numOrdinates ] = cartesian_D[indexes[m] + n*numOrdinates];
                }
                else
                {
                    D[ m + n*numOrdinates ] = 0;
                }
            }
        }

    }
    

    return D;
}

std::vector< double >
QuadServices_GQ::computeD(std::vector<Ordinate> const &ordinates,
                          std::vector< lk_index > const &n2lk,
                          std::vector<double> const &mM,
                          unsigned const,
                          double const)
{
    int n = n2lk.size();
    int m = ordinates.size();

    Require( n == m );

    std::vector< double > M(mM);
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
    Remember(int result = ) gsl_linalg_LU_decomp( &gsl_M.matrix, p, &signum );
    Check( result == 0 );
    Check( diagonal_not_zero( M, n, m ) );

    // Compute the inverse of the matrix LU from its LU decomposition (LU,p),
    // storing the results in the matrix Dmatrix.  The inverse is computed by
    // solving the system (LU) x = b for each column of the identity matrix.

    Remember(result = ) gsl_linalg_LU_invert( &gsl_M.matrix, p, &gsl_D.matrix );

    Check( result == 0 );

    // Free the space reserved for the permutation matrix.
    gsl_permutation_free( p );

    return D;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 *
 * This function computes the mapping as specified by Morel in "A Hybrid
 * Collocation-Galerkin-Sn Method for Solving the Boltzmann Transport
 * Equation." 
 */

std::vector< QuadServices::lk_index >
QuadServices_GQ::compute_n2lk( unsigned const expansionOrder,
                              unsigned const dim)
{
    unsigned const L( expansionOrder+1 );

    if( dim == 3 )
    {
        return compute_n2lk_3D(L);
    }

    Check( dim < 3 );
    if( dim == 2 )
    {
        return compute_n2lk_2D(L);
    }
    
    Check( dim == 1 );
    return this->compute_n2lk_1D(L);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index >
QuadServices_GQ::compute_n2lk_3D(unsigned const L)
{
    std::vector< lk_index > result;

    // Choose: l= 0, ..., L-1, k = -l, ..., l
    for( unsigned ell=0; ell< L; ++ell )
	for( int k = -ell; k <= static_cast<int>(ell); ++k )
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
std::vector< QuadServices::lk_index >
QuadServices_GQ::compute_n2lk_2D( unsigned const L )
{
    std::vector< lk_index > result;
    
    // Choose: l= 0, ..., N-1, k = 0, ..., l
    for( unsigned ell=0; ell<L; ++ell )
	for( int k=0; k<=static_cast<int>(ell); ++k )
	    result.push_back( lk_index(ell,k) );

    // Add ell=N and k>0, k odd
    int ell( L );
    for( int k=1; k<=ell; k+=2 )
        result.push_back( lk_index(ell,k) );
    
    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Multiply M and D and compare the result to the identity matrix.
 * \return true if M = D^(-1), otherwise false.
 */
bool QuadServices_GQ::D_equals_M_inverse() const
{
    using rtt_dsxx::soft_equiv;

    unsigned n( n2lk.size() );
    int m( this->getOrdinates().size() );

    // create non-const versions of M and D.
    std::vector< double > Marray( M );
    std::vector< double > Darray( D );
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


} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of QuadServices_GQ.cc
//---------------------------------------------------------------------------//

