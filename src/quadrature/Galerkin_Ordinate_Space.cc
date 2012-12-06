//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Galerkin_Ordinate_Space.cc
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Define methods of class Galerkin_Ordinate_Space
 * \note   Copyright (C) 2006-2012 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------------------//
// $Id: Galerkin_Ordinate_Space.cc 6855 2012-11-06 16:39:27Z kellyt $
//---------------------------------------------------------------------------------------//

// Vendor software
#include <gsl/gsl_linalg.h>
// #include <gsl/gsl_blas.h>
// #include <gsl/gsl_sf_legendre.h>

#include "Galerkin_Ordinate_Space.hh"

#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"

using namespace rtt_units;

namespace rtt_quadrature
{

//---------------------------------------------------------------------------------------//
vector< Moment >
Galerkin_Ordinate_Space::compute_n2lk_1D_( Quadrature_Class,
                                           unsigned const L )
{
    vector< Moment > result;

    // Choose: l= 0, ..., L-1, k = 0
    int k(0); // k is always zero for 1D.

    for( unsigned ell=0; ell<L; ++ell )
        result.push_back( Moment(ell,k) );

    return result;
}

//---------------------------------------------------------------------------------------//
vector< Moment >
Galerkin_Ordinate_Space::compute_n2lk_1Da_( Quadrature_Class,
                                            unsigned const L )
{
    std::vector< Moment > result;
    
    // Choose: l= 0, ..., N, k = 0, ..., l to eliminate moments even in xi
    // which are identically zero by symmetry.
    for( int ell=0; ell<static_cast<int>(L); ++ell )
        for( int k=0; k<=ell; ++k )
            if ((ell-k)%2 == 0)
                // Eliminate moments even in eta which are identically zero by
                // symmetry.
                result.push_back( Moment(ell,k) );

    // Add ell=N and k>0, k odd
    int ell = L;
    for( int k=1; k<=ell; k+=2 )
        if ((ell-k)%2 == 0)
            result.push_back( Moment(ell,k) );
    
    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
vector< Moment >
Galerkin_Ordinate_Space::compute_n2lk_2D_( Quadrature_Class,
                                           unsigned const L )
{
    std::vector< Moment > result;
    
    // Choose: l= 0, ..., N-1, k = 0, ..., l  to eliminate moments even in xi,
    // which are identically zero by symmetry.
    for( unsigned ell=0; ell<L; ++ell )
	for( int k=0; k<=static_cast<int>(ell); ++k )
	    result.push_back( Moment(ell,k) );

    // Add ell=N and k>0, k odd
    int ell = L;
    for( int k=1; k<=ell; k+=2 )
        result.push_back( Moment(ell,k) );
    
    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
vector< Moment >
Galerkin_Ordinate_Space::compute_n2lk_3D_(Quadrature_Class,
                                          unsigned const L)
{
    vector< Moment > result;

    // Choose: l= 0, ..., L-1, k = -l, ..., l
    for( unsigned ell=0; ell<L; ++ell )
	for( int k = -ell; k <= static_cast<int>(ell); ++k )
	    result.push_back( Moment(ell,k) );

    // Add ell=L and k<0
    {
	unsigned ell( L );
	for( int k(-1*static_cast<int>(ell)); k<0; ++k )
	    result.push_back( Moment(ell,k) );
    }

    // Add ell=L, k>0, k odd
    {
	int ell( L );
	for( int k=1; k<=ell; k+=2 )
	    result.push_back( Moment(ell,k) );
    }

    // Add ell=L+1 and k<0, k even
    {
	unsigned ell( L+1 );
	for( int k(-1*static_cast<int>(ell)+1); k<0; k+=2 )
	    result.push_back( Moment(ell,k) );
    }

    return result;
}

//---------------------------------------------------------------------------------------//
void
Galerkin_Ordinate_Space::compute_M()
{
    using rtt_sf::Ylm;

    std::vector< double > M;

    rtt_mesh_element::Geometry const geometry = this->geometry();
    vector<Moment> const &n2lk = this->moments();
    
    if (geometry == rtt_mesh_element::CARTESIAN)
    {
        M = compute_M_GQ(this->ordinates(),
                         n2lk,
                         this->dimension(),
                         this->norm());
    }
    else
    {
        // First construct a vector of ordinates without starting directions

        std::vector<Ordinate> const ordinates(this->ordinates());
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

        std::vector< double > cartesian_M =
            compute_M_GQ(cartesian_ordinates, n2lk, this->dimension(), this->norm());

        // Now augment the matrix and store it appropriately

        unsigned const numMoments = n2lk.size();
        M.resize(numMoments*ordinates.size()); 

        for( unsigned n=0; n<numMoments; ++n )
        {
            unsigned const ell ( n2lk[n].L() );
            int      const k   ( n2lk[n].M() );

            for( unsigned m=0; m<numOrdinates; ++m )
            {
                if (ordinates[m].wt() != 0)
                {
                    M[ n + m*numMoments ] = cartesian_M[n + indexes[m]*numMoments ];
                }
                else
                {
                    if (this->dimension() == 1)
                    {
                        double mu ( ordinates[m].mu() );
                        M[ n + m*numMoments ] = Ylm( ell, k, mu, 0.0, this->norm());
                    }
                    else
                    {
                        double mu ( ordinates[m].mu() );
                        double eta( ordinates[m].eta() );
                        double xi(  ordinates[m].xi() );
                        
                        double phi( compute_azimuthalAngle(mu, xi, eta) );
                        M[ n + m*numMoments ] = Ylm( ell, k, eta, phi, this->norm());
                    }
                }
            }
        }

    }
    
    M_.swap( M );
}

//---------------------------------------------------------------------------------------//

vector< double >
Galerkin_Ordinate_Space::compute_M_GQ(vector<Ordinate> const &ordinates,
                             vector< Moment > const &n2lk,
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
            unsigned const ell ( n2lk[n].L() );
            int      const k   ( n2lk[n].M() ); 
        
            if( dim == 1 ) // 1D mesh, 1D quadrature
            { 
                double mu ( ordinates[m].mu() );
                Mmatrix[ n + m*numMoments ] = Ylm( ell, k, mu, 0.0, sumwt );
            }
            else 
            {
                double mu ( ordinates[m].mu() );
                double eta( ordinates[m].eta() );
                double xi ( ordinates[m].xi() );

                double phi( compute_azimuthalAngle(mu, xi, eta) );
                Mmatrix[ n + m*numMoments ] = Ylm( ell, k, eta, phi, sumwt );
            }
        } // n: end moment loop
    } // m: end ordinate loop

    return Mmatrix;
}

//---------------------------------------------------------------------------------------//

vector<double>
Galerkin_Ordinate_Space::compute_D_GQ(vector<Ordinate> const &ordinates,
                             vector< Moment > const &n2lk,
                             vector<double> const &mM,
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
    // Check( diagonal_not_zero( M, n, m ) );

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
 * \brief Compute the discrete-to-moment matrix. 
 *
 * Computes \f$ \mathbf{D} \equiv \mathbf{M}^{-1} \f$.  This private function
 * is called by the constuctor.
 *
 * Normally, M will not be square because we only have isotropic scatting.
 * For isotropic scattering M will be (numOrdinates x 1 moment).  We will use the
 * Moore-Penrose Pseudo-Inverse Matrix, \f$ D = (M^T * M)^-1 * M^T.\f$
 */
void Galerkin_Ordinate_Space::compute_D()
{
    std::vector<double> D;

    rtt_mesh_element::Geometry const geometry = this->geometry();
    vector<Moment> const &n2lk = this->moments();

    if (geometry == rtt_mesh_element::CARTESIAN)
    {
        D = compute_D_GQ(this->ordinates(), n2lk, M_, this->dimension(), this->norm());
    }
    else
    {
        // First construct a vector of ordinates without starting directions

        std::vector<Ordinate> const ordinates(this->ordinates());
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
        unsigned numCartOrdinates = cartesian_ordinates.size();

        // And compute the operators for these ordinates only

        std::vector< double > cartesian_M =
            compute_M_GQ(cartesian_ordinates, n2lk, this->dimension(), this->norm());
        
        std::vector< double > cartesian_D =
            compute_D_GQ(cartesian_ordinates,
                         n2lk, cartesian_M,
                         this->dimension(),
                         this->norm());

        // Now augment the matrix and store it appropriately

        unsigned const numMoments = n2lk.size();
        D.resize(numMoments*numOrdinates); 
        
        for( unsigned m=0; m<numOrdinates; ++m )
        {
            for( unsigned n=0; n<numMoments; ++n )
            {
                if (ordinates[m].wt() != 0)
                {
                    D[ m + n*numOrdinates ] = cartesian_D[indexes[m] + n*numCartOrdinates];
                }
                else
                {
                    D[ m + n*numOrdinates ] = 0;
                }
            }
        }

    }

    D_.swap(D);
}


//---------------------------------------------------------------------------------------//
/*!
 *
 * \param dimension Dimension of the physical problem space (1, 2, or 3)
 *
 * \param geometry Geometry of the physical problem space (spherical,
 * axisymmetric, Cartesian)
 *
 * \param ordinates Set of ordinate directions
 *
 * \param quadrature_class Class of the quadrature used to generate the
 * ordinate set. At presente, only TRIANGLE_QUADRATURE is supported.
 *
 * \param sn_order Order of the quadrature. This is equal to the number of
 * levels for triangular and square quadratures.
 *
 * \param expansion_order Expansion order of the desired scattering moment
 * space.
 *
 * \param extra_starting_directions Add extra directions to each level set. In most
 * geometries, an additional ordinate is added that is opposite in direction
 * to the starting direction. This is used to implement reflection exactly in
 * curvilinear coordinates. In 1D spherical, that means an additional angle is
 * added at mu=1. In axisymmetric, that means additional angles are added that
 * are oriented opposite to the incoming starting direction on each level.
 *
 * \param ordering Ordering into which to sort the ordinates.
 */

Galerkin_Ordinate_Space::Galerkin_Ordinate_Space( unsigned const  dimension,
                                                  Geometry const  geometry,
                                                  vector<Ordinate> const &ordinates,
                                                  Quadrature_Class quadrature_class,
                                                  unsigned sn_order,
                                                  unsigned const  expansion_order,
                                                  bool const  extra_starting_directions,
                                                  Ordering const ordering)
    : Ordinate_Space(dimension,
                         geometry,
                         ordinates,
                         expansion_order,
                         extra_starting_directions,
                         ordering)
{
    Require(dimension>0 && dimension<4);
    Require(geometry!=rtt_mesh_element::END_GEOMETRY);
    Require(expansion_order<=sn_order);  // May be relaxed in the future
    Require(quadrature_class == TRIANGLE_QUADRATURE || dimension==1);
    Require(sn_order>0 && sn_order%2==0);

    compute_moments_(quadrature_class,
                     sn_order);

    Ensure(check_class_invariants());
}


//---------------------------------------------------------------------------------------//
bool Galerkin_Ordinate_Space::check_class_invariants() const
{
    return
        D_.size() == ordinates().size() * moments().size() &&
        M_.size() == ordinates().size() * moments().size();
}

//---------------------------------------------------------------------------------------//
QIM Galerkin_Ordinate_Space::quadrature_interpolation_model() const
{
    return GQ;
}


//---------------------------------------------------------------------------------------//
/*!
 * In the future, this function will allow the client to specify the maximum
 * order to include, but for the moment, we include all full orders, leaving
 * out any Galerkin augments.
 */

vector<double> Galerkin_Ordinate_Space::D() const
{
    unsigned const number_of_ordinates = ordinates().size();
    unsigned const number_of_moments = this->number_of_moments();

    vector<double> Result(number_of_ordinates * number_of_moments);

    for (unsigned a=0; a<number_of_ordinates; ++a)
    {
        for (unsigned m=0; m<number_of_moments; ++m)
        {
            Result[a + number_of_ordinates*m] = D_[a + number_of_ordinates*m];
        }
    }
    return Result;
}
//---------------------------------------------------------------------------------------//
/*!
 * In the future, this function will allow the client to specify the maximum
 * order to include, but for the moment, we include all full orders, leaving
 * out any Galerkin augments.
 */
vector<double> Galerkin_Ordinate_Space::M() const
{
    unsigned const number_of_ordinates = ordinates().size();
    unsigned const number_of_moments = this->number_of_moments();
    unsigned const total_number_of_moments = moments().size();

    vector<double> Result(number_of_ordinates * number_of_moments);

    for (unsigned a=0; a<number_of_ordinates; ++a)
    {
        for (unsigned m=0; m<number_of_moments; ++m)
        {
            Result[m + number_of_moments*a] = M_[m + total_number_of_moments*a];
        }
    }
    return Result;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                 end of Galerkin_Ordinate_Space.cc
//---------------------------------------------------------------------------------------//
