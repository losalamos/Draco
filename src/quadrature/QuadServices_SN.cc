//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices_SN.cc
 * \author Kelly Thompson
 * \date   Mon Nov  8 11:17:12 2004
 * \brief  Provide Moment-to-Discrete and Discrete-to-Moment operators.
 * \note   © Copyright 2006 LANSLLC All rights reserved. 
 */
//---------------------------------------------------------------------------//
// $Id: QuadServices_SN.cc 6499 2012-03-15 20:19:33Z kgbudge $
//---------------------------------------------------------------------------//

#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>

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

#include "QuadServices_SN.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor builds square D and M operators using Morel's
 * Galerkin-Sn heuristic. 
 * \param spQuad_ a smart pointer to a Quadrature object.
 * \post \f$ \mathbf{D} = \mathbf{M}^{-1} \f$.
 */
QuadServices_SN::QuadServices_SN( vector<Ordinate> const &ordinates,
                                  double const norm,
                                  unsigned const dimension,
                                  unsigned const expansionOrder,
                                  rtt_mesh_element::Geometry const geometry)
    : QuadServices(ordinates, norm, dimension, expansionOrder),
      geometry(geometry),
      n2lk( compute_n2lk( expansionOrder, dimension) ),
      maxExpansionOrder(  max_available_expansion_order(n2lk) ),
      moments( compute_moments(maxExpansionOrder, n2lk) ),
      M( computeM() ),
      D( computeD() )
{ 
    std::vector< unsigned > dimsM;
    dimsM.push_back( getNumMoments() );
    dimsM.push_back( getNumOrdinates() );
    print_matrix( "M", M, dimsM );

    std::vector< unsigned > dimsD;
    dimsD.push_back( getNumOrdinates() );
    dimsD.push_back( getNumMoments() );
    print_matrix( "D", D, dimsD );

    //double s=0;
    //for (unsigned i=0; i<dimsD[0]; ++i)
    //    s += D[i];
    //std::cout << " sum over D " << s << std::endl;
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS SPECIFIC TO THE SN INTERPOLATION METHOD
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the standard SN moment-to-discrete matrix.
 */

std::vector< double >
QuadServices_SN::computeM() 
{
    return computeM(this->getOrdinates(), n2lk, this->getDimension(), this->getNorm());
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Create the M array (moment-to-discrete matrix).
 * \return The moment-to-discrete matrix.
 *
 * This static member function may be used by clients who provide their own
 * ordinate sets.
 * 
 * The moment-to-discrete matrix will be num_moments by num_ordinates in size.
 */

std::vector< double >
QuadServices_SN::computeM(std::vector<Ordinate> const &ordinates,
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
        
            if( dim == 1 && geometry != rtt_mesh_element::AXISYMMETRIC) // 1D mesh, 1D quadrature
            { 
                double mu ( ordinates[m].mu() );
                Mmatrix[ n + m*numMoments ] = Ylm( ell, k, mu, 0.0, sumwt );
            }
            else if ( dim == 1 ) // 1D mesh, 2D quadrature
            {
                double mu ( ordinates[m].mu() );
                double eta( ordinates[m].eta() );
                double xi(  ordinates[m].xi() );

                double phi( compute_azimuthalAngle(mu, eta, xi) );
                if ((ell-k)%2 == 0)
                    Mmatrix[ n + m*numMoments ] = Ylm( ell, k, xi, phi, sumwt );
                else
                    Mmatrix[ n + m*numMoments ] = 0;
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
 * \brief Compute the standard SN discrete-to-moment matrix.
 *
 */
//---------------------------------------------------------------------------//
std::vector< double >
QuadServices_SN::computeD() 
{
    vector<double> D;

    //D = computeD(this->getOrdinates(), n2lk, this->getDimension(), this->getNorm());
    D = computeD(this->getOrdinates(), n2lk, getM());

    return D;
}

std::vector< double >
QuadServices_SN::computeD(std::vector<Ordinate> const &ordinates,
                          std::vector< lk_index > const &n2lk,
                          unsigned const dim,
                          double const sumwt)
{
    using rtt_sf::Ylm;

    unsigned const numOrdinates = ordinates.size();
    unsigned const numMoments = n2lk.size();
    
    std::vector< double > D( numOrdinates*numMoments );
    
    for( unsigned m=0; m<numOrdinates; ++m )
    {
        double mu( ordinates[m].mu() );
        double wt( ordinates[m].wt() );

        for( unsigned n=0; n<numMoments; ++n )
        {
            unsigned const ell ( n2lk[n].first  );
            int      const k   ( n2lk[n].second );  

            if( dim == 1 ) // 1D mesh, 1D quadrature
            {   // for 1D, phi = 0 and k==0
                D[ m + n*numOrdinates ] = wt * Ylm( ell, k, mu, 0.0, sumwt ); 
            }
            else 
            { 
                double eta( ordinates[m].eta() );
                double xi ( ordinates[m].xi() );

                double phi( compute_azimuthalAngle( mu, eta, xi) );
                D[ m + n*numOrdinates ] = wt * Ylm( ell, k, xi, phi, sumwt );
            }
            
        } // n: end moment loop
    } // m: end ordinate loop

    return D;
}

//---------------------------------------------------------------------------//

std::vector< double >
QuadServices_SN::computeD(std::vector<Ordinate> const &ordinates,
                          std::vector< lk_index > const &n2lk,
                          std::vector<double> const &Mm)
{
    unsigned const numOrdinates = ordinates.size();
    unsigned const numMoments = n2lk.size();

    // ---------------------------------------------------
    // Create diagonal matrix of quadrature weights
    // ---------------------------------------------------

    gsl_matrix *gsl_W = gsl_matrix_alloc(numOrdinates, numOrdinates);
    gsl_matrix_set_identity(gsl_W);

    for( unsigned m=0; m<numOrdinates; ++m )
        gsl_matrix_set(gsl_W, m, m, ordinates[m].wt());

    // ---------------------------------------------------
    // Create the discrete-to-moment matrix 
    // ---------------------------------------------------

    std::vector< double > M( Mm );
    gsl_matrix_view gsl_M = gsl_matrix_view_array( &M[0], numOrdinates, numMoments );

    std::vector< double > D( numMoments*numOrdinates );  // rows x cols
    gsl_matrix_view gsl_D = gsl_matrix_view_array( &D[0], numMoments, numOrdinates);

    unsigned ierr = gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, &gsl_M.matrix, gsl_W, 0.0, &gsl_D.matrix);
    Insist(!ierr, "GSL blas interface error");

    gsl_matrix_free( gsl_W);

    return D;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 *
 * This function computes the mapping for the standard SN method.
 */
std::vector< QuadServices::lk_index >
QuadServices_SN::compute_n2lk( unsigned const expansionOrder,
                               unsigned const dim)
{
    unsigned const L( expansionOrder+1 );

    if( dim == 3 )
    {
        return compute_n2lk_3D(L);
    }
    else if( dim == 2 || geometry == rtt_mesh_element::AXISYMMETRIC)
    {
        return compute_n2lk_2D(L);
    }
    else
    {
        Check( dim == 1 );
        return this->compute_n2lk_1D(L);
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
std::vector< QuadServices::lk_index >
QuadServices_SN::compute_n2lk_3D( unsigned const L )
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
std::vector< QuadServices::lk_index >
QuadServices_SN::compute_n2lk_2D( unsigned const L )
{
    std::vector< lk_index > result;
    
    // Choose: l= 0, ..., N, k = 0, ..., l
    for( int ell=0; ell<static_cast<int>(L); ++ell )
	for( int k=0; k<=ell; ++k )
	    result.push_back( lk_index(ell,k) );

    return result;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of QuadServices_SN.cc
//---------------------------------------------------------------------------//

