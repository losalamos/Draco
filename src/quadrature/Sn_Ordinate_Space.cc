//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Sn_Ordinate_Space.cc
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Define methods of class Sn_Ordinate_Space
 * \note   Copyright (C) 2007-2013 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------------------//
// $Id: Sn_Ordinate_Space.cc 6855 2012-11-06 16:39:27Z kellyt $
//---------------------------------------------------------------------------------------//

#include "Sn_Ordinate_Space.hh"

#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"

// Vendor software
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <iomanip>

/*
static
void print_matrix( std::string const & matrix_name,
                   std::vector<double> const & x,
                   std::vector<unsigned> const & dims )
{
    using std::cout;
    using std::endl;
    using std::string;

    Require( dims[0]*dims[1] == x.size() );

    unsigned pad_len( matrix_name.length()+2 );
    string padding( pad_len, ' ' );
    cout << matrix_name << " =";
    // row
    for( unsigned i=0; i<dims[1]; ++i )
    {
        if( i != 0 ) cout << padding;

        cout << "{ ";

        for( unsigned j=0; j<dims[0]-1; ++j )
            cout << std::setprecision(10) << x[j+dims[0]*i] << ", ";

        cout << std::setprecision(10) << x[dims[0]-1+dims[0]*i] << " }." << endl;
    }
    cout << endl;
    return;
}
*/

namespace rtt_quadrature
{

//---------------------------------------------------------------------------------------//
vector< Moment >
Sn_Ordinate_Space::compute_n2lk_1D_( Quadrature_Class,
                                     unsigned /*N*/)
{
    vector< Moment > result;

    unsigned const L = expansion_order();

    // Choose: l= 0, ..., L-1, k = 0
    int k(0); // k is always zero for 1D.

    for( unsigned ell=0; ell<=L; ++ell )
        result.push_back( Moment(ell,k) );

    return result;
}

//---------------------------------------------------------------------------------------//
vector< Moment >
Sn_Ordinate_Space::compute_n2lk_1Da_( Quadrature_Class,
                                      unsigned /*N*/)
{
    vector< Moment > result;

    unsigned const L = expansion_order();

    // Choose: l= 0, ..., N, k = 0, ..., l
    for( int ell=0; ell<=static_cast<int>(L); ++ell )
        for( int k=0; k<=ell; ++k )
            if ((ell+k)%2 == 0)
                result.push_back( Moment(ell,k) );

    return result;
}

//---------------------------------------------------------------------------------------//
vector< Moment >
Sn_Ordinate_Space::compute_n2lk_2D_( Quadrature_Class,
                                     unsigned /*N*/)
{
    vector< Moment > result;

    unsigned const L = expansion_order();

    // Choose: l= 0, ..., N, k = 0, ..., l
    for( int ell=0; ell<=static_cast<int>(L); ++ell )
        for( int k=0; k<=ell; ++k )
            result.push_back( Moment(ell,k) );

    return result;
}

//---------------------------------------------------------------------------------------//
vector< Moment >
Sn_Ordinate_Space::compute_n2lk_2Da_( Quadrature_Class quadrature_class,
                                      unsigned N)
{
    // This is the same as the X-Y moment mapping
    vector< Moment > result=compute_n2lk_2D_(quadrature_class, N);
    return result;
}

//---------------------------------------------------------------------------------------//
vector< Moment >
Sn_Ordinate_Space::compute_n2lk_3D_( Quadrature_Class,
                                     unsigned /*N*/)
{
    vector< Moment > result;

    unsigned const L = expansion_order();

    // Choose: l= 0, ..., L, k = -l, ..., l
    for( int ell=0; ell<=static_cast<int>(L); ++ell )
        for( int k(-static_cast<int>(ell)); k <= ell; ++k )
            result.push_back( Moment(ell,k) );

    return result;
}


//---------------------------------------------------------------------------------------//
void
Sn_Ordinate_Space::compute_M()
{
    using rtt_sf::Ylm;

    vector<Ordinate> const &ordinates = this->ordinates();
    unsigned const numOrdinates = ordinates.size();

    vector<Moment> const &n2lk = this->moments();
    unsigned const numMoments = n2lk.size();

    unsigned const dim = dimension();
    Geometry const geometry(this->geometry());
    double const sumwt(norm());

    // resize the M matrix.
    M_.resize( numMoments*numOrdinates );
    
//    double polar, azimuthal;
    for( unsigned n=0; n<numMoments; ++n )
    {
        for( unsigned m=0; m<numOrdinates; ++m )
        {
            unsigned const ell ( n2lk[n].L() );
            int      const k   ( n2lk[n].M() ); 
        
            if( dim == 1 && geometry != rtt_mesh_element::AXISYMMETRIC) // 1D mesh, 1D quadrature 
            { 
                double mu ( ordinates[m].mu() );
                M_[ n + m*numMoments ] = Ylm( ell, k, mu, 0.0, sumwt );

//                polar = mu;
//                azimuthal = 0.0;
            }
            else 
            {
                double mu ( ordinates[m].mu() );
                double eta( ordinates[m].eta() );
                double xi(  ordinates[m].xi() );
                
                if (geometry == rtt_mesh_element::AXISYMMETRIC)
                {
                    // R-Z coordinate system
                    //
                    // It is important to remember here that the positive mu axis points to the
                    // left and the positive eta axis points up, when the unit sphere is
                    // projected on the plane of the mu- and eta-axis in R-Z. In this case, phi is
                    // measured from the mu-axis counterclockwise.
                    //
                    // This accounts for the fact that the aziumuthal angle is discretized
                    // on levels of the xi-axis, making the computation of the azimuthal angle
                    // here consistent with the discretization by using the eta and mu
                    // ordinates to define phi.

                    double phi( compute_azimuthalAngle(mu, xi) );
                    M_[ n + m*numMoments ] = Ylm( ell, k, eta, phi, sumwt );

//                    polar = eta;
//                    azimuthal = phi;
                }
                else if (geometry == rtt_mesh_element::CARTESIAN)
                {
                    // X-Y coordinate system
                    //
                    // Note that we choose the same moments and spherical
                    // harmonics as for R-Z in this case, unlike the Galerkin
                    // method.
                    //
                    // This is because we choose the "front" of the
                    // hemisphere, here, so that the spherical harmoincs
                    // chosen are even in the azimuthal angle (symmetry from
                    // front to back) and not even in the polar angle.
                    // Thus, in this case, the polar angle is measured from the
                    // eta-axis [0, Pi], and the azimuthal angle is measured
                    // from the mu-axis [0,Pi].
                    //
                    // In contrast, the Galerkin methods chooses the "top"
                    // hemisphere, and projects down onto the x-y plane.
                    // Hence the polar angle in that case is xi and extends from
                    // [0,Pi/2] while the azimuthal angle is on [0, 2 Pi].
                    // Therefore, in that case, the spherical harmonics must
                    // be those that are even in the polar angle.
                    // That may be determined by considering the even-ness
                    // of the associated legendre polynomials.

                    double phi( compute_azimuthalAngle(mu, xi) );
                    M_[ n + m*numMoments ] = Ylm( ell, k, eta, phi, sumwt );

//                    polar = mu;
//                    azimuthal = phi;
                }
            }        
/*
            if (n == 0)
                    std::cout << "   " << m
                              << "   " << ordinates[m].mu() 
                              << "   " << ordinates[m].eta() 
                              << "   " << ordinates[m].xi() 
                              << "   " << ordinates[m].wt() 
                              << "   " << polar
                              << "   " << azimuthal*180.0/3.141592653589793238462643383279
                              << std::endl;
*/
        } // ordinate loop
    } // moment loop

}

//---------------------------------------------------------------------------------------//
/*! This computation requires that the moment-to-discrete matrix M
 *  is already created.
 */

void
Sn_Ordinate_Space::compute_D()
{

    Insist(!M_.empty(), "The SN ordinate space computation for D requires that M be available.");

    vector<Ordinate> const &ordinates = this->ordinates();
    unsigned const numOrdinates = ordinates.size();
    vector<Moment> const &n2lk = this->moments();
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
    
    std::vector< double > M( M_ );
    gsl_matrix_view gsl_M = gsl_matrix_view_array( &M[0], numOrdinates, numMoments );
    
    std::vector< double > D( numMoments*numOrdinates );  // rows x cols
    gsl_matrix_view gsl_D = gsl_matrix_view_array( &D[0], numMoments, numOrdinates);
    
    unsigned ierr = gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, &gsl_M.matrix, gsl_W, 0.0, &gsl_D.matrix);
    Insist(!ierr, "GSL blas interface error");
    
    gsl_matrix_free( gsl_W);
    
    D_.swap(D);
}


//---------------------------------------------------------------------------------------//
/*!
 *
 * The computation of the tau and alpha coefficients is described by Morel in
 * various technical notes on the treatment of the angle derivatives in the
 * streaming operator.
 *
 * \param dimension Dimension of the physical problem space (1, 2, or 3)
 *
 * \param geometry Geometry of the physical problem space (spherical,
 * axisymmetric, Cartesian)
 *
 * \param ordinates Set of ordinate directions
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

Sn_Ordinate_Space::Sn_Ordinate_Space( unsigned const  dimension,
                                Geometry const  geometry,
                                vector<Ordinate> const &ordinates,
                                unsigned const  expansion_order,
                                bool const  extra_starting_directions,
                                Ordering const ordering)
    : Ordinate_Space(dimension,
                     geometry,
                     ordinates,
                     expansion_order,
                     extra_starting_directions,
                     ordering),
      D_(),
      M_()
{
    Require(dimension>0 && dimension<4);
    Require(geometry!=rtt_mesh_element::END_GEOMETRY);

    compute_moments_(END_QUADRATURE,   // not used by Sn
                     expansion_order); // also not actually used

    // compute the operators; MUST be called in this order
    compute_M();
    compute_D();

/*
    unsigned const numOrdinates(this->ordinates().size());
    unsigned const numMoments(this->moments().size());

    for( unsigned n=0; n<numMoments; ++n )
    {
        unsigned const ell ( moments()[n].L() );
        int      const k   ( moments()[n].M() ); 
        
        std::cout << " moment " << n
                  << "     l = " << ell << " k = " << k
                  << std::endl;
    }

    std::vector< unsigned > dimsM;
    dimsM.push_back( numMoments );
    dimsM.push_back( numOrdinates );
    print_matrix( "M", M_, dimsM );
        
    std::vector< unsigned > dimsD;
    dimsD.push_back( numOrdinates );
    dimsD.push_back( numMoments );
    print_matrix( "D", D_, dimsD );

    std::cout << " Ordinate Set (may differ from quadrature) " << std::endl; 
    for (unsigned i=0; i<numOrdinates; ++i)
        std::cout << "   " << i
                  << "   " << ordinates[i].mu() 
                  << "   " << ordinates[i].eta() 
                  << "   " << ordinates[i].xi() 
                  << "   " << ordinates[i].wt() 
                  << std::endl;
*/

    Ensure(check_class_invariants());
}

//---------------------------------------------------------------------------------------//
bool Sn_Ordinate_Space::check_class_invariants() const
{
    return
        D_.size() == ordinates().size() * moments().size() &&
        M_.size() == ordinates().size() * moments().size();
}

//---------------------------------------------------------------------------------------//
QIM Sn_Ordinate_Space::quadrature_interpolation_model() const
{
    return SN;
}

//---------------------------------------------------------------------------------------//
/*!
 * In the future, this function will allow the client to specify the maximum
 * order to include, but for the moment, we include all orders.
 */

vector<double> Sn_Ordinate_Space::D() const
{
    return D_;
}
//---------------------------------------------------------------------------------------//
/*!
 * In the future, this function will allow the client to specify the maximum
 * order to include, but for the moment, we include all orders.
 */
vector<double> Sn_Ordinate_Space::M() const
{
    return M_;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                 end of Sn_Ordinate_Space.cc
//---------------------------------------------------------------------------------------//
