//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide Moment-to-Discrete and Discrete-to-Moment operations.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_QuadServices_hh
#define quadrature_QuadServices_hh

#include "ds++/SP.hh"
#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class QuadServices
 * \brief Provide Moment-to-Discrete and Discrete-to-Moment operations.
 *
 * This class provides the Moment-to-Discrete and Discrete-to-Moment matrices
 * as described by Jim Morel in "A Hybrid Collocation-Galerkin-Sn Method for
 * Solving the Boltzmann Transport Equation," Nuclear Science and
 * Engineering, Vol. 101, pp. 72-87 (2001).
 *
 * Some code was derived from the Partisn source file "galerkq.f."  Other
 * pieces were derived from the Snac2 soure file "sncYnm.f90."
 *
 * Multiplying the moment representation of the solution vector by the matrix
 * M() provides a discrete (SN) representation of the solution vector.  The
 * matrix M() is the moment-to-discrete operator, c * Y_(n,m).  Where c *
 * Y_(n,m) is the nth spherical harmonic evaluated at the mth quadrature
 * direction, multiplied by a constant.
 *
 * If L is the order of anisotropic scattering, then we need the (k,l)th
 * spherical harmonic for k=0..L and l=-k..k.
 * 
 * \f$ \mathbf{M}_{m,(k,\ell)} 
 * = \frac{2k+1}{\sum{w_m}}
 * c_{k,\ell} Y_{k,\ell}(\mathbf\Omega_m),\,\,\,\,\,m=1,\ldots,M,
 * \,\,\,\,\,k=0,\ldots, N,\,\,\,\,\,\ell=-k,\ldots,k \f$
 *
 * The constant c_{k,l} has been defined by Morel to be
 *
 * \f$ c_{k,\ell} = \sqrt{ (2-\delta_{\ell 0})
 * \frac{(k-\abs\ell)!}{(k+\abs\ell)!}} \f$
 *
 * \sa QuadServices.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     SP<Quadrature> spQuad = QuadCreator().quadCreate( 
 *                                 QuadCreator::LevelSym, 4 );
 *
 *     \\ Create a QuadServices object that is attached to a specific
 *     \\ quadrature set. 
 *     QuadServices myQS( spQuad );
 *
 * \endcode
 */
/*! 
 * \example quadrature/test/tQuadServices.cc 
 * 
 * description of example
 */
//===========================================================================//

class QuadServices 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! Default constructor assumes that only isotropic scattering is used. 
    explicit QuadServices( rtt_dsxx::SP< const Quadrature > spQuad_);

    //! Copy constructor (the long doxygen description is in the .cc file).
    QuadServices( QuadServices const & rhs );

    //! Destructor.
    virtual ~QuadServices() { /* empty */ };

    // MANIPULATORS
    
    //! Assignment operator for QuadServices.
    QuadServices& operator=( QuadServices const & rhs );

    // ACCESSORS

    std::vector< double > getM() const { return Mmatrix; }

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    //! \bug move to ds++
    template< typename T > T kronecker_delta( T const, T const ) const;
    
    //! \bug move to ds++
    template< typename T > T factorial( T const ) const;
    
    //! Compute the (k,ell) spherical harmonic evaluated at Omega_m.
    double spherical_harmonic( unsigned const ell,
			       int      const k,
			       unsigned const m ) const;    

    //! Compute the (k,ell) Legendre polynomial evaluated at mu_m.
    double legendre_polynomial( unsigned const ell, 
				unsigned const k,
				double   const x ) const;

    //! Build the Mmatrix.
    std::vector< double > computeM() const;
    void computeM_3D( std::vector< double > & Mmat ) const;
    void computeM_2D( std::vector< double > & Mmat ) const;
    void computeM_1D( std::vector< double > & Mmat ) const;
    std::vector< double > computeD() const;

    //! Helper functions to compute coefficients
    inline double gamma( int k, unsigned ell ) const;
    inline double beta(  unsigned ell )        const;
    inline void   cYnm(  std::vector< double > & Mmatrix,
			 unsigned n, int k, unsigned ell ) const;

    // DATA
    rtt_dsxx::SP< const Quadrature > const spQuad;
    unsigned                         const numMoments;
    std::vector< double >            const Mmatrix;
    std::vector< double >            const Dmatrix;

};

} // end namespace rtt_quadrature

#include "QuadServices.i.hh"

#endif // quadrature_QuadServices_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices.hh
//---------------------------------------------------------------------------//
