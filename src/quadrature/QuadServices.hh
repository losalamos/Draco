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

class QuadServices 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    /*!
     * \brief The typedef specifies how the moment index \f$ n \f$ maps
     *        to the double index \f$ (\ell,k) \f$.
     */
    typedef std::pair< unsigned, int > lk_index;

    // CREATORS
    
    //! Default constructor assumes that only isotropic scattering is used. 
    explicit QuadServices( rtt_dsxx::SP< const Quadrature > const spQuad_ );
    
    //! Constructor that allows the user to pick the (k,l) moments to use.
    //! \todo This still needs to be defined.
    QuadServices( rtt_dsxx::SP< const Quadrature > const   spQuad_,
		  std::vector< lk_index >          const & lkMoments_ );

    //! Copy constructor (the long doxygen description is in the .cc file).
    QuadServices( QuadServices const & rhs );

    //! Destructor.
    virtual ~QuadServices() { /* empty */ }

    // MANIPULATORS
    
    //! Assignment operator for QuadServices.
    QuadServices& operator=( QuadServices const & rhs );

    // ACCESSORS

    //! \brief Return the moment-to-discrete operator.
    std::vector< double > getM() const { return Mmatrix; }

    //! \brief Return the discrete-to-moment operator.
    std::vector< double > getD() const { return Dmatrix; }

    std::vector< double > applyM( std::vector< double > const & phi ) const;
    std::vector< double > applyD( std::vector< double > const & psi ) const;

    //! \brief Provide the number of moments used by QuadServices.
    unsigned getNumMoments() const { return numMoments; }

    //! \brief Pretty print vector<T> in a 2D format.
    template< typename T > 
    void print_matrix( std::string           const & matrix_name,
		       std::vector<T>        const & x,
		       std::vector<unsigned> const & dims ) const;
    

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    //! \bug move to ds++
    template< typename T > T kronecker_delta( T const, T const ) const;
    
    //! \bug move to ds++
    template< typename T > T factorial( T const ) const;

    //! Compute the (k,ell) spherical harmonic evaluated at Omega_m.
    double spherical_harmonic( unsigned const m,
			       unsigned const k,
			       int      const ell ) const;    

    //! \brief constuct maps between moment index n and the tuple (k,l).
    // This can optionally be provided by the user.
    std::vector< lk_index > compute_n2lk()    const;
    std::vector< lk_index > compute_n2lk_1D() const;
    std::vector< lk_index > compute_n2lk_2D() const;
    std::vector< lk_index > compute_n2lk_3D() const;

    //! Build the Mmatrix.
    std::vector< double > computeM() const;
    std::vector< double > computeD() const;

    //! Helper functions to compute coefficients
   inline double compute_clk( unsigned k, int ell ) const;
   inline double compute_azimuthalAngle( double const mu,
	  			         double const eta,
				         double const xi ) const;

    //! Checks
    bool D_equals_M_inverse() const;

    // DATA
    rtt_dsxx::SP< const Quadrature > const spQuad;
    unsigned                         const numMoments;
    std::vector< lk_index >          const n2lk;
    std::vector< double >            const Mmatrix;
    std::vector< double >            const Dmatrix;

};

} // end namespace rtt_quadrature

#include "QuadServices.i.hh"

#endif // quadrature_QuadServices_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices.hh
//---------------------------------------------------------------------------//
