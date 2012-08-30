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

#include <iostream>
#include <iomanip>

#include "ds++/SP.hh"
#include "Ordinate.hh"
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
    QuadServices( std::vector<Ordinate> const &ordinates,
                  double                const norm,
                  unsigned              const dimension,
                  unsigned              const expansionOrder = 0);

    //! Copy constructor (the long doxygen description is in the .cc file).
    QuadServices( QuadServices const & rhs );

    //! Destructor.
    virtual ~QuadServices() { /* empty */ }

    // MANIPULATORS

    std::vector< double > applyM( std::vector< double > const & phi) const;
    std::vector< double > applyD( std::vector< double > const & psi) const;
    
    //! Assignment operator for QuadServices.
    QuadServices& operator=( QuadServices const & rhs );

    //! Compute extra "moment-to-discrete" entries that are needed for starting direction ordinates.
    double augmentM( unsigned moment, Ordinate const & Omega, std::vector< lk_index > const & n2lk ) const;

    // ACCESSORS

    //! \brief Provide the number of moments used by QuadServices.
    unsigned getExpansionOrder() const { return expansionOrder_; }

    //! \brief Provide the quadrature normalization used by QuadServices.
    double getNorm() const { return norm_; }

    //! \brief Provide the dimensionality used by QuadServices.
    unsigned getDimension() const { return dimension_; }

    //! \brief Provide the dimensionality used by QuadServices.
    std::vector<Ordinate> const &getOrdinates() const { return ordinates_; }

    //! \brief Provide the number of moments used by QuadServices.
    unsigned getNumOrdinates() const { return ordinates_.size(); }

    //! \brief Pretty print vector<T> in a 2D format.
    void print_matrix( std::string           const & matrix_name,
		       std::vector<double>   const & x,
		       std::vector<unsigned> const & dims ) const
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
            
            std::cout << "{ ";
            
            for( unsigned j=0; j<dims[0]-1; ++j )
                std::cout << std::setprecision(10) << x[j+dims[0]*i] << ", ";
            
            std::cout << std::setprecision(10) << x[dims[0]-1+dims[0]*i] << " }." << std::endl;
        }
        std::cout << std::endl;
        return;
    }
    
    //! \brief Provide the maximum available expansion order.
    virtual unsigned getMaxExpansionOrder() const = 0;

    //! \brief Provide the number of moments used by QuadServices.
    virtual unsigned getNumMoments() const = 0;

    //! \brief Provide the number of moments in each expansion order
    virtual std::vector<unsigned> getMoments() const = 0;

    //! \brief Provide the number of moments used by QuadServices.
    virtual std::vector< lk_index > const &get_n2lk() const = 0;

    //! \brief Return the (l,k) index pair associated with moment index n.
    virtual lk_index lkPair( unsigned n ) const = 0;

    //! \brief Return the moment-to-discrete operator.
    std::vector< double > getM() const;

    //! \brief Return the discrete-to-moment operator.
    std::vector< double > getD() const;

    unsigned max_available_expansion_order(std::vector< lk_index > const &n2lk);

    std::vector<unsigned> compute_moments(unsigned const L,
                                          std::vector< lk_index > const &n2lk);

    double compute_azimuthalAngle( double const mu,
                                   double const eta,
                                   double const xi ) const;
    
    // STATICS

    //! Helper function to check validity of LU matrix
    static bool diagonal_not_zero( std::vector<double> const & vec,
                                   int m, int n ) ;
    
    static void moment_to_flux(double Phi_10,
                               double &Fz)
    {
        Fz = Phi_10;
    }
    
    static void moment_to_flux(double Phi_1p1,
                               double Phi_10,
                               double &Fx,
                               double &Fz)
    {
        Fx = Phi_1p1;
        Fz = -Phi_10;
    }
    
    static void moment_to_flux(double Phi_1m1,
                               double Phi_10,
                               double Phi_1p1,
                               double &Fx,
                               double &Fy,
                               double &Fz)
    {
        Fx = -Phi_1p1;
        Fy = -Phi_1m1;
        Fz = Phi_10;
    }

    //! \brief Return the moment-to-discrete operator.
    virtual std::vector< double > getM_() const = 0;

    //! \brief Return the discrete-to-moment operator.
    virtual std::vector< double > getD_() const = 0;

  protected:

    std::vector< lk_index > compute_n2lk_1D( unsigned L);

  private:
    
    
    // legacy calculation, not currently implemented
    std::vector< double > computeD_SVD(std::vector<Ordinate> const &ordinates,
                                       std::vector< lk_index > const &n2lk,
                                       std::vector<double> const &M,
                                       unsigned const dim,
                                       double const sumwt);

    
    // DATA
    std::vector<Ordinate>     const ordinates_;
    double                    const norm_;
    unsigned                  const dimension_;
    unsigned                  const expansionOrder_;
};

} // end namespace rtt_quadrature

#endif // quadrature_QuadServices_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices.hh
//---------------------------------------------------------------------------//
