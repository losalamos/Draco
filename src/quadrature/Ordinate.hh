//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Ordinate.hh
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Ordinate_hh
#define quadrature_Ordinate_hh

#include <vector>
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"
#include "mesh_element/Geometry.hh"
#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Ordinate
 * \brief Struct containing angle cosines and weights for an element of
 * an angle quadrature set.
 *
 * Provides a container that represents \f$ \mathbf\Omega_m = \mu_m \mathbf e_x +
 * \eta_m \mathbf e_y + \xi_m \mathbf e_z \f$ plus the associated point weight,
 * $\f w_m \f$.
 */
//===========================================================================//

class Ordinate
{
  public:

    // CREATORS

    //! Create an uninitialized Ordinate.  This is required by the
    //! constructor for vector<Ordinate>.
    Ordinate() : x(0), y(0), z(0), weight(0) {}

    //! Construct an Ordinate from the specified vector and weight.
    inline
    Ordinate(double x, double y, double z, double weight)
        : x(x), y(y), z(z), weight(weight)
    {
        // Ensure(rtt_dsxx::soft_equiv(x*x+y*y+z*z, 1.0));
    }

    //! Construct a 1D Ordinate from the specified angle and weight.
    inline
    Ordinate(double xx, double weight)
    : x(xx), y(0.0), z(0.0), weight(weight)
    {
        Require(xx != 0.0 && xx<=1.0);
    }
    
    // Accessors
    
    double mu()  const { return x; };
    double eta() const { return y; };
    double xi()  const { return z; };
    double wt()  const { return weight; };
    
    //! STL-compatible comparator predicate to sort ordinates by xi
    //! then mu. 
    static
    bool SnCompare(const Ordinate &, const Ordinate &);

    //! Compute a real representation of the spherical harmonics.
    static
    double Y(unsigned l,
             int m,
             Ordinate const &ordinate,
             double const sumwt );

  private:

    // DATA

    // The data must be kept private in order to preserve the norm invariant.
    
    //! Angle cosines for the ordinate.
    double x, y, z;
    //! Quadrature weight for the ordinate.
    double weight;
};

//===========================================================================//
/*!
 * \class OrdinateSet
 * \brief A collection of Ordinates that make up a complete quadrature set.
 *
 * Calculate an ordinate set from a Quadrature.
 *
 */
//===========================================================================//
class OrdinateSet
{
  public:

    // CREATORS

    //! default creator
    OrdinateSet( rtt_dsxx::SP< Quadrature const > const quadrature,
                 rtt_mesh_element::Geometry geometry,
                 unsigned const dimension );

    //! destructor
    virtual ~OrdinateSet(){}

    //! Return the ordinates.
    std::vector<Ordinate> const &getOrdinates() const { return ordinates; }

    //! Return the quadrature on which this OrdinateSet is based.
    rtt_dsxx::SP< const Quadrature > getQuadrature( void ) const
    {
        return quadrature;
    }

    //! Return the geometry.
    rtt_mesh_element::Geometry getGeometry() const { return geometry; }

    //! Return the dimension.
    unsigned getDimension() const { return dimension; }
    
  private:

    // Helper functions called by the constructor.
    void create_set_from_1d_quadrature( void );
    void create_set_from_2d_quadrature_for_2d_mesh( void );
    void create_set_from_2d_quadrature_for_1d_mesh( void );

    // DATA
    std::vector<Ordinate> ordinates;
    rtt_dsxx::SP< Quadrature const > const quadrature;
    rtt_mesh_element::Geometry const geometry;
    unsigned const dimension;
    
};

} // end namespace rtt_quadrature

#endif // quadrature_Ordinate_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Ordinate.hh
//---------------------------------------------------------------------------//
