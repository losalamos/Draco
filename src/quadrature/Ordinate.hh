//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Ordinate.hh
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   Copyright ©  2006-2010 Los Alamos National Security, LLC. All rights
 *         reserved. 
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
 * \brief Class containing angle cosines and weights for an element of
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
    Ordinate() : x_(0), y_(0), z_(0), weight_(0) {}

    //! Construct an Ordinate from the specified vector and weight.
    Ordinate(double const x,
             double const y,
             double const z,
             double const weight)
        :
        x_(x), y_(y), z_(z), weight_(weight)
    {
    }

    //! Construct a 1D Ordinate from the specified angle and weight.
    inline
    Ordinate(double const xx, double const weight)
    : x_(xx), y_(0.0), z_(0.0), weight_(weight)
    {
        Require(xx != 0.0 && xx<=1.0);
    }
    
    // Accessors
    
    double mu()  const { return x_; };
    double eta() const { return y_; };
    double xi()  const { return z_; };
    double wt()  const { return weight_; };
    
    //! STL-compatible comparator predicate to sort ordinates by xi
    //! then mu. 
    static
    bool SnCompare(const Ordinate &, const Ordinate &);
    
    //! STL-compatible comparator predicate to sort ordinates into PARTISN 2-D
    //! axisymmetric ordering.
    static
    bool SnComparePARTISN2a(const Ordinate &, const Ordinate &);
    
    //! STL-compatible comparator predicate to sort ordinates into PARTISN 3-D
    //! ordering.
    static
    bool SnComparePARTISN3(const Ordinate &, const Ordinate &);

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
    double x_, y_, z_;
    //! Quadrature weight for the ordinate.
    double weight_;
};

//---------------------------------------------------------------------------//
//! Test ordinates for equality
inline bool operator==(Ordinate const &a, Ordinate const &b)
{
    return
        a.mu()==b.mu() &&
        a.eta()==b.eta() &&
        a.xi()==b.xi() &&
        a.wt()==b.wt();
}

//! Typedef for ordinate comparator functions
typedef bool (*comparator_t)(Ordinate const &, Ordinate const &);

//===========================================================================//
/*!
 * \class OrdinateSet
 * \brief A collection of Ordinates that make up a complete quadrature set.
 */
//===========================================================================//
class OrdinateSet
{
  public:

    // CREATORS

    //! default creator
    OrdinateSet( rtt_dsxx::SP< Quadrature const > const quadrature,
                 rtt_mesh_element::Geometry geometry,
                 unsigned const dimension,
                 bool const extra_starting_directions=false);


    //! destructor
    virtual ~OrdinateSet(){}

    // ACCESSORS

    //! Return the ordinates.
    std::vector<Ordinate> const &getOrdinates() const { return ordinates_; }

    //! Return the quadrature on which this OrdinateSet is based.
    rtt_dsxx::SP< const Quadrature > getQuadrature() const
    {
        Ensure(quadrature_ != rtt_dsxx::SP<Quadrature>());
        
        return quadrature_;
    }

    //! Return the geometry.
    rtt_mesh_element::Geometry getGeometry() const { return geometry_; }

    //! Return the dimension.
    unsigned getDimension() const { return dimension_; }

    //! Return the norm.
    double getNorm() const { return norm_; }

    //! Return the ordering operator.
    comparator_t getComparator() const { return comparator_; }

    //! Return the ordering operator.
    bool extra_starting_directions() const { return extra_starting_directions_; }

    bool check_class_invariants() const;
    
  private:

    // Helper functions called by the constructor.
    void create_set_from_1d_quadrature();
    void create_set_from_2d_quadrature_for_2d_mesh();
    void create_set_from_2d_quadrature_for_1d_mesh();
    void create_set_from_2d_quadrature_for_3d_mesh();
    void create_set_from_3d_quadrature_for_3d_mesh();

    // DATA
    std::vector<Ordinate> ordinates_;
    rtt_dsxx::SP< Quadrature const > quadrature_;
    rtt_mesh_element::Geometry geometry_;
    unsigned dimension_;
    double norm_;
    comparator_t comparator_;
    bool const extra_starting_directions_;
    
};

} // end namespace rtt_quadrature

#endif // quadrature_Ordinate_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Ordinate.hh
//---------------------------------------------------------------------------//
