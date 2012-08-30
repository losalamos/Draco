//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Ordinate.hh
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   Copyright (C)  2006-2012 Los Alamos National Security, LLC.
 *         All rights reserved. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Ordinate_hh
#define quadrature_Ordinate_hh

#include <vector>

#include "mesh_element/Geometry.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

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

    // STATIC
    
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

} // end namespace rtt_quadrature

#endif // quadrature_Ordinate_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Ordinate.hh
//---------------------------------------------------------------------------//
