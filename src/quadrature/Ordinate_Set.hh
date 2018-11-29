//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Ordinate_Set.hh
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   Copyright (C)  2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Ordinate.hh 6607 2012-06-14 22:31:45Z kellyt $
//---------------------------------------------------------------------------//

#ifndef quadrature_OrdinateSet_hh
#define quadrature_OrdinateSet_hh

#include <vector>

#include "Ordinate.hh"
#include "mesh_element/Geometry.hh"

namespace rtt_quadrature {
using rtt_mesh_element::Geometry;
using std::vector;

//===========================================================================//
/*!
 * \class Ordinate_Set
 *
 * \brief An ordered collection of Ordinates that make up a complete
 * quadrature set for a specified geometry.
 */
//===========================================================================//

class DLL_PUBLIC_quadrature Ordinate_Set {
public:
  // ENUMERATIONS

  //! Ordering of ordinates.
  enum Ordering {
    LEVEL_ORDERED, // Ordered xi, then mu, then eta

    OCTANT_ORDERED // Ordered by xi halfsphere, then eta halfsphere, then
                   // mu halfsphere, then absolute value of xi, then
                   // absolute value of eta, then absolute value of mu
  };

  // CREATORS

  Ordinate_Set(unsigned const dimension, Geometry const geometry,
               vector<Ordinate> const &ordinates,
               bool const has_starting_directions = true,
               bool const has_extra_starting_directions = false,
               Ordering const ordering = LEVEL_ORDERED);

  //! destructor
  virtual ~Ordinate_Set() {}

  // ACCESSORS

  //! Return the ordinates.
  vector<Ordinate> const &ordinates() const { return ordinates_; }

  //! Return the ordinates.
  vector<Ordinate> &ordinates() { return ordinates_; }

  //! Return the geometry.
  Geometry geometry() const { return geometry_; }

  //! Return the dimension.
  unsigned dimension() const { return dimension_; }

  //! Return the norm.
  double norm() const { return norm_; }

  //! Does this set include starting directions?
  bool has_starting_directions() const { return has_starting_directions_; }

  //! Does this set include extra starting directions?
  bool has_extra_starting_directions() const {
    return has_extra_starting_directions_;
  }

  Ordering ordering() const { return ordering_; }

  // SERVICES

  //! Print a description of the quadrature set to the terminal.
  void display() const;

  bool check_class_invariants() const;

  // STATIC

  //! Comparator for enforcing LEVEL_ORDERING
  static bool level_compare(Ordinate const &a, Ordinate const &b);

  //! Comparator for enforcing OCTANT_ORDERING
  static bool octant_compare(Ordinate const &a, Ordinate const &b);

private:
  // DATA

  // initialized
  Geometry geometry_;
  unsigned dimension_;
  bool has_starting_directions_;
  bool has_extra_starting_directions_;
  Ordering ordering_;
  double norm_;
  vector<Ordinate> ordinates_;
};

} // end namespace rtt_quadrature

#endif // quadrature_OrdinateSet_hh

//---------------------------------------------------------------------------//
// end of quadrature/OrdinateSet.hh
//---------------------------------------------------------------------------//
