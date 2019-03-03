//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Sn_Ordinate_Space.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of class Sn_Ordinate_Space
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC
 */
//---------------------------------------------------------------------------------------//
// $Id: Sn_Ordinate_Space.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef quadrature_Sn_Ordinate_Space_hh
#define quadrature_Sn_Ordinate_Space_hh

#include "Ordinate_Space.hh"

namespace rtt_quadrature {
using std::ostream;

//=======================================================================================//
/*!
 * \class Sn_Ordinate_Space
 * \brief Represents ordinate operators for a conventional Sn moment space.
 *
 * The moment space contains all moments up to the specified scattering order,
 * and the moment to discrete and discrete to moment operators are calculated
 * in a straightforward manner from the Ylm and the weight associated with
 * each ordinate direction.
 */
//=======================================================================================//

class Sn_Ordinate_Space : public Ordinate_Space {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Specify the ordinate quadrature with defaults.
  Sn_Ordinate_Space(unsigned dimension, Geometry geometry,
                    vector<Ordinate> const &, int expansion_order,
                    bool extra_starting_directions = false,
                    Ordering ordering = LEVEL_ORDERED);

  // MANIPULATORS

  // ACCESSORS

  bool check_class_invariants() const;

  // SERVICES

  virtual QIM quadrature_interpolation_model() const;

  //! Return the discrete to moment transform matrix
  virtual vector<double> D() const;

  //! Return the moment to discrete transform matrix
  virtual vector<double> M() const;

  // STATICS

protected:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  virtual vector<Moment> compute_n2lk_1D_(Quadrature_Class, unsigned sn_order);

  virtual vector<Moment> compute_n2lk_1Da_(Quadrature_Class, unsigned sn_order);

  virtual vector<Moment> compute_n2lk_2D_(Quadrature_Class, unsigned sn_order);

  virtual vector<Moment> compute_n2lk_2Da_(Quadrature_Class, unsigned sn_order);

  virtual vector<Moment> compute_n2lk_3D_(Quadrature_Class, unsigned sn_order);

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  void compute_M();
  void compute_D();

  // DATA

  vector<double> D_;
  //! Moment to discrete matrix
  vector<double> M_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Sn_Ordinate_Space_hh

//---------------------------------------------------------------------------------------//
//              end of quadrature/Sn_Ordinate_Space.hh
//---------------------------------------------------------------------------------------//
