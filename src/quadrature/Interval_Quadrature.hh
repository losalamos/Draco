//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   quadrature/Interval_Quadrature.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Legendre quadrature set.
 * \note   Copyright 2016-2019 Triad National Security, LLC. All rights
 *         reserved.  */
//----------------------------------------------------------------------------//

#ifndef __quadrature_Interval_Quadrature_hh__
#define __quadrature_Interval_Quadrature_hh__

#include "Quadrature.hh"

namespace rtt_quadrature {

//============================================================================//
/*!
 * \class Interval_Quadrature
 *
 * \brief A class representing an interval quadrature set.
 *
 * This is an abstraction of all interval (e.g. 1D) angle quadrature sets.
 */
//============================================================================//

class Interval_Quadrature : public Quadrature {
public:
  // CREATORS

  Interval_Quadrature(unsigned const sn_order);

  // ACCESSORS

  virtual Quadrature_Class quadrature_class() const;

  virtual bool has_axis_assignments() const;

  double mu(unsigned const m) const {
    Require(mu_.size() == sn_order());
    Require(m < mu_.size());

    return mu_[m];
  };

  double wt(unsigned const m) const {
    Require(wt_.size() == sn_order());
    Require(m < wt_.size());

    return wt_[m];
  };

  // STATICS

protected:
  using Quadrature::create_ordinates_;

  //! Virtual hook for create_ordinate_set
  virtual vector<Ordinate>
  create_ordinates_(unsigned dimension, Geometry, double norm, unsigned mu_axis,
                    unsigned eta_axis, bool include_starting_directions,
                    bool include_extra_directions) const;

  //! Virtual hook for create_ordinate_set
  virtual vector<Ordinate>
  create_ordinates_(unsigned dimension, Geometry, double norm,
                    bool include_starting_directions,
                    bool include_extra_directions) const;

  //! Virtual hook for create_ordinate_set
  virtual vector<Ordinate> create_level_ordinates_(double norm) const = 0;

  // DATA

  vector<double> mu_, wt_;
};

} // end namespace rtt_quadrature

#endif // __quadrature_Quadrature_hh__

//----------------------------------------------------------------------------//
// end of quadrature/Interval_Quadrature.hh
//----------------------------------------------------------------------------//
