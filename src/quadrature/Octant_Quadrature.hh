//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Octant_Quadrature.hh
 * \author Kent Budge
 * \date   Friday, Nov 30, 2012, 08:28 am
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef quadrature_Octant_Quadrature_hh
#define quadrature_Octant_Quadrature_hh

#include "Quadrature.hh"

namespace rtt_quadrature {

//===========================================================================//
/*!
 * \class Octant_Quadrature
 * \brief This is an abstract class representing all quadratures over the unit
 * sphere.
 *
 * At present, all our unit sphere quadratures are symmetric in octants, though
 * we will likely relax this restriction in the future.
 *
 * For level quadratures, the levels must be in the xi direction cosine. The
 * user may override the default axis assignments when he constructs an
 * Ordinate_Set or an Ordinate_Space from the Octant_Quadrature.
 */
//===========================================================================//

class Octant_Quadrature : public Quadrature {
public:
  // CREATORS

  Octant_Quadrature(unsigned const sn_order)
      : Quadrature(sn_order), has_axis_assignments_(false), mu_axis_(),
        eta_axis_() { /* empty */
  }

  Octant_Quadrature(unsigned const sn_order, unsigned const mu_axis,
                    unsigned const eta_axis)
      : Quadrature(sn_order), has_axis_assignments_(true), mu_axis_(mu_axis),
        eta_axis_(eta_axis) { /* empty */
  }

  // ACCESSORS

  // SERVICES
  DLL_PUBLIC_quadrature virtual bool has_axis_assignments() const;

protected:
  virtual string as_text(string const &indent) const = 0;

  // IMPLEMENTATION

  //! Virtual hook for create_ordinate_set
  DLL_PUBLIC_quadrature virtual void
  create_octant_ordinates_(vector<double> &mu, vector<double> &eta,
                           vector<double> &wt) const = 0;

  // STATICS

  static void parse(Token_Stream &tokens, bool &has_axis_assignments,
                    unsigned &mu_axis, unsigned &eta_axis);

private:
  // IMPLEMENTATION

  using Quadrature::create_ordinates_;

  //! Virtual hook for create_ordinates
  DLL_PUBLIC_quadrature virtual vector<Ordinate>
  create_ordinates_(unsigned dimension, Geometry, double norm,
                    bool include_starting_directions,
                    bool include_extra_directions) const;

  //! Virtual hook for create_ordinate_set
  DLL_PUBLIC_quadrature virtual vector<Ordinate>
  create_ordinates_(unsigned dimension, Geometry, double norm, unsigned mu_axis,
                    unsigned eta_axis, bool include_starting_directions,
                    bool include_extra_directions) const;

  // DATA

  bool has_axis_assignments_;
  unsigned mu_axis_, eta_axis_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Octant_Quadrature_hh

//---------------------------------------------------------------------------//
// end of quadrature/Octant_Quadrature.hh
//---------------------------------------------------------------------------//
