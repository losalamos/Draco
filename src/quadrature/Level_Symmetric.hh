//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Level_Symmetric.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef quadrature_Level_Symmetric_hh
#define quadrature_Level_Symmetric_hh

#include "Octant_Quadrature.hh"

namespace rtt_quadrature {

//===========================================================================//
/*!
 * \class Level_Symmetric
 * \brief A class to encapsulate a 3D Level Symmetric quadrature set.
 */
//===========================================================================//

class Level_Symmetric : public Octant_Quadrature {
public:
  // CREATORS

  // The default values for snOrder_ and norm_ were set in QuadCreator.
  explicit Level_Symmetric(unsigned sn_order)
      : Octant_Quadrature(sn_order)

  {
    Require(sn_order > 0 && sn_order % 2 == 0);
  }

  Level_Symmetric(); // disable default construction

  // ACCESSORS

  // SERVICES

  // These functions override the virtual member functions specifed in the
  // parent class Quadrature.

  DLL_PUBLIC_quadrature string name() const;

  DLL_PUBLIC_quadrature string parse_name() const;

  DLL_PUBLIC_quadrature Quadrature_Class quadrature_class() const;

  DLL_PUBLIC_quadrature unsigned number_of_levels() const;

  DLL_PUBLIC_quadrature string as_text(string const &indent) const;

  // STATICS

  static std::shared_ptr<Quadrature> parse(Token_Stream &tokens);

private:
  // IMPLEMENTATION

  //! Virtual hook for create_ordinate_set
  DLL_PUBLIC_quadrature virtual void
  create_octant_ordinates_(vector<double> &mu, vector<double> &eta,
                           vector<double> &wt) const;

  // DATA
};

} // end namespace rtt_quadrature

#endif // quadrature_Level_Symmetric_hh

//---------------------------------------------------------------------------//
// end of quadrature/Level_Symmetric.hh
//---------------------------------------------------------------------------//
