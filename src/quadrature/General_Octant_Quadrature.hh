//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/General_Octant_Quadrature.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#ifndef quadrature_General_Octant_Quadrature_hh
#define quadrature_General_Octant_Quadrature_hh

#include "Octant_Quadrature.hh"

namespace rtt_quadrature {

//===========================================================================//
/*!
 * \class General_Octant_Quadrature
 * \brief A class to encapsulate a client-defined ordinate set.
 */
//===========================================================================//

class General_Octant_Quadrature : public Octant_Quadrature {
public:
  // CREATORS
  DLL_PUBLIC_quadrature
  General_Octant_Quadrature(unsigned const sn_order, vector<double> const &mu,
                            vector<double> const &eta, vector<double> const &xi,
                            vector<double> const &wt, unsigned number_of_levels,
                            Quadrature_Class);

  General_Octant_Quadrature(); // disable default construction

  // ACCESSORS

  vector<double> const &mu() const { return mu_; }
  vector<double> const &eta() const { return eta_; }
  vector<double> const &xi() const { return xi_; }
  vector<double> const &wt() const { return wt_; }

  // SERVICES

  // These functions override the virtual member functions specifed in the
  // parent class Quadrature.

  string name() const;

  string parse_name() const;

  Quadrature_Class quadrature_class() const;

  unsigned number_of_levels() const;

  string as_text(string const &indent) const;

  bool is_open_interval() const;

  bool check_class_invariants() const;

  // STATICS

  static std::shared_ptr<Quadrature> parse(Token_Stream &tokens);

private:
  // IMPLEMENTATION

  //! Virtual hook for create_ordinate_set
  virtual void create_octant_ordinates_(vector<double> &mu, vector<double> &eta,
                                        vector<double> &wt) const;

  // DATA
  vector<double> mu_, eta_, xi_, wt_;
  unsigned number_of_levels_;
  Quadrature_Class quadrature_class_;
  bool is_open_interval_;
};

} // end namespace rtt_quadrature

#endif // quadrature_General_Octant_Quadrature_hh

//---------------------------------------------------------------------------//
// end of quadrature/General_Octant_Quadrature.hh
//---------------------------------------------------------------------------//
