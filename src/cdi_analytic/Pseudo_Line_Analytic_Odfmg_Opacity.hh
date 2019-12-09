//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Pseudo_Line_Analytic_Odfmg_Opacity class definition.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Pseudo_Line_Analytic_Odfmg_Opacity_hh__
#define __cdi_analytic_Pseudo_Line_Analytic_Odfmg_Opacity_hh__

#include "Analytic_Odfmg_Opacity.hh"
#include "Pseudo_Line_Base.hh"

namespace rtt_cdi_analytic {
using std::pair;

//===========================================================================//
/*!
 * \class Pseudo_Line_Analytic_Odfmg_Opacity
 *
 * \brief Derived rtt_cdi::OdfmgOpacity class for analytic opacities.
 *
 * Primarily code from Analytic_Multigroup_Opacity.
 */
//===========================================================================//

class Pseudo_Line_Analytic_Odfmg_Opacity : public Analytic_Odfmg_Opacity,
                                           public Pseudo_Line_Base {
private:
  Averaging averaging_;
  unsigned qpoints_;
  vector<pair<double, pair<double, double>>> baseline_;

  void precalculate(vector<double> const &groups, vector<double> const &bands,
                    double Tref);

public:
  // Constructor.
  Pseudo_Line_Analytic_Odfmg_Opacity(
      const sf_double &groups, const sf_double &bands,
      rtt_cdi::Reaction reaction_in,
      std::shared_ptr<Expression const> const &cont, int number_of_lines,
      double line_peak, double line_width, int number_of_edges,
      double edge_ratio, double Tref, double Tpow, double emin, double emax,
      Averaging averaging, unsigned qpoints, unsigned seed);

  /*!
   * \brief Constructor 2
   * \bug No doxygen documentation
   * \bug No unit test in Draco but used by Capsaicin in
   *      thermal_data/microphysics_parser.cc.
   */
  Pseudo_Line_Analytic_Odfmg_Opacity(
      const sf_double &groups, const sf_double &bands,
      rtt_cdi::Reaction reaction_in, string const &cont_file,
      int number_of_lines, double line_peak, double line_width,
      int number_of_edges, double edge_ratio, double Tref, double Tpow,
      double emin, double emax, Averaging averaging, unsigned qpoints,
      unsigned seed);

  /*!
   * \brief Constructor 3
   * \bug No doxygen documentation
   * \bug No unit test in Draco but used by Capsaicin in
   *      thermal_data/microphysics_parser.cc.
   */
  Pseudo_Line_Analytic_Odfmg_Opacity(
      const sf_double &groups, const sf_double &bands,
      rtt_cdi::Reaction reaction_in, double nu0, double C, double Bn, double Bd,
      double R, int number_of_lines, double line_peak, double line_width,
      int number_of_edges, double edge_ratio, double Tref, double Tpow,
      double emin, double emax, Averaging averaging, unsigned qpoints,
      unsigned seed);

  // Constructor for packed Pseudo_Line_Analytic_Odfmg_Opacities
  explicit Pseudo_Line_Analytic_Odfmg_Opacity(const sf_char &);

  std::vector<std::vector<double>> getOpacity(double targetTemperature,
                                              double /*targetDensity*/) const;

  std::vector<std::vector<std::vector<double>>>
  getOpacity(const std::vector<double> &targetTemperature,
             double targetDensity) const;

  std::vector<std::vector<std::vector<double>>>
  getOpacity(double targetTemperature,
             const std::vector<double> &targetDensity) const;

  // Get the data description of the opacity.
  std_string getDataDescriptor() const;

  // Pack the Pseudo_Line_Analytic_Odfmg_Opacity into a character string.
  sf_char pack() const;
};

} // end namespace rtt_cdi_analytic

#endif // __cdi_analytic_Pseudo_Line_Analytic_Odfmg_Opacity_hh__

//---------------------------------------------------------------------------//
// end of cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.hh
//---------------------------------------------------------------------------//
