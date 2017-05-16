//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Pseudo_Line_Analytic_Odfmg_Opacity class member definitions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Pseudo_Line_Analytic_Odfmg_Opacity.hh"
#include "Pseudo_Line_Analytic_MultigroupOpacity.hh"
#include "cdi/CDI.hh"
#include "ds++/DracoMath.hh"
#include "ds++/Packing_Utils.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"
#include "units/PhysicalConstantsSI.hh"
#include <fstream>

namespace rtt_cdi_analytic {
using namespace std;
using namespace rtt_dsxx;
using namespace rtt_ode;
using namespace rtt_cdi;

rtt_parser::Unit const keV = {2, 1, -2, 0, 0,
                              0, 0, 0,  0, 1e3 * rtt_units::electronChargeSI};

//---------------------------------------------------------------------------//
void Pseudo_Line_Analytic_Odfmg_Opacity::precalculate(
    vector<double> const &groups, vector<double> const &bands,
    double const Tref) {
  // Precalculate basic opacities

  unsigned const number_of_groups = groups.size() - 1U;
  unsigned const N = qpoints_;
  baseline_.resize(N * number_of_groups);

#if 1
  ofstream out("pseudo.dat");
#endif

  double g1 = groups[0];
  for (unsigned g = 0; g < number_of_groups; ++g) {
    double const g0 = g1;
    g1 = groups[g + 1];
    double const delt = (g1 - g0) / N;
    double x = g0 - 0.5 * delt;
    for (unsigned iq = 0; iq < N; ++iq) {
      x += delt;
      baseline_[iq + N * g].first = monoOpacity(x, Tref);
      baseline_[iq + N * g].second.first = x - 0.5 * delt;
      baseline_[iq + N * g].second.second = x + 0.5 * delt;

#if 1
      out << x << ' ' << baseline_[iq + N * g].first << endl;
#endif
    }
    if (bands.size() > 2) {
      sort(baseline_.begin() + N * g, baseline_.begin() + N * (g + 1));
    }
  }
}

//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_Odfmg_Opacity::Pseudo_Line_Analytic_Odfmg_Opacity(
    const sf_double &groups, const sf_double &bands,
    rtt_cdi::Reaction reaction_in,
    std::shared_ptr<Expression const> const &continuum, int number_of_lines,
    double line_peak, double line_width, int number_of_edges, double edge_ratio,
    double Tref, double Tpow, double emin, double emax, Averaging averaging,
    unsigned qpoints, unsigned seed)
    : Analytic_Odfmg_Opacity(groups, bands, reaction_in),
      Pseudo_Line_Base(continuum, number_of_lines, line_peak, line_width,
                       number_of_edges, edge_ratio, Tref, Tpow, emin, emax,
                       seed),
      averaging_(averaging), qpoints_(qpoints), baseline_() {
  Require(qpoints > 0);

  precalculate(groups, bands, Tref);
}

//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_Odfmg_Opacity::Pseudo_Line_Analytic_Odfmg_Opacity(
    const sf_double &groups, const sf_double &bands,
    rtt_cdi::Reaction reaction_in, string const &cont_file, int number_of_lines,
    double line_peak, double line_width, int number_of_edges, double edge_ratio,
    double Tref, double Tpow, double emin, double emax, Averaging averaging,
    unsigned qpoints, unsigned seed)
    : Analytic_Odfmg_Opacity(groups, bands, reaction_in),
      Pseudo_Line_Base(cont_file, number_of_lines, line_peak, line_width,
                       number_of_edges, edge_ratio, Tref, Tpow, emin, emax,
                       seed),
      averaging_(averaging), qpoints_(qpoints), baseline_() {
  Require(qpoints > 0);

  precalculate(groups, bands, Tref);
}

//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_Odfmg_Opacity::Pseudo_Line_Analytic_Odfmg_Opacity(
    const sf_double &groups, const sf_double &bands,
    rtt_cdi::Reaction reaction_in, double nu0, double C, double Bn, double Bd,
    double R, int number_of_lines, double line_peak, double line_width,
    int number_of_edges, double edge_ratio, double Tref, double Tpow,
    double emin, double emax, Averaging averaging, unsigned qpoints,
    unsigned seed)
    : Analytic_Odfmg_Opacity(groups, bands, reaction_in),
      Pseudo_Line_Base(nu0, C, Bn, Bd, R, number_of_lines, line_peak,
                       line_width, number_of_edges, edge_ratio, Tref, Tpow,
                       emin, emax, seed),
      averaging_(averaging), qpoints_(qpoints), baseline_() {
  Require(qpoints > 0);

  precalculate(groups, bands, Tref);
}

//---------------------------------------------------------------------------//
// OPACITY INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the group opacities given a scalar temperature and density.
 *
 * Given a scalar temperature and density, return the group opacities
 * (vector<double>) for the reaction type specified by the constructor.  The
 * analytic opacity model is specified in the constructor
 * (Pseudo_Line_Analytic_Odfmg_Opacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return group opacities (coefficients) in cm^2/g
 *
 */
std::vector<std::vector<double>>
Pseudo_Line_Analytic_Odfmg_Opacity::getOpacity(double T,
                                               double /* rho */) const {
  sf_double const &group_bounds = this->getGroupBoundaries();
  sf_double const &bands = this->getBandBoundaries();
  unsigned const number_of_groups = group_bounds.size() - 1U;
  unsigned const bands_per_group = bands.size() - 1U;
  vector<vector<double>> Result(number_of_groups,
                                vector<double>(bands_per_group));

  unsigned const N = qpoints_;

  double const Tf = pow(T / Tref(), Tpow());

  switch (averaging_) {
  case NONE:
    if (N > 1) {
      for (unsigned g = 0; g < number_of_groups; ++g) {
        double b1 = bands[0];
        for (unsigned b = 0; b < bands_per_group; ++b) {
          double const b0 = b1;
          b1 = bands[b + 1];
          // little bit of an elaborate interpolation dance here
          double f = 0.5 * (b0 + b1) * N - 0.5;
          unsigned const i = static_cast<unsigned int>(f);
          f -= i;
          Check(i < N);
          Check(f >= -0.5 && f <= 0.5);
          Result[g][b] = Tf * ((1.0 - f) * baseline_[i + N * g].first +
                               f * baseline_[i + 1 + N * g].first);
        }
      }
    } else {
      Check(bands_per_group == 1);
      for (unsigned g = 0; g < number_of_groups; ++g) {
        Result[g][0] = Tf * baseline_[g].first;
      }
    }
    break;

  case ROSSELAND:
    for (unsigned g = 0; g < number_of_groups; ++g) {
      double b1 = bands[0];
      for (unsigned b = 0; b < bands_per_group; ++b) {
        double const b0 = b1;
        b1 = bands[b + 1];
        double t = 0.0, w = 0.0;

        unsigned const q0 = static_cast<unsigned>(b0 * N);
        unsigned const q1 = static_cast<unsigned>(b1 * N);
        for (unsigned q = q0; q < q1; ++q) {
          double const x0 = baseline_[q + N * g].second.first;
          double const x1 = baseline_[q + N * g].second.second;

          double weight =
              CDI::integrateRosselandSpectrum(x0 / keV.conv, x1 / keV.conv, T) +
              numeric_limits<double>::min();

          t += weight / baseline_[q + N * g].first;
          w += weight;
        }

        Result[g][b] = w * Tf / t;
      }
    }
    break;

  case PLANCK:
    for (unsigned g = 0; g < number_of_groups; ++g) {
      double b1 = bands[0];
      for (unsigned b = 0; b < bands_per_group; ++b) {
        double const b0 = b1;
        b1 = bands[b + 1];
        double t = 0.0, w = 0.0;

        unsigned const q0 = static_cast<unsigned>(b0 * N);
        unsigned const q1 = static_cast<unsigned>(b1 * N);
        for (unsigned q = q0; q < q1; ++q) {
          double const x0 = baseline_[q + N * g].second.first;
          double const x1 = baseline_[q + N * g].second.second;

          double weight =
              CDI::integratePlanckSpectrum(x0 / keV.conv, x1 / keV.conv, T) +
              numeric_limits<double>::min();

          t += weight * baseline_[q + N * g].first;
          w += weight;
        }

        Result[g][b] = t * Tf / w;
      }
    }
    break;

  default:
    Insist(false, "bad case");
  }

  return Result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband
 *     opacity 2-D vectors that correspond to the provided vector of
 *     temperatures and a single density value.
 */
std::vector<std::vector<std::vector<double>>>
Pseudo_Line_Analytic_Odfmg_Opacity::getOpacity(
    const std::vector<double> &targetTemperature, double targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(
      targetTemperature.size());

  for (size_t i = 0; i < targetTemperature.size(); ++i) {
    opacity[i] = getOpacity(targetTemperature[i], targetDensity);
  }
  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband
 *     opacity 2-D vectors that correspond to the provided
 *     temperature and a vector of density values.
 */
std::vector<std::vector<std::vector<double>>>
Pseudo_Line_Analytic_Odfmg_Opacity::getOpacity(
    double targetTemperature, const std::vector<double> &targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(targetDensity.size());

  //call our regular getOpacity function for every target density
  for (size_t i = 0; i < targetDensity.size(); ++i) {
    opacity[i] = getOpacity(targetTemperature, targetDensity[i]);
  }
  return opacity;
}

//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_Odfmg_Opacity::std_string
Pseudo_Line_Analytic_Odfmg_Opacity::getDataDescriptor() const {
  std_string descriptor;

  rtt_cdi::Reaction const reaction = getReactionType();

  if (reaction == rtt_cdi::TOTAL)
    descriptor = "Pseudo Line Odfmg Total";
  else if (reaction == rtt_cdi::ABSORPTION)
    descriptor = "Pseudo Line Odfmg Absorption";
  else if (reaction == rtt_cdi::SCATTERING)
    descriptor = "Pseudo Line Odfmg Scattering";
  else {
    Insist(0, "Invalid Pseudo Line Odfmg model opacity!");
  }

  return descriptor;
}

//---------------------------------------------------------------------------//
// Packing function

Analytic_MultigroupOpacity::sf_char
Pseudo_Line_Analytic_Odfmg_Opacity::pack() const {
  sf_char const pdata = Analytic_Odfmg_Opacity::pack();
  sf_char const pdata2 = Pseudo_Line_Base::pack();

  sf_char Result(pdata.size() + pdata2.size());
  copy(pdata.begin(), pdata.end(), Result.begin());
  copy(pdata2.begin(), pdata2.end(), Result.begin() + pdata.size());
  return pdata;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
// end of Pseudo_Line_Analytic_Odfmg_Opacity.cc
//---------------------------------------------------------------------------//
