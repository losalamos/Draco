//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/F_eta.cc
 * \author Kent Budge
 * \date   Mon Sep 20 15:01:53 2004
 * \brief  Implementation of F_eta.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * These routines are based on C routines from Numerical Recipes.
 */
//---------------------------------------------------------------------------//

#include "F_eta.hh"
#include "Factorial.hh"
#include "ds++/DracoMath.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"
#include "units/PhysicalConstants.hh"

namespace rtt_sf {
using namespace std;
using namespace rtt_dsxx;
using namespace rtt_ode;
using rtt_units::PI;

// Parametrization of integrand and inverse functions
static double leta, lgamma;

//---------------------------------------------------------------------------//
/*!
 * \brief Integrand of integral representation of F_eta
 *
 * \param x[in] Argument
 *
 * \return Value of integrand at the argument,
 *         \f$\frac{(x^2+2x)^{3/2}}{e^\frac{x-\eta}{\gamma}+1}\f$
 *
 * \post \c Result>=0
 */
static double Feta_integrand(double x) {
  double const y = x * x + 2 * x;
  double const d = (2 * lgamma * lgamma * sqrt(2 * lgamma));
  double const expp1 = exp((x - leta) / lgamma) + 1;
  double const dexpp1 = -exp((x - leta) / lgamma) / lgamma;
  double const Result = -dexpp1 * cube(sqrt(y)) / (square(expp1) * d);

  Ensure(Result >= 0.0);
  return Result;
}

static double Feta_brute(double const eta, double const gamma) {
  // Partial degenerate: Sommerfeld expansion not sufficiently accurate.  Must
  // integrate explicitly.
  leta = eta;
  lgamma = gamma;
  double const max1 = (eta > 0 ? Feta_integrand(eta) : 0);
  double const max2 = Feta_integrand(1.5 * gamma);
  double const max3 = Feta_integrand(3 * gamma);
  double tol = numeric_limits<double>::epsilon() * max(max1, max(max2, max3)) *
               (max(eta, 0.0) + gamma);

  // help the compiler out by telling it that rkqs is a function pointer that
  // returns void and has the following argument list.  We have added this
  // typedef because cxx needs help parsing the call to quad(...).
  typedef void (*fpv)(std::vector<double> &, std::vector<double> const &,
                      double &, double, double, std::vector<double> const &,
                      double &, double &, Quad_To_ODE<double (*)(double)>);
  fpv rkqs_fpv = &rkqs<double, Quad_To_ODE<double (*)(double)>>;
  return quad(&Feta_integrand, 0.0, max(eta, 0.0) + 50 * gamma, tol, rkqs_fpv);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Evaluate the relativistic Fermi-Dirac integral
 *
 * The relativistic Fermi-Dirac integral is defined as
 * \f[
 * F_{3/2}(\eta, \gamma) = \frac{1}{2^{3/2}\gamma^{5/2}} \int_0^\infty
 * \frac{(x^2+2x)^{3/2}}{e^\frac{x-\eta}{\gamma}+1} dx
 * \f]
 * The dimensionless number density is its partial derivative with eta.
 *
 * \param eta[in] Dimensionless chemical potential \f$\eta=\frac{\mu}{kT}\f$
 * \param gamma[in] Dimensionless temperature \f$\gamma=\frac{kT}{mc^2}\f$
 *
 * \return Dimensionless number density \f$\frac{\partial
 *         F_{3/2}}{\partial\eta}=\frac{3Nh^3}{8\pi m^3c^3}\f$.
 */
double F_eta(double const eta, double const gamma) {
  Require(gamma > 0.0);

  double const TOL = 1.0e-5;

  double const e = exp(eta / gamma);
  double const de = e / gamma;
  if (e <= 0.5 && gamma < 0.2) {
    // Classical regime.  Expand in powers of exp(eta/gamma).
    double sum = 0.0;
    double dsum = 0.0;
    double fac = 1;
    for (int i = 0; i < 9; i++) {
      double si = 1;
      double dsi = 0;
      double ep = 1;
      double dep = 0;
      double sign = 1;
      for (int j = 2;; j++) {
        dep = dep * e + ep * dep;
        ep *= e;
        sign *= -1;
        double srt = pow((double)j, i + 1.5);
        double term = sign * ep / srt;
        if (fabs(term) < fabs(si) * std::numeric_limits<double>::epsilon())
          break;
        double dterm = sign * dep / srt;
        si += term;
        dsi += dterm;
      }
      double const term = 0.75 * sqrt(PI) * e * fac * si / factorial(i);
      if (fabs(term) < fabs(sum) * std::numeric_limits<double>::epsilon())
        break;
      double const dterm =
          0.75 * sqrt(PI) * fac * (de * si + e * dsi) / factorial(i);
      sum += term;
      dsum += dterm;
      fac *= 0.25 * gamma * (3 - 2 * i) * (2.5 + i);
    }
    return dsum;
  } else {
    // Degenerate regime?
    double const e = square(eta + 1) - 1;

    if (e <= 0) {
      return Feta_brute(eta, gamma);
    } else {
      double const de = 2 * (eta + 1);
      double const x = sqrt(e);
      double const dx = 0.5 * de / x;
      double const rad = sqrt(x * x + 1);
      double const drad = x * dx / rad;

      double n1, dn1;
      // Two versions of Sommerfeld expansion with different roundoff
      // properties.
      if (x < 0.07) {
        n1 = 0.2 * x * x * x * x * x *
             (1 - (5. / 14.) * x * x + (5. / 24.) * x * x * x * x -
              (25. / 176.) * x * x * x * x * x * x +
              (35. / 1024.) * x * x * x * x * x * x * x * x);
        dn1 = dx *
              (5 * n1 / x +
               0.2 * x * x * x * x * x * x *
                   (-5. / 7. + (5. / 6.) * x * x - (75. / 88.) * x * x * x * x +
                    (35. / 128.) * x * x * x * x * x * x));
      } else {
        n1 = (x * (2 * x * x - 3) * rad + 3 * log(rad + x)) / 8;
        dn1 = (dx * (2 * x * x - 3) * rad + 4 * x * x * dx * rad +
               x * (2 * x * x - 3) * drad + 3 * (drad + dx) / (rad + x)) /
              8;
      }
      double const n3 = 7 * square(square(PI * gamma)) * (2 * x * x - 1) * rad /
                        (120 * x * x * x);
      double const dn3 = 7 * square(square(PI * gamma / x)) * dx / (40 * rad);
      if (fabs(n3 / n1) < TOL) {
        // Sommerfeld expansion is sufficiently accurate
        double const dn2 =
            square(PI * gamma) * (2 * x * x + 1) * dx / (2 * rad);
        return (dn1 + dn2 + dn3) / (2 * gamma * gamma * sqrt(2 * gamma));
      } else {
        return Feta_brute(eta, gamma);
      }
    }
  }
}

} // end namespace rtt_sf

//---------------------------------------------------------------------------//
// end of F_eta.cc
//---------------------------------------------------------------------------//
