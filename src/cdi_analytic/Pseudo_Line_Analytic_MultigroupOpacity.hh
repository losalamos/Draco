//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.hh
 * \author Kent G. Budge
 * \date   Tue Apr  5 08:36:13 MDT 2011
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Pseudo_Line_Analytic_MultigroupOpacity_hh__
#define __cdi_analytic_Pseudo_Line_Analytic_MultigroupOpacity_hh__

#include "parser/Expression.hh"
#include "Analytic_MultigroupOpacity.hh"

namespace rtt_cdi_analytic
{
using rtt_dsxx::SP;
using rtt_parser::Expression;

//---------------------------------------------------------------------------//
/*!
 * \class Pseudo_Line_Analytic_Opacity_Model
 * \brief Derived Analytic_Opacity_Model class that defines a random line
 * spectrum for the opacity.
 *
 * The opacity function is a continuum on which is superimposed a number of
 * lines of the specified peak and width. The line locations are chosen at
 * random.
 *
 * The mass opacity coefficient is assumed independent of temperature or
 * density, which allows precalculation of the opacity structure, an important
 * time saver.
 *
 */
class Pseudo_Line_Analytic_MultigroupOpacity
    : public Analytic_MultigroupOpacity
{
  public:

    enum Averaging
    {
        NONE,      //!< evaluate opacity at band center
        ROSSELAND, //!< form a Rosseland (transparency) mean
        PLANCK,     //!< form a Planck (extinction) mean

        END_AVERAGING //!< sentinel value
    };
    
  private:
    // Coefficients
    SP<Expression const> continuum_;  // continuum opacity [cm^2/g]
    unsigned seed_;
    unsigned number_of_lines_;
    double line_peak_; // peak line opacity [cm^2/g]
    double line_width_;  // line width as fraction of line frequency.
    unsigned number_of_edges_;
    double edge_ratio_;
    Averaging averaging_;

    double Tref_;  // reference temperature for temperature dependence
    double Tpow_;  // temperature dependence exponent

    sf_double center_; // line centers for this realization
    sf_double edge_;   // edges for this realization
    sf_double edge_factor_; // opacity at threshold

    friend class PLR_Functor; // used in calculation of Rosseland averages
    friend class PLP_Functor; // used in calculation of Planck averages

  public:

    Pseudo_Line_Analytic_MultigroupOpacity(sf_double const &group_bounds,
                                           rtt_cdi::Reaction reaction,
                                           SP<Expression const> const &cont,
                                           unsigned number_of_lines,
                                           double line_peak,
                                           double line_width,
                                           unsigned number_of_edges,
                                           double edge_ratio,
                                           double Tref,
                                           double Tpow,
                                           double emin,
                                           double emax,
                                           Averaging averaging,
                                           unsigned seed);
    
    //! Constructor for packed state.
    explicit  Pseudo_Line_Analytic_MultigroupOpacity(const sf_char &packed);

    // Get the group opacities.
    virtual sf_double getOpacity(double, double) const;

    // Get the group opacity fields given a field of temperatures.
    virtual vf_double getOpacity(const sf_double &, double) const;

    // Get the group opacity fields given a field of densities.
    virtual vf_double getOpacity(double, const sf_double &) const;

    // Get the data description of the opacity.
    virtual std_string getDataDescriptor() const;

    //! Pack up the class for persistence.
    sf_char pack() const;

    //! Compute a monochromatic opacity
    double monoOpacity(double nu, double T) const;
};

} // end namespace rtt_cdi_analytic

#endif  // __cdi_analytic_Pseudo_Line_Analytic_MultigroupOpacity_hh__

//---------------------------------------------------------------------------//
//       end of cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.hh
//---------------------------------------------------------------------------//
