//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/nGray_Analytic_Odfmg_Opacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  nGray_Analytic_Odfmg_Opacity class definition.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_nGray_Analytic_Odfmg_Opacity_hh__
#define __cdi_analytic_nGray_Analytic_Odfmg_Opacity_hh__

#include "Analytic_Odfmg_Opacity.hh"

namespace rtt_cdi_analytic {

//===========================================================================//
/*!
 * \class nGray_Analytic_Odfmg_Opacity
 *
 * \brief Derived rtt_cdi::OdfmgOpacity class for analytic opacities.
 *
 * Primarily code from Analytic_Multigroup_Opacity.
 */
//===========================================================================//

class nGray_Analytic_Odfmg_Opacity : public Analytic_Odfmg_Opacity {
public:
  // Useful typedefs.
  typedef std::shared_ptr<Analytic_Opacity_Model> SP_Analytic_Model;
  typedef std::shared_ptr<const Analytic_Opacity_Model> const_Model;
  typedef std::vector<SP_Analytic_Model> sf_Analytic_Model;
  typedef std::vector<double> sf_double;
  typedef std::vector<sf_double> vf_double;
  typedef std::string std_string;
  typedef std::vector<char> sf_char;

private:
  // Analytic models for each group.
  sf_Analytic_Model group_models;

public:
  // Constructor.
  nGray_Analytic_Odfmg_Opacity(const sf_double &groups, const sf_double &bands,
                               const sf_Analytic_Model &models,
                               rtt_cdi::Reaction reaction_in,
                               rtt_cdi::Model model_in = rtt_cdi::ANALYTIC);

  // Constructor for packed nGray_Analytic_Odfmg_Opacities
  explicit nGray_Analytic_Odfmg_Opacity(const sf_char &);

  // >>> ACCESSORS
  const_Model get_Analytic_Model(int g) const { return group_models[g - 1]; }

  //! right now, all bands have same model (same opacity)
  const_Model get_Analytic_Model(int g, int /*b*/) const {
    return group_models[g - 1];
  }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::OdfmgOpacity

  /*!
   * \brief Opacity accessor that returns a 2-D vector of opacities
   *        (groups*bands) that correspond to the provided temperature and
   *        density.
   *
   * \param targetTemperature The temperature value for which an opacity value
   *             is being requested.
   * \param targetDensity The density value for which an opacity value is being
   *             requested.
   * \return A vector of opacities.
   */
  std::vector<std::vector<double>> getOpacity(double targetTemperature,
                                              double targetDensity) const;

  std::vector<std::vector<std::vector<double>>>
  getOpacity(const std::vector<double> &targetTemperature,
             double targetDensity) const;

  std::vector<std::vector<std::vector<double>>>
  getOpacity(double targetTemperature,
             const std::vector<double> &targetDensity) const;

  // Get the data description of the opacity.
  inline std_string getDataDescriptor() const;

  // Pack the nGray_Analytic_Odfmg_Opacity into a character string.
  sf_char pack() const;
};

//---------------------------------------------------------------------------//
//! Return a string describing the opacity model.
nGray_Analytic_Odfmg_Opacity::std_string
nGray_Analytic_Odfmg_Opacity::getDataDescriptor() const {
  std_string descriptor;

  rtt_cdi::Reaction const rxn = getReactionType();

  if (rxn == rtt_cdi::TOTAL)
    descriptor = "Analytic Odfmg Total";
  else if (rxn == rtt_cdi::ABSORPTION)
    descriptor = "Analytic Odfmg Absorption";
  else if (rxn == rtt_cdi::SCATTERING)
    descriptor = "Analytic Odfmg Scattering";
  else
    Insist(0, "Invalid analytic multigroup model opacity!");

  return descriptor;
}

} // end namespace rtt_cdi_analytic

#endif // __cdi_analytic_nGray_Analytic_Odfmg_Opacity_hh__

//---------------------------------------------------------------------------//
// end of cdi_analytic/nGray_Analytic_Odfmg_Opacity.hh
//---------------------------------------------------------------------------//
