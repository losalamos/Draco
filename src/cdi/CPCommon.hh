//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CPCommon.hh
 * \author Kendra Long
 * \date   Mon Jan 19 13:41:01 2001
 * \brief  Datatypes needed in CPEloss classes
 * \note   Copyright (C) 2016-2019 Los Alamos National Securty, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_CPCommon_hh__
#define __cdi_CPCommon_hh__

#include <cstdint>

namespace rtt_cdi {

namespace constants {
//! Number of models contained in rtt_cdi::CPModel.
unsigned int const num_CPModels(3);
} // namespace constants

enum class CPModel {
  ELOSS = 0,        /*!< energy loss model */
  LA_SCATTERING = 1 /*!< large-angle scattering model */
};

enum class CPModelType {
  UNASSIGNED_ETYPE =
      0, /*!< unassigned type; used as a placeholder before deciding type */
  ANALYTIC_ETYPE = 1, /*!< an Analytic eloss model */
  TABULAR_ETYPE = 2   /*!< a Tabular eloss model */
};

// Basic data required to describe a charged particle
class CParticle {
private:
  int32_t zaid;
  double mass;

public:
  int32_t get_zaid() const { return zaid; }
  double get_mass() const { return mass; }

  CParticle(const int32_t zaid_in, const double mass_in)
      : zaid(zaid_in), mass(mass_in) {}
};

} // namespace rtt_cdi

#endif
