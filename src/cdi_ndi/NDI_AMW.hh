//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_AMW.hh
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AMW class definition.
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_AMW_hh
#define cdi_ndi_NDI_AMW_hh

#include "units/PhysicalConstexprs.hh"
#include <string>

namespace rtt_cdi_ndi {

//============================================================================//
/*!
 * \class NDI_AMW
 *
 * \brief Class for getting atomic mass weights by ZAID from NDI data using a
 *        method due to T. Saller. For more details on NDI, see
 *        https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *        Currently only multigroup data is supported, continuous energy data
 *        is probably best added through a refactor.
 */
//============================================================================//
class NDI_AMW {

public:
  //! Constructor (default gendir path)
  NDI_AMW();

  //! Constructor (overridden gendir path)
  NDI_AMW(const std::string gendir_path_in);

  //! Retrieve atomic mass weight for isotope with given ZAID
  double get_amw(const int zaid) const;

private:
  //! Path to gendir file
  std::string gendir_path;

  //! Unit system
  rtt_units::PhysicalConstexprs<rtt_units::CGS> pc;
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_AMW_hh

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_AMW.hh
//----------------------------------------------------------------------------//
