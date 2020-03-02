//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_TN.hh
 * \author Ben R. Ryan
 * \date   2019 Nov 4
 * \brief  NDI_TN class definition.
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_TN_hh
#define cdi_ndi_NDI_TN_hh

#include "NDI_Base.hh"
#include "ndi.h"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>
#include <memory>
#include <string>

namespace rtt_cdi_ndi {

//============================================================================//
/*!
 * \class NDI_TN
 *
 * \brief Class for wrapping NDI calls to NDI tn data. Reads data, constructs
 *        internal storage amenable to radiation calculations, and provides
 *        accessors.
 */
//============================================================================//

class NDI_TN : public NDI_Base {

public:
  //! Constructor
  NDI_TN(const std::string &gendir_in, const std::string &library_in,
         const std::string &reaction_in, const MG_FORM mg_form_in);

  //! Disable default constructor
  NDI_TN() = delete;

  //! Disable copy constructor (meaning no implicit move assignment operator
  //! or move constructor)
  NDI_TN(const NDI_TN &) = delete;

  //! Return spectrum PDF at a given temperature
  std::vector<double> get_PDF(const int product_zaid,
                              const double temperature) const;
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_TN_hh

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_TN.hh
//----------------------------------------------------------------------------//
