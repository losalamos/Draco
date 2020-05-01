//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_TNReaction.hh
 * \author Ben R. Ryan
 * \date   2019 Nov 4
 * \brief  NDI_TNReaction class definition.
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_TNReaction_hh
#define cdi_ndi_NDI_TNReaction_hh

#include "NDI_Base.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>
#include <memory>
#include <string>

namespace rtt_cdi_ndi {

//============================================================================//
/*!
 * \class NDI_TNReaction
 *
 * \brief Class for wrapping NDI calls to NDI tn data. Reads data, constructs
 *        internal storage amenable to radiation calculations, and provides
 *        accessors.
 * \example cdi_ndi/test/tstNDI_TNReaction.cc
 */
//============================================================================//

class NDI_TNReaction : public NDI_Base {

public:
  //! Constructor (default gendir path)
  NDI_TNReaction(const std::string &library_in, const std::string &reaction_in,
                 const std::vector<double> mg_e_bounds_in);

  //! Constructor (overridden gendir path)
  NDI_TNReaction(const std::string &gendir_in, const std::string &library_in,
                 const std::string &reaction_in,
                 const std::vector<double> mg_e_bounds_in);

  //! Disable default constructor
  NDI_TNReaction() = delete;

  //! Disable copy constructor (meaning no implicit move assignment operator
  //! or move constructor)
  NDI_TNReaction(const NDI_TNReaction &) = delete;

  //! Return spectrum PDF at a given temperature
  std::vector<double> get_PDF(const int product_zaid,
                              const double temperature) const;

private:
  void load_ndi();
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_TNReaction_hh

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_TNReaction.hh
//----------------------------------------------------------------------------//
