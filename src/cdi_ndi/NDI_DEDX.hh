//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/Tabular_CP_Eloss.hh
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  Tabular_CP_Eloss class definition.
 * \note   Copyright (C) 2019-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef cdi_ndi_Tabular_CP_Eloss_hh
#define cdi_ndi_Tabular_CP_Eloss_hh

namespace rtt_cdi_ndi {

//============================================================================//
/*!
 * \class Tabular_CP_Eloss
 *
 * \brief Derived rtt_cdi::CPEloss class for tabular eloss.
 *        This class implements the interface found in cdi/CPEloss.hh for
 *        the case where CP energy loss data is in tabular form, retrieved
 *        from NDI.
 */
//============================================================================//

class Tabular_CP_Eloss : public rtt_cdi::CPEloss {

}

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_Tabular_CP_Eloss_hh

//----------------------------------------------------------------------------//
// End cdi_ndi/Tabular_CP_Eloss.hh
//----------------------------------------------------------------------------//
