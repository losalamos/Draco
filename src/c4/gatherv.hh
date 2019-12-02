//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/gatherv.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:42:03 2002
 * \brief  Data gatherv functions
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 *
 * This file contains the declarations for determinate and indeterminate
 * variable length gather functions.
 */
//---------------------------------------------------------------------------//

#ifndef c4_gatherv_hh
#define c4_gatherv_hh

#include "C4_Traits.hh"
#include "C4_sys_times.h"
#include <vector>

namespace rtt_c4 {
//---------------------------------------------------------------------------//
/*!
 * \brief Gather messages of known but processor-dependent size
 *
 * This subroutine handles the case where the lengths of each processor's
 * message are known in advance.
 *
 * \param outgoing_data Data to be send to root processor.
 *
 * \param incoming_data Ignored on any processor but the root processor. On
 * the root processor, the size of each subarray must be set to the expected
 * size of the incoming message. On return, contains the gathered data.
 */
template <class T>
void determinate_gatherv(std::vector<T> &outgoing_data,
                         std::vector<std::vector<T>> &incoming_data);

//---------------------------------------------------------------------------//
/*!
 * \brief Gather messages of unknown size
 *
 * This subroutine handles the case where the lengths of each processor's
 * message are not known in advance.
 *
 * \param outgoing_data Data to be send to root processor.
 *
 * \param incoming_data Ignored on any processor but the root processor. On
 * the root processor, on return, contains the gathered data.
 */
template <class T>
void indeterminate_gatherv(std::vector<T> &outgoing_data,
                           std::vector<std::vector<T>> &incoming_data);

void indeterminate_gatherv(std::string &outgoing_data,
                           std::vector<std::string> &incoming_data);

} // end namespace rtt_c4

#endif // c4_gatherv_hh

//---------------------------------------------------------------------------//
// end of c4/gatherv.hh
//---------------------------------------------------------------------------//
