//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.hh
 * \author Mike Buksas, Rob Lowrie
 * \date   Mon Nov 19 10:09:10 2007
 * \brief  Implementation of Invert_Comm_Map
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#ifndef c4_Invert_Comm_Map_hh
#define c4_Invert_Comm_Map_hh

#include "C4_Functions.hh"
#include <map>

namespace rtt_c4 {

//! Map type for invert_comm_map
typedef std::map<int, size_t> Invert_Comm_Map_t;

//---------------------------------------------------------------------------//
/**
 * \brief Invert the contents of a one-to-many mapping between nodes.
 *
 * \param[in] to_map On input, a map from processor number to the size of
 *        information to be sent to (or received from) that processor
 *        by the current processor.
 * \param[out] from_map On output, a map from processor number to the size of
 *        information to be received from (or sent to) that processor by
 *        the current processor.  On input, ignored and deleted.
 *
 * Here, the units of the "size of information" is up to the caller.  For
 * example, it might be the number of bytes, or the number of elements in an
 * array.  The size must be positive (specifically, nonzero).
 */
DLL_PUBLIC_c4 void invert_comm_map(Invert_Comm_Map_t const &to_map,
                                   Invert_Comm_Map_t &from_map);

} // end namespace rtt_c4

#endif // c4_Invert_Comm_Map_hh

//---------------------------------------------------------------------------//
// end of c4/Invert_Comm_Map.hh
//---------------------------------------------------------------------------//
