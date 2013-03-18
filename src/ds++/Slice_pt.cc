//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   container/Slice.cc
 * \author Kent Budge
 * \date   Thu Jul  8 08:06:53 2004
 * \brief  Specializations of Slice template
 * \note   Copyright (C) 2004-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Slice.hh"
#include <vector>

namespace rtt_dsxx
{

template<> class Slice<double*>;
template<> class Slice<std::vector<double>::iterator>;

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of container/Slice.cc
//---------------------------------------------------------------------------//
