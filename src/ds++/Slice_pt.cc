//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   container/Slice.cc
 * \author Kent Budge
 * \date   Thu Jul  8 08:06:53 2004
 * \brief  Specializations of Slice template
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>

#include "Slice.hh"

namespace rtt_dsxx
{

template<> class Slice<double*>;
template<> class Slice<std::vector<double>::iterator>;

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
//              end of container/Slice.cc
//---------------------------------------------------------------------------//
