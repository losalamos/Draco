//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shared_lib/test/Foo_Base.cc
 * \author Rob Lowrie
 * \date   Wed Dec 29 11:24:17 2004
 * \brief  Implementatin of Foo_Base.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#include "Foo_Base.hh"

namespace rtt_shared_lib_test {

// This should be linked in with the executable, NOT the shared lib.
Foo_Base::~Foo_Base() {}

} // end namespace rtt_shared_lib_test

//----------------------------------------------------------------------------//
//              end of shared_lib/test/Foo_Base.cc
//----------------------------------------------------------------------------//
