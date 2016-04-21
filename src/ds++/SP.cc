//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/SP.cc
 * \author Kent G. Budge
 * \date   Wed Jan  4 16:51:37 2012
 * \brief  Explicit template instatiations for class SP.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "SP.hh"

namespace rtt_dsxx
{

void incompatible(std::type_info const &X, std::type_info const &T)
{
    std::string msg =
        std::string("Incompatible dumb pointer conversion between ")
        + X.name() + " and SP<" + T.name() + ">.";

    Insist(false, msg.c_str());
}

// Explicit template instantiations go here.

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of SP_pt.cc
//---------------------------------------------------------------------------//
