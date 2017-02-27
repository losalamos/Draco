//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Subrandom_Generator.cc
 * \author Kent Budge
 * \date   Thu Dec 22 13:38:35 2005
 * \brief
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Subrandom_Generator.hh"
#include "ds++/Assert.hh"

namespace rtt_rng {

//---------------------------------------------------------------------------//
/*!
 * \param count Dimension of the vector calculated by this generator.
 */

Subrandom_Generator::Subrandom_Generator(unsigned const count)
    : count_(count), element_(0) {}

} // end namespace rtt_rng

//---------------------------------------------------------------------------//
//                 end of Subrandom_Generator.cc
//---------------------------------------------------------------------------//
