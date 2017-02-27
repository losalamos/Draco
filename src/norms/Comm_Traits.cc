//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/Comm_Traits.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 12:45:49 2005
 * \brief  Implementations for Comm_Traits.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Comm_Traits.hh"

namespace rtt_norms {

// Index_Labeled send/receive

void Comm_Traits<Index_Labeled>::send(const Index_Labeled &x, const size_t n) {
  rtt_c4::send(&(x.index), 1, n);
  rtt_c4::send(&(x.processor), 1, n);

  size_t const s = x.label.size();
  rtt_c4::send(&s, 1, n);
  rtt_c4::send(x.label.c_str(), s, n);
}

void Comm_Traits<Index_Labeled>::receive(Index_Labeled &x, const size_t n) {
  rtt_c4::receive(&(x.index), 1, n);
  rtt_c4::receive(&(x.processor), 1, n);

  size_t s(0);
  rtt_c4::receive(&s, 1, n);
  x.label.resize(s);
  rtt_c4::receive(&(x.label[0]), s, n);
}

// Index_Proc send/receive

void Comm_Traits<Index_Proc>::send(const Index_Proc &x, const size_t n) {
  rtt_c4::send(&(x.index), 1, n);
  rtt_c4::send(&(x.processor), 1, n);
}

void Comm_Traits<Index_Proc>::receive(Index_Proc &x, const size_t n) {
  rtt_c4::receive(&(x.index), 1, n);
  rtt_c4::receive(&(x.processor), 1, n);
}

} // end namespace rtt_norms

//---------------------------------------------------------------------------//
//              end of norms/Comm_Traits.hh
//---------------------------------------------------------------------------//
