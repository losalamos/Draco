//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/qr_unpack_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of qr_unpack
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "qr_unpack.i.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

//---------------------------------------------------------------------------//
// RandomContainer = vector<double>
//---------------------------------------------------------------------------//

template DLL_PUBLIC_linear void qr_unpack(vector<double> &r, const unsigned n,
                                          const vector<double> &c,
                                          const vector<double> &d,
                                          vector<double> &qt);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
// end of qr_unpack.cc
//---------------------------------------------------------------------------//
