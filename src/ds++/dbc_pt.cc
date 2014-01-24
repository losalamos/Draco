//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/dbc_pt.cc
 * \author Kent Budge
 * \date   Tue Feb 19 14:28:59 2008
 * \brief  Explicit template instatiations for class dbc.
 * \note   Copyright (C) 2008-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>
#include "dbc.i.hh"

namespace rtt_dsxx
{
using namespace std;

template
bool is_strict_monotonic_increasing(vector<double>::iterator first,
                                    vector<double>::iterator last);

template
bool is_strict_monotonic_increasing(vector<double>::const_iterator first,
                                    vector<double>::const_iterator last);


} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of dbc_pt.cc
//---------------------------------------------------------------------------//
