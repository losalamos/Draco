//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/global_containers_pt.cc
 * \author Kent Budge
 * \date   Mon Mar 24 10:17:40 2008
 * \brief  Explicit template instatiations for class global_containers.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 * 
 * 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4/config.h"
#ifdef C4_MPI
#include "global_containers.i.hh"

namespace rtt_c4
{
using namespace std;

template
void global_merge(set<unsigned> &);

template
void global_merge(map<unsigned, double> &);

template
void global_merge(map<unsigned, bool> &);

} // end namespace rtt_c4
#endif // C4_MPI 
//---------------------------------------------------------------------------//
//                 end of global_containers_pt.cc
//---------------------------------------------------------------------------//
