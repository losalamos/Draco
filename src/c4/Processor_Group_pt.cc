//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Processor_Group_pt.hh
 * \author Kent Budge
 * \date   Fri Oct 20 13:49:10 2006
 * \brief  Member definitions of class Processor_Group
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4/config.h"
#ifdef C4_MPI

#include "Processor_Group.hh"
#include "Processor_Group.i.hh"

namespace rtt_c4
{
using namespace std;

template void Processor_Group::sum(vector<double> &values);

} // end namespace rtt_c4

#endif  //C4_MPI

//---------------------------------------------------------------------------//
//              end of c4/Processor_Group_pt.hh
//---------------------------------------------------------------------------//
