//----------------------------------*-C++-*----------------------------------//
// tstTally_pt.cc
// Thomas M. Evans
// Wed Apr 28 16:09:30 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Instantiation of Tally
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Tally.t.hh"
#include <iostream>

namespace rtt_imc
{

using rtt_mc::OS_Mesh;
using std::iostream;

template class Tally<OS_Mesh>;
template ostream& operator<<(ostream &, const Tally<OS_Mesh> &);

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstTally_pt.cc
//---------------------------------------------------------------------------//
