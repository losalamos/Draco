//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Rep_Topology.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 30 17:01:44 1999
 * \brief  Rep_Topology implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Rep_Topology.hh"

namespace rtt_mc
{

using std::cout;
using std::endl;
using std::ostream;

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Rep_Topology constructor.
 *
 * The Rep_Topology constructor only requires the number of global cells in
 * the problem to build an instance of Rep_Topology.  All other functions are
 * trivially derived based on the simplicity of a full replication topology.
 *
 * \param num_global_cells number of global cells in the problem 
 */
Rep_Topology::Rep_Topology(int num_global_cells)
    : Topology("replication"), global_cells(num_global_cells)
{
    Ensure ( get_parallel_scheme() == "replication");
}

//---------------------------------------------------------------------------//
// PRINT FOR DIAGNOSTICS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out data in the Rep_Topology class for diagnostic purposes.
 */
void Rep_Topology::print(ostream &out) const
{
    out << endl;
    out << "*** REP_TOPOLOGY DATA ***" << endl;
    out << "-------------------------" << endl;

    out << endl;

    out << "Parallel Topology: " << get_parallel_scheme() << endl;
    out << "-------------------------" << endl;
}

} // end of rtt_mc

//---------------------------------------------------------------------------//
//                              end of Rep_Topology.cc
//---------------------------------------------------------------------------//
