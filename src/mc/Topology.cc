//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Topology.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 17 14:56:26 1999
 * \brief  Topology class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Topology.hh"
#include "ds++/Assert.hh"

namespace rtt_mc
{

using std::ostream;

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*! 
 * \brief Topology abstract base class constructor.
 *
 * \param scheme string stating the parallel scheme, must be "replication",
 * "DD", or "DD/replication"
*/
Topology::Topology(const std_string &scheme)
    : parallel_scheme(scheme)
{
    Ensure (parallel_scheme == "replication" ||
	    parallel_scheme == "DD" ||
	    parallel_scheme == "DD/replication");
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

ostream& operator<<(ostream &out, const Topology &object)
{
    object.print(out);
    return out;
}

} // end of rtt_mc

//---------------------------------------------------------------------------//
//                              end of Topology.cc
//---------------------------------------------------------------------------//
