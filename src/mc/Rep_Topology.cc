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

//===========================================================================//
// REP_TOPOLOGY CLASS DEFINITIONS
//===========================================================================//

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
// PACK
//---------------------------------------------------------------------------//
/*!

 * \brief Pack up the Rep_Topology and return a smart pointer to a
 * Topology::Pack base class.

 * Data in the Rep_Topology is packed into a Topology::Pack base class
 * pointer.  The actual data resides in the Rep_Topology::Pack derived class.
 * The pointer can be used to access pointer to the data for message passing
 * and output operations.

 */
Topology::SP_Pack Rep_Topology::pack() const
{
    // we are return a smart pointer to the base class pack object
    Topology::SP_Pack return_pack;
    
    // allocate a pointer the size of the global mesh
    int *data = new int(global_cells);

    // build the pack object, no need to clean up the allocated pointer, the
    // pack object will do it in its destructor
    return_pack = new Rep_Topology::Pack(data);
    Ensure (*return_pack->begin() == global_cells);

    return return_pack;
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

//===========================================================================//
// REP_TOPOLOGY::PACK CLASS DEFINITIONS
//===========================================================================//
/*!
 * \brief Unpack a Rep_Topology::Pack object into a smart pointer to a
 * rtt_mc::Topology.

 * The unpack function creates a smart pointer (rtt_dsxx::SP) to the Topology
 * base class from the data stored inside the Rep_Topology::Pack object.
 * Obviously, the Rep_Topology::Pack class creates a Topology base class
 * pointer to a Rep_Topology object.

 * \return rtt_dsxx::SP to a Topology base class pointer that points to a
 * Rep_Topology derived class instance

 */
Topology::SP_Topology Rep_Topology::Pack::unpack() const
{
    // base class smart pointer
    SP_Topology topology;

    // build the replication topology
    topology = new Rep_Topology(*data);
    Ensure (topology);

    return topology;
}

} // end of rtt_mc

//---------------------------------------------------------------------------//
//                              end of Rep_Topology.cc
//---------------------------------------------------------------------------//
