//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Comm_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Thu Jan  3 17:00:09 2002
 * \brief  Comm_Builder member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Comm_Builder.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// BUILD COMMUNICATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Build a Communicator object using the rtt_mc::Topology.

 * This is the primary service of the Comm_Builder class.  Given a valid
 * rtt_mc::Topology, this function returns a rtt_dsxx::SP to a Communicator.
 * If the Topology is full replication then a Communicator is unnecessary and
 * a null SP is returned.  For full DD and DD/replication topologies a fully
 * defined Communicator is returned.

 * \param topology smart pointer to a rtt_mc::Topology

 * \return smart pointer (null if the topology is replication) to a
 * Communicator

 */
template<class PT>
Comm_Builder<PT>::SP_Communicator 
Comm_Builder<PT>::build_Communicator(SP_Topology topology)
{
    Require (topology);

    // return communicator argument
    SP_Communicator communicator;

    // do work depending on the topology
    if (topology->get_parallel_scheme() == "replication")
    {
	Insist (!communicator, "Defined communicator for replication!");
    }
    else if (topology->get_parallel_scheme() == "DD")
    {
        communicator = build_DD_Comm(topology);
	Ensure (communicator);
    }
    else if (topology->get_parallel_scheme() == "DD/replication")
    {
	Insist (0, "Don't support DD/replication yet!");
    }
    else 
    {
	Insist (0, "Invalid topology given in Comm_Builder!");
    }

    // return the communicator
    return communicator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Deprecated build function for when the boundary nodes and cells are
 * known.

 * This function builds a communicator when the boundary nodes and boundary
 * cells are known as simple 1-D vectors (with one-to-one matching).  This
 * can be used as a builder to hosts that provide this information in a flat
 * interface.  Additionally, this function \b only works for full DD
 * topologies.

 * Clients should use the build_Communicator(SP_Topology) instead.

 */
template<class PT> Comm_Builder<PT>::SP_Communicator 
Comm_Builder<PT>::build_Communicator(const sf_int &b_node,
				     const sf_int &b_cell)
{
    using std::vector;

    // requirements
    Require (b_node.size() == b_cell.size());
    
    // return communicator
    SP_Communicator return_com;

    // make objects needed by the Communicator
    int num_cells = b_node.size();
    sf_int com_nodes;
    vf_int boundary_node(num_cells);
    vf_int boundary_cell(num_cells);
    vector<bool> procs(C4::nodes(), false);

    // now assign to boundary_data
    for (int i = 0; i < num_cells; i++)
    {
	// dimension the number of processors that each boundary cell resides 
	// on, because this is a full DD problem, each boundary cell lives on 
	// only one processor
	boundary_node[i].resize(1);
	boundary_cell[i].resize(1);

	// get the processor and local cell on processor for a boundary cell
	int send_proc  = b_node[i];
	int local_cell = b_cell[i];
	Check (send_proc < C4::nodes());
	
	// assign the data
	boundary_node[i][0] = send_proc;
	boundary_cell[i][0] = local_cell;

	// let this processor now that it communicates with send_proc
	procs[send_proc] = true;
    }

    // loop over the nodes on the problem and see if we communicate with it-- 
    // lets reserve some memory in the com_nodes vector for efficiency,
    // following Lippman, we reserve only one space for performance
    com_nodes.reserve(1);
    for (int n = 0; n < C4::nodes(); n++)
	if (procs[n])
	    com_nodes.push_back(n);

    // return the communicator, because we are in the DD realm, we have
    // one-to-one communication
    return_com = new Communicator<PT>(com_nodes, com_nodes, 
				      boundary_node, boundary_cell); 

    // some final assurances
    Ensure (boundary_node.size() == boundary_cell.size());
    Ensure (return_com);

    return return_com;
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION 
//---------------------------------------------------------------------------//
/*!
 * \brief Build a communicator in DD topologies.
 */
template<class PT>
Comm_Builder<PT>::SP_Communicator 
Comm_Builder<PT>::build_DD_Comm(SP_Topology topology)
{
    using std::vector;

    Require (topology->get_parallel_scheme() == "DD");

    // return communicator
    SP_Communicator communicator;

    // number of boundary cells
    int num_cells = topology->get_boundary_cells(C4::node());

    // make data fields that will be required for this communicator

    // list of nodes that this processor sends and receives from; because
    // this is a full DD communicator the receive and send nodes data is
    // identical; there is one-to-one communication
    sf_int com_nodes;

    // boundary cell and node fields, each boundary cell's real cell resides
    // on only one processor
    vf_int bnd_node(num_cells, vector<int>(1));
    vf_int bnd_cell(num_cells, vector<int>(1));

    // list of processors to determine if we communicate with them or not
    vector<bool> procs(C4::nodes(), false);

    // loop through boundary cells and determine what the local cell is on
    // the send-to processor
    int processor;
    int local_cell;
    int global_cell;
    for (int bc = 1; bc <= num_cells; bc++)
    {
	// determine the global cell 
	global_cell = topology->boundary_to_global(bc, C4::node());
	
	// determine the processor it lives on, remember each global cell
	// only lives on one processor in full DD
	processor = topology->get_procs(global_cell).front();
	Check (processor >= 0 && processor < C4::nodes());
	Check (processor != C4::node());

	// determine the local cell index of the global cell on that
	// processor
	local_cell = topology->local_cell(global_cell, processor);
	Check (local_cell);

	// indicate that we communicate with this processor
	procs[processor] = true;

	// fill up the boundary node and cell data
	bnd_node[bc-1][0] = processor;
	bnd_cell[bc-1][0] = local_cell;
    }

    // loop over the nodes on the problem and see if we communicate with it-- 
    // lets reserve some memory in the com_nodes vector for efficiency,
    // following Lippman, we reserve only one space for performance
    com_nodes.reserve(1);
    for (int n = 0; n < C4::nodes(); n++)
	if (procs[n])
	    com_nodes.push_back(n);

    // build and return communicator, remember the recv and send node data is 
    // the same for full DD because of one-to-one communication
    communicator = new Communicator<PT>(com_nodes, com_nodes, bnd_node, 
					bnd_cell);
    Ensure (communicator);
    return communicator;
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                        end of mc/Comm_Builder.t.hh
//---------------------------------------------------------------------------//
