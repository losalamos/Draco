//----------------------------------*-C++-*----------------------------------//
// Comm_Builder.t.hh
// Thomas M. Evans
// Tue Jun  1 17:01:44 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Comm_Builder implementation
//---------------------------------------------------------------------------//

#include "Comm_Builder.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"

namespace rtt_imc
{

using std::vector;

//---------------------------------------------------------------------------//
// Build Communicators
//---------------------------------------------------------------------------//
// build a communicator for full dd stuff

template<class PT> Comm_Builder<PT>::SP_Comm 
Comm_Builder<PT>::build_Communicator(const intvec &b_node,
				     const intvec &b_cell)
{
    // requirements
    Require (b_node.size() == b_cell.size());
    
    // return communicator
    SP_Comm return_com;

    // make objects needed by the Communicator
    int num_cells = b_node.size();
    intvec com_nodes;
    dintvec boundary_node(num_cells);
    dintvec boundary_cell(num_cells);
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
// build a general DD/replication communicator

template<class PT>
Comm_Builder<PT>::SP_Comm Comm_Builder<PT>::build_Communicator()
{
    Insist(0, "Have not implemented this service yet!");
    return SP_Comm();
}

}  // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Comm_Builder.t.hh
//---------------------------------------------------------------------------//
