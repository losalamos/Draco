//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Topology_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Mon Nov 29 14:33:42 1999
 * \brief  Topology_Builder template implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Topology_Builder.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"

#include <algorithm>

namespace rtt_imc
{

// namespace declarations
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using dsxx::SP;
using std::vector;
using std::string;
using std::find;

//---------------------------------------------------------------------------//
// Constructors
//---------------------------------------------------------------------------//
/*!
 * \brief Serial Topology_Builder constructor.
 *
 * The Topology_Builder constructor makes a Topology_Builder object that can
 * be used to build rtt_mc::Topology objects.  The Topology_Builder
 * constructor sets one piece of data, the processor per cell capacity.  It
 * gets this information from the interface object which is required to have
 * a get_capacity() member function.  The capacity must be greater than
 * zero. 
 *
 * This constructor requires calling on the Host node (node 0).  Thus, it can
 * be thought of as the "serial" constructor for the Topology_Builder class.
 *
 * This function must be explicitly instantiated with the appropriate
 * interfaces that are to be used in a given problem.
 *
 * \param interface Smart Pointer (dsxx::SP) to an appropriate interface
 * class.  
 */
template<class MT>
template<class IT>
Topology_Builder<MT>::Topology_Builder(SP<IT> interface)
    : capacity()
{
    Require (C4::node() == 0);

    // get the processor/cell capacity from the interface
    capacity = interface->get_capacity();

    Ensure (capacity > 0);
}

//---------------------------------------------------------------------------//
// Build Topology Functions
//---------------------------------------------------------------------------//
/*!
 * \brief Serial topology building function.
 *
 * This function takes a global instantiation of the mesh and does a parallel
 * topology reconstruction based on the cell/processor capacity.  This
 * function either builds a replication topology (all cells replicated on all
 * meshes) or a full DD topology (each cell uniquely placed on a processor).
 *
 * This function can only be run on the host processor (node 0); thus, it is
 * the "serial" topology build function.
 *
 * The parallel topology is built in one of two ways:
 * \arg \b replication: capacity > number of cells in the mesh;
 * \arg \b DD: capacity * nodes = number of cells in the mesh.
 * More advanced topologies require more input than the capacity and will be
 * added at a later date.  If the capacity * nodes > number of cells in the
 * mesh, and capacity < number of cells in the mesh a failure is noted.  This
 * situation indicates a DD/replication problem that will be dealt with at a
 * later date.
 *
 * \param mesh smart pointer (dsxx::SP) to a \b global mesh
 * \return smart pointer (dsxx::SP) to a rtt_mc::Topology 
 */
template<class MT>
SP<Topology> Topology_Builder<MT>::build_Topology(SP<MT> mesh)
{
    // this function can only be run on the host node
    Require (C4::node() == 0);

    // make sure we have a mesh and that it is not a submesh
    Require (mesh);
    Require (mesh->full_Mesh());

    // number of cells on the global mesh
    int num_cells = mesh->num_cells();

    // calculate the total capacity
    int total_capacity = capacity * C4::nodes();
    Check (total_capacity >= num_cells);

    // smart pointer to return topology
    SP<Topology> return_top;

    // do full replication / full DD / failure (if rep-DD)
    if (num_cells <= capacity)
    {
	// build a full replication topology
	return_top = new Rep_Topology(num_cells);
    }
    else if (num_cells == total_capacity)
    {
	// data necessary for a General_Topology instance
	vf_int cells_per_proc(C4::nodes());
	vf_int procs_per_cell(num_cells);
	vf_int bound_cells;
	string parallel_scheme;

	// build DD topology data
	build_DD(cells_per_proc, procs_per_cell, bound_cells,
		 parallel_scheme, mesh); 

	// build DD topology
	return_top = new General_Topology(cells_per_proc, procs_per_cell,
					  bound_cells, parallel_scheme);
    }
    else
    {
	// we don't handle the DD/replication case here yet
	Insist(0, "Don't handle DD/replication yet");
    }

    // return the Topology
    Ensure (return_top);
    return return_top;
}

//---------------------------------------------------------------------------//
// Private Functions
//---------------------------------------------------------------------------//
/*!
 * \brief Build full DD topology data.
 *
 * \param mesh dsxx::SP to a mesh, this must be a global mesh.
 */
template<class MT>
void Topology_Builder<MT>::build_DD(vf_int &cells_per_proc, 
				    vf_int &procs_per_cell, 
				    vf_int &bound_cells,
				    string &parallel_scheme,
				    SP<MT> mesh)
{
    // get the number of cells 
    int num_cells = procs_per_cell.size();
    Require (num_cells > 0);
    Require (num_cells == mesh->num_cells());
    Require (capacity * C4::nodes() == num_cells);
    Require (mesh->full_Mesh());

    // set the parallel_scheme
    parallel_scheme = "DD";

    // sweep through processors and parse out cells to each processor until
    // the processor capacity is reached.  This form of DD is simple, cells
    // are parsed in numerical order to be placed on a mesh.
    
    // start at cell 1
    int cell = 1;
    for (int proc = 0; proc < C4::nodes(); proc++)
    {
	// on-processor cell counter
	int local_cell = 1;
	while (local_cell <= capacity)
	{
	    procs_per_cell[cell-1].push_back(proc);
	    cells_per_proc[proc].push_back(cell);
	    Ensure (procs_per_cell[cell-1].size() == 1);

	    // update the cell and local cell counter
	    local_cell++;
	    cell++;
	}
	Ensure (cells_per_proc[proc].size() == capacity);
    }

    // now calculate the boundary cells on each processor
    bound_cells.resize(C4::nodes());
    
    // we sweep through processors and the cells on processor, for each cell
    // we check its neighbors, if a neighbor does not sit on the processor we 
    // define it as a boundary cell in the bound_cell container.  NOTE: the
    // cell indices in use here are global cell indices.  We do not use local 
    // (on-processor) cell indices at this point.  FYI, the local cell index
    // on processor for each global cell is simple the index+1 of the
    // cells_per_proc container.

    // we need a const_iterator to a vector<int> to store the results of
    // std::find 
    vector<int>::const_iterator itr;

    // markers for our sweep through the mesh
    int this_cell;
    int next_cell;

    // sweep through each processor
    for (int proc = 0; proc < C4::nodes(); proc++)
    {
	// sweep through the cells on processor
	for (int i = 0; i < cells_per_proc[proc].size(); i++)
	{
	    // get the cell index of a local (on-processor) cell
	    this_cell = cells_per_proc[proc][i];

	    // determine the neighbors of that cell-->this requires a
	    // global mesh
	    vector<int> neighbors = mesh->get_neighbors(this_cell);

	    // sweep through neighboring cells and see if they: (a) exist on
	    // another processor (b) are not off the problem boundary (have
	    // an index of 0)
	    for (int face = 0; face < neighbors.size(); face++)
	    {
		// neighboring cell index
		next_cell = neighbors[face];

		// check to see if the cell is on processor
		itr = find(cells_per_proc[proc].begin(),
			   cells_per_proc[proc].end(), next_cell);

		// if the cell is off processor and is not a boundary then
		// add it to the bound_cells vector
		if (itr == cells_per_proc[proc].end() && next_cell)
		    bound_cells[proc].push_back(next_cell);
	    }
	}
    }
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Topology_Builder.t.hh
//---------------------------------------------------------------------------//
