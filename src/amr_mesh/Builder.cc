//----------------------------------*-C++-*----------------------------------//
// Builder.cc
// B.T. Adams (bta@lanl.gov)
// 18 May 99
//
//
//                                 >>>>>>
//                                   >>>>>>
//                                     >>>>>>
//                                       >>>>>>
//                                     >>>>>>
//                                   >>>>>>
//                                 >>>>>>     +1
//                                    |     -----
//                                    |     |
//                                    |     |
//                                    |     |
//                             _____________|
//                             |      |
//                             |      |
//                          -1 |      |
//                         ____|      |
//                                    |
//                                    |
//                                    |
//                                    |
//                                    |
//                                   / \
//                                  /   \
//                                 /     \
//                                /       \
//                               /         \
//                              /           \
//                             /             \
//                            /               \
//                           /                 \
//                          /                   \
//                      ----                     ----
//---------------------------------------------------------------------------//
// @> Builder class implementation file (developed from OS_Builder.cc)
/*! 
 * \file   amr_mesh/Builder.cc
 * \author B.T. Adams
 * \date   Tue May 18 10:33:26 1999
 * \brief  Implementation file for CAR_CU_Builder class library.
 */
//---------------------------------------------------------------------------//

#include "Builder.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <algorithm>
#include <iterator>
#include <functional>

namespace rtt_amr 
{
using std::cout;
using std::endl;
using std::set;
using std::merge;
using std::lower_bound;
using std::insert_iterator;
using std::sort;
using std::multimap;
using std::make_pair;
using std::max;
using std::transform;
using std::plus;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public Mesh build member functions
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// The mesh can be constructed using four alternative architechtures: 
//     unstructured (US)   - obtained with a single call to assign_Generations
//                           and no calls to the member functions assign_Meshes
//                           and seek_Adoption.
//     linked child (LC)   - this non-hiearchial mesh structure is obtained by
//                           calling assign_Generations followed by a single 
//                           call to assign_Meshes. This forms a series of
//                           generation-specific Cartesian meshes that are 
//                           bounded by either the problem physical geometry
//                           or other meshes composed of cells at a different
//                           generation level.
//     nested complex (NC) - this hiearchial mesh structure is obtained by 
//                           sequentially calling assign_Generations, 
//                           assign_Meshes, and seek_Adoption. Cells that have
//                           been refined contain a "complex" composed of the
//                           8 (3D) or 4 (2D) child cells that were grouped to
//                           form the parent cell.
//     hiearchial linked child (HLC) - This heirachial extension of the linked
//                           mesh architechture is obtained by first forming a 
//                           linked child mesh as described above, and then 
//                           sequently calling seek_Adoption and assign_Meshes.
//                           The second call to assign_Meshes generates sets
//                           of generation-specific Cartesian mesh composed of
//                           both child and parent cells. 
// Currently this is not a user option, but must be manually configured prior
// to compilation. It may be deemed desirable to make this a user selectable
// feature in the future (but my life is complicated enough as is right now).
// While some capability for the hierarchial mesh structures has already been
// developed, this feature can not yet be fully implemented herein. 
//---------------------------------------------------------------------------//

SP<CAR_CU_Mesh> CAR_CU_Builder::build_Mesh(const SP<RTT_Format> & rttMesh)
{
  // declare smart pointers
    SP<Layout> layout;
    SP<CAR_CU_Mesh> return_mesh;
    
    // Set up some useful RTT Format variables
    // total number of nodes.
    int nnodes = rttMesh->get_dims_nnodes();
    // total number of cells.
    int ncells = rttMesh->get_dims_ncells();
    // problem dimension (2D or 3D).
    int ndim = rttMesh->get_dims_ndim();
    // maximum number of nodes per cell.
    int nnodes_max = rttMesh->get_dims_nnodes_max();
    // maximum number of sides per cell.
    int nsides_max = rttMesh->get_dims_nsides_max();

    // initialization variables for Mesh
    CAR_CU_Mesh::ccsf_i generation(ncells);
    // set of cells in each generation
    vector<set<int> > gen_cells;
    // assign each cell to a generation and receive a vector of sets with the
    // cells per generation.
    gen_cells = assign_Generations(generation, rttMesh);

    // multimap containing sets composed of the cells that are used to create
    // a new parent cell. The parent cell number is taken to be the same as
    // that of the first cell used in it's creation, and this cell number is 
    // also the multimap key. The proud_parents multimap is not used in the
    // first call to assign_Meshes, but it must be declared prior to the call.
    multimap<int, set<int> > proud_parents;

    // assign each cell to a generation-specific mesh in child_cells. There 
    // can be more than one mesh at a given generation level if the cells are
    // not connected (e.g., regions of cells that are separated by cells that
    // are at either a higher or lower generation level). The proud_parents 
    // multimap is null for this first call
    vector<set<int> > child_cells = assign_Meshes(gen_cells, proud_parents, 
						  rttMesh);

    // Determine the recursive nature of the AMR mesh. Create parent cells for
    // the refined regions, and assign the child cells to the parents in the
    // multimap proud_parents. Note that none of the pre-existing cells in the
    // RTT Format file can have children because a hiearchy does not exist in
    // the data. A "gap" will always exist between subsequent cell numbers at
    // a given generation level where refined cells are connected, so the new 
    // cell numbers can be taken to be the same as that of the first refined 
    // cell that is being "grouped" to create the parent cell. Add these new 
    // cell numbers to the respective sets in gen_cells.
    proud_parents = seek_Adoption(gen_cells, rttMesh);

    // assign both the child and parent cells to a generation-specific mesh 
    // in parent_cells. The parents present in these sets can join separate 
    // child regions and thus allow consideration of multigrid acceleration 
    // techniques.
    vector<set<int> > parent_cells = assign_Meshes(gen_cells, proud_parents, 
						   rttMesh);

    // renormalize the nodes specific to this mesh. This has the advantage
    // that, at the present time, ICEM writes several nodes that occur outside
    // of the problem's physical geometry. These nodes are not used in the cell
    // definitions, but they do appear in the RTT_Format node data. These 
    // unused nodes are eliminated due to the renormalization, which could 
    // save considerable computational effort in algorithms that search for a 
    // node number. The index of vector node_map is the renormalized node 
    // number and the original node number is stored.
    vector<int> node_map = map_Nodes(rttMesh);
    nnodes = node_map.size();

    // input the existing node coordinate values from the RTT Format data.
    node_coords.resize(ndim);
    for (int d = 0 ; d < ndim; d++)
    {
        // resize node_coords to accomodate the existing cell-corner nodes (the
        // vector could also be resized here to accomodate the face-centered 
        // nodes that will be subsequently created if FC_Nodes is called, but
        // some programs may not need this data so we don't do it here).
        node_coords[d].resize(nnodes);
        for (int node = 0; node < nnodes; node++)
            node_coords[d][node] = rttMesh->get_nodes_coords(node_map[node],d);
    }

    // input the existing cell_nodes data (nodes that make up a cell) from the 
    // RTT Format data.
    cell_nodes.resize(ncells);
    for (int cell = 0; cell < ncells; cell++)
    {
        // size the cell pair vector to accomodate the cell-corner nodes that
        // are defined in the RTT_Format file (will need to be resized again
        // if face-centered nodes are subsequently created, but better to 
        // force a second resize than accomodate a lot of unneccessary storage
        // if these nodes are not needed).
        cell_nodes[cell].resize(nnodes_max);
        for (int node = 0; node < nnodes_max; node++)
	{
	    int * new_node = lower_bound(node_map.begin(), node_map.end(),
					 rttMesh->get_cells_nodes(cell, node));
	    // cell_nodes numbers are referenced relative to a start at one,
	    // while the RTT Format data is referenced relative to a start at
	    // zero;
	    cell_nodes[cell][node] = 1 + (new_node - node_map.begin());
	}
    }
    // done with the node_map vector;
    node_map.resize(0);

    // generate nodes that are centered on the cell faces, add them to the 
    // cell_nodes vector, calculate their coordinates, and add these to the 
    // node_coords vector.
    FC_Nodes(nnodes, rttMesh);

    // build mesh-independent objects
    layout = build_Layout(rttMesh);

    // create mesh
    SP<CAR_CU_Mesh> mesh_return(new CAR_CU_Mesh(* layout, node_coords, 
       cell_nodes, generation));

    // clean up
    gen_cells.resize(0);
    proud_parents.clear();
    child_cells.resize(0);
    parent_cells.resize(0);

  // return mesh to builder
    return mesh_return;
}

//---------------------------------------------------------------------------//
// Layout build member functions
//---------------------------------------------------------------------------//

SP<Layout> CAR_CU_Builder::build_Layout(const SP<RTT_Format> & rttMesh)
{
    // Set up some useful RTT Format variables
    // total number of cells.
    int ncells = rttMesh->get_dims_ncells();
    // maximum number of sides per cell.
    int nsides_max = rttMesh->get_dims_nsides_max();
    // side flags number for the boundary conditions
    int bndf = rttMesh->get_side_flags_boundary_flag_number();

    // set size of new Layout for the required number of cells
    SP<Layout> layout(new Layout(ncells));
    Layout & local_layout = * layout;

    // set the mesh layout from the RTT_Format data.
    for (int cell = 1; cell <= ncells; cell++)
    {
        // set size of new Layout for the required number of faces
	local_layout.set_size(cell, nsides_max);
        for (int side = 1; side <= nsides_max; side++)
	{
	    // determine the number of cells adjacent to this face and resize
	    // the layout accordingly.
	    int num_adj = rttMesh->get_adjCell_size(cell - 1, side - 1);
	    local_layout.set_adj_size(cell, side, num_adj);

	    // going from a coarse region to a refined region if there is more
	    // than one cell across from this face.
	    if (num_adj > 1)
	    {
	        for (int adjCell = 1; adjCell <= num_adj; adjCell++)
		    local_layout(cell, side, adjCell) = 1 + 
		        rttMesh->get_adjCell(cell-1, side-1, adjCell-1);
	    }
	    else
	    {
	        int adjCell = rttMesh->get_adjCell(cell - 1,side - 1);
	        // a positive number for the adjacent cell at this face 
	        // indicates either a cell of the same generation level or 
	        // going from a refined region to a coarse region.
	        if (adjCell >= 0)
		{
	            local_layout(cell, side, 1) = adjCell + 1;
		}
	        // a negative number indicates the flag number for a cell face
		// on the problem outer boundary. 
	        else
	        {
		    int flag  = - adjCell;
	            // have to subtract one from the cell flag number to get 
		    // the corresponding index.
		    --flag;
	            string name = rttMesh->get_side_flags_flag_name(bndf,flag);
	            // vacuum boundary
	            if ((name[0] == 'v' || name[0] == 'V') && 
	                 name.find_first_not_of("vacumVACUM") == string::npos)
		        local_layout(cell, side, 1) = 0;
	            // reflection boundary
	            else if ((name[0] == 'r' || name[0] == 'R') && 
	                      name.find_first_not_of("reflctionREFLCTION")
			      == string::npos)
		        local_layout(cell, side, 1) = cell;
		    else
		        Insist(0,"Illegal boundary condition found!");
	        }
	    }
	}
    }
    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;
    if (debugging)
    {
        int num_cells = local_layout.get_num_cells();
	for (int cell = 1; cell <= ncells; cell++)
	{
	    int num_faces = local_layout.get_num_cell_faces(cell);
	    for (int side = 1; side <= num_faces; side++)
	    {
	        cout << "cell " << cell << " side " << side << " adjacent ";
		int num_adj = local_layout.get_num_adj_cells(cell, side);
		for (int n = 1; n <= num_adj; n++)
		    cout << " " << local_layout(cell, side, n);
		cout <<  endl;
	    }
	}	
    }
    // return built Layout
    return layout;
}

//---------------------------------------------------------------------------//
// Assign each cell to a generation - modifies the gen_cells ccsf_i directly.
//
// This algorithm originally required two linear traverses through all of the 
// cells, due to the inherently serial nature of initially determining the 
// number of generations present in the mesh based upon the refinement levels
// in previously assigned adjacent cells at the negative x boundary. This was
// 2ds, so the boundary cells were retained in the RTT_Format class for use 
// herein. This reduced the initial traverse of the cells to only include those
// cells that had already been identified to be "boundary" cells in the mesh
// connectivity algorithm, but this implementation is still somewhat inherently
// serial. If this is still 2ds, a parallel scheme can be implemented by 
// creating a routine to first calculate a vector consisting of the relative
// generation levels on the negative x face (currently an integer referred to
// as min_gen) and then calculating minimum and maximum generation levels
// that occur within each "row and column" of cells in the positive x 
// direction in parallel. The second loop to actually assign renormalized 
// generation levels to each cell can then also be executed in parallel.
// Toodles. 
//---------------------------------------------------------------------------//

vector<set<int> > CAR_CU_Builder::assign_Generations(CAR_CU_Mesh::ccsf_i & 
				  cell_gens, const SP<RTT_Format> & rttMesh)
{
    // Set up some useful RTT Format variables
    // problem dimension (2D or 3D).
    int ndim = rttMesh->get_dims_ndim();
    // total number of cells.
    int ncells = rttMesh->get_dims_ncells();
    // maximum number of sides per cell.
    int nsides_max = rttMesh->get_dims_nsides_max();
    // set of boundary cells (either between generation levels or on the 
    // problem physical boundary) in both the negative (face index = ndim - 1)
    // and positive (face index = ndim) x directions.
    set<int> bndryCells;
    set<int> lox = rttMesh->get_bndryCells(ndim - 1);
    set<int> hix = rttMesh->get_bndryCells(ndim);
    insert_iterator<set<int> >  bciiter(bndryCells,bndryCells.begin());
    merge(lox.begin(), lox.end(), hix.begin(), hix.end(), bciiter);
    set<int>::iterator bciter = bndryCells.begin();
    // clear some memory
    lox.clear();
    hix.clear();

    // determine the number of generations present in the mesh by looking for 
    // -x and +x faces that ajoin four (3D) or two (2D) cells, and assign a 
    // preliminary generation level to each cell that has such a boundary.
    int min_gen = 0;      // minimum generation level found in the mesh.
    int max_gen = 0;      // maximum generation level found in the mesh.
    int generation = 0;   // current generation level

    while (bciter != bndryCells.end())
    {
        int cell_number = * bciter;
        // check that this cell is not on the the negative x boundary (as 
        // indicated by the fact that the adjacent cell number on the 
	// ndim - 1 face < 0), or if it is either the first or last cell.
        if (rttMesh->get_adjCell(cell_number,ndim - 1) >= 0 || 
	    cell_number == 0 || cell_number == ncells -1)
	{
            for (int face = ndim - 1; face <= ndim; face++)
	    {
	        // check if the adjacent face connects to multiple cells
	        int adjCell = rttMesh->get_adjCell(cell_number,face);
	        if (adjCell >= 0 && 
		     rttMesh->get_adjCell_size(adjCell,nsides_max-1-face) > 1)
	        {
		    // going from a coarse region to a refined region. 
		    // Increment the current generation level and assign to
		    // this cell.
		    if (face == ndim - 1)
		    {
		        generation = 1 + cell_gens[adjCell];
		        cell_gens[cell_number] = generation ;
		        if (generation > max_gen)
		            max_gen = generation;
		    }
		    // going from a refined region to a coarse region. Assign
		    // the current generation level to this cell and then 
		    // decrement the current generation level.
	            else
		    {
		        cell_gens[cell_number] = generation;
		        --generation;
		        if (generation < min_gen)
		            min_gen = generation;
		    }
	        }
		// this cell is at the same generation level as the previous
		// cell.
	        else
	            cell_gens[cell_number] = generation;
	    }
	}
	// At the negative x boundary. Set the current generation number equal
	// to that of an adjacent cell that has already been assigned (either
	// the -y or -z face assuming serial processing).
	else
	{
	    // check the -y face first, then -z.
	    for (int face = ndim - 2; face >= 0; face--)
	    {

	        int adjCell = rttMesh->get_adjCell(cell_number,face);
		// find an adjoining cell face (the faces on the system 
		// boundary return negative values for the adjacent cell).
		if (adjCell >= 0)
		{
		    // Set the current generation level equal to the adjacent
		    // cell generation level
		    generation = cell_gens[adjCell];
		    // decrement the generation level number if this face 
		    // adjoins multiple cells. 
		    if (rttMesh->get_adjCell_size(cell_number,face) > 1)
		    {
		        --generation;
			if (generation < min_gen)
		            min_gen = generation;
		    }
		    // increment the generation level number if the selected
		    // adjacent cell face connects to multiple cells. 
		    if (rttMesh->get_adjCell_size(adjCell,nsides_max-1-face) 
			> 1)
		    {
		        ++generation;
			if (generation > max_gen)
		            max_gen = generation;
		    }
		    // assign the current generation level to this cell and
		    // end the face loop
		    cell_gens[cell_number] = generation;
		    face = -1;
		}
	    }
	    // decrement the current generation number if the cell in the 
	    // positive x direction is a coarser grid.
	    int adjCell = rttMesh->get_adjCell(cell_number,ndim);
	    if (rttMesh->get_adjCell_size(adjCell,ndim-1) > 1)
	    {
	        --generation;
		if (generation < min_gen)
		    min_gen = generation;
	    }
	}
	// increment the iterator used to establish the boundary cell number
	++bciter;
    }

    // renormalize the cell generation levels to start at zero, and create
    // a vector containing sets of all of the cell numbers at a given 
    // generation level.
    vector<set<int> > gen_cells(max_gen - min_gen + 1);
    bciter = bndryCells.begin();
    while (bciter != bndryCells.end())
    {
        int cell_number = * bciter;
	// renormalize the boundary cell generation level and add this cell 
	// to the generation vector
        if (min_gen != 0)
	    cell_gens[cell_number] -= min_gen;
	gen_cells[cell_gens[cell_number]].insert(cell_number);
	// assign this generation level to all of the cells up to the next 
	// boundary cell
	++bciter;
	++cell_number;
	while (bciter != bndryCells.end() && cell_number < * bciter)
	{
	    cell_gens[cell_number] = cell_gens[cell_number - 1];
	    gen_cells[cell_gens[cell_number]].insert(cell_number);
	    ++cell_number;
	}
    }
    // make sure every cell was assigned to a generation and echo the number
    // of cells per generation to the screen (the unstructured mesh)
    cout << endl;
    cout << "--- Resolving unstructured mesh ---" << endl;
    int gen_check = 0;
    for (int gen = 0; gen < gen_cells.size(); gen ++)
    {
        gen_check += gen_cells[gen].size();
	cout << "cells for mesh at generation " << gen 
	     << " = " << gen_cells[gen].size() << endl;	
    }
    Insist(gen_check == ncells,"Cell generation assignment failure!");

    return gen_cells;
}

//---------------------------------------------------------------------------//
// Assign children to parents - establish a hiearchy for the mesh.
//---------------------------------------------------------------------------//

multimap<int, set<int> > CAR_CU_Builder::seek_Adoption(vector<set< int> > & 
    gen_cells, const SP<RTT_Format> & rttMesh)
{

    // Set up some useful RTT Format variables
    // problem dimension (2D or 3D).
    int ndim = rttMesh->get_dims_ndim();
    // returned cells that were generated.
    multimap<int, set<int> > new_cells;

    // start at the highest generation level (most refined) and work down.
    for (int gen = gen_cells.size() - 1; gen > 0; gen--)
    {
        // local set of cell numbers in a generation that we can destroy.
        set<int> local_gen = gen_cells[gen];
	set<int>::const_iterator citer;
	int parent_check = 0;
	int orphan_count = 0;
	while (!local_gen.empty())
	{
	    citer = local_gen.begin();
	    set<int> grouped_cells;
	    int cell_number = * citer;
	    grouped_cells.insert(cell_number);
	    // look for cell groupings. Recursive loop finds cells in the
	    // positive x,y, and z directions, increments the cell number to
	    // the adjacent cell in the positive z direction and finds the 
	    // connected cells in the positive x and y directions, and then 
	    // increments the cell number to the adjacent cell in the positive
	    // y direction and finds the adjacent cell in the positive x 
	    // direction. Works similiarly for 2D.
	    for (int dir = ndim; dir < 2 * ndim; dir++)
	    {
	        int adjCell;
	        for (int face = ndim; face < 3 * ndim - dir; face++)
		{
	            adjCell = rttMesh->get_adjCell(cell_number,face);
		    // have to skip over cells in the same direction to get
		    // the correct adjacent cell if this is a parent cell.
		    if (new_cells.count(cell_number) > 0)
		    {
			while(adjCell >= 0 && local_gen.count(adjCell) == 0 &&
			      gen_cells[gen - 1].count(adjCell) == 0)
			    adjCell = rttMesh->get_adjCell(adjCell,face);
		    }
		    if (adjCell >= 0 && local_gen.count(adjCell) != 0)
		        grouped_cells.insert(adjCell);
		}
		cell_number = adjCell;
		// exit the loop if the next cell number is actually just a 
		// boundary condition. This is yet another problem caused by
		// the hopeless orphan cells (see comments below)!
		if (cell_number < 0)
		    dir = 2 * ndim;
	    }
	    // have to pick up the adjacent cell in the negative z direction
	    // for 3D. The cell number currently sits on top of this last cell 
	    // unless it is less than zero (which means that we are sitting on
	    // top of a hopeless orphan) or if it is a parent cell.
	    if (ndim == 3 && cell_number > 0)
	    {
	        int adjCell = rttMesh->get_adjCell(cell_number,0);
		// have to skip over a cell in the same direction to get the 
		// right adjacent cell if the adjacent cell has been grouped 
		// into a parent cell (note that this differs from the check
		// in the positive directions).
		while(adjCell >= 0 && local_gen.count(adjCell) == 0 &&
		      gen_cells[gen - 1].count(adjCell) == 0)
		    adjCell = rttMesh->get_adjCell(adjCell,0);

		if (adjCell >= 0 && local_gen.count(adjCell) != 0)
		    grouped_cells.insert(adjCell);
	    }

	    // assign this group of children to the parent generation and 
	    // create the parent cell. Assign the cell number of the new
	    // parent cell to be equal to the first cell in the group
	    cell_number = * grouped_cells.begin();
	    if (grouped_cells.size() == static_cast<int>(pow(2.0,ndim)))
	    {
		gen_cells[gen - 1].insert(cell_number);
		new_cells.insert(make_pair(cell_number,grouped_cells));
		for (set<int>::const_iterator citer = grouped_cells.begin(); 
		     citer != grouped_cells.end(); citer++)
		    local_gen.erase(* citer);
		parent_check += static_cast<int>(pow(2.0,ndim));
	    }
	    // flag a hopeless orphan (see comments below)
	    else
	    {
	        local_gen.erase(cell_number);
		++orphan_count;
	    }
	    grouped_cells.clear();
	}
	// ICEM currently produces some meshes that are not "legal" CAR meshes
	// in that children cells exist that cannot be grouped into parent 
	// cells. These illegitimate cells occur on the outer boundary of the
	// problem. There is no way to adopt such unlovable kids, so we can 
	// only flag the condition (which may challenge multigrid acceleration 
	// schemes since the entire mesh cannot be represent in coarse levels).
	// There is typically only a single orphan cell in any given direction
	// on the boundary. We insist that the number of orphan cells plus the
	// number of parents generated * 2^dim must equal the size of the 
	// generation to trap adoption failures.
	if (orphan_count != 0)
	    cout << "Warning: " << orphan_count << " orphan cells were found "
		 << "in generation " << gen << endl;
	Insist(parent_check + orphan_count == gen_cells[gen].size(),
	       "Illegitimate orphan children found in the mesh!");
    }

    return new_cells;
}

//---------------------------------------------------------------------------//
// Create hierachial meshes. Generates sets of linked child mesh cells if the
// parents multimap is null and sets of hiearchial mesh cells otherwise. 
//---------------------------------------------------------------------------//

vector<set< int> > CAR_CU_Builder::assign_Meshes(const vector<set< int> > & 
    gen_cells, const multimap<int, set<int> > & parents, 
    const SP<RTT_Format> & rttMesh)
{
    // Assign each cell to a mesh. There will be at least one mesh for each 
    // generation level, but there can be more than one mesh per generation
    // if all of the cells at a generation level are not interconnected.
    vector<set< int> > mesh_cells;

    // Set up some useful RTT Format variables
    // problem dimension (2D or 3D).
    int ndim = rttMesh->get_dims_ndim();
    // total number of cells.
    int ncells = rttMesh->get_dims_ncells();
    // maximum number of sides per cell.
    int nsides_max = rttMesh->get_dims_nsides_max();
    // current mesh number (index to mesh_cells)
    int mesh_number = 0;
    // check to ensure that every cell is assigned to a mesh
    int mesh_check = 0;

    // screen echo to distinguish between a child mesh and a parent mesh
    cout << endl;
    if (parents.empty())
      cout << "--- Resolving child meshes ---" << endl;
    else
      cout << "--- Creating parent meshes ---" << endl;

    for (int gen = 0; gen < gen_cells.size(); gen++)
    {
        // local set of cells in this generation that we can destroy
        set<int> local_gen = gen_cells[gen];
        // Assign the first cell at this generation level to the current mesh.
	mesh_cells.resize(mesh_number + 1);
	mesh_cells[mesh_number].insert(* local_gen.begin());
	vector<int> cell_set;
	cell_set.push_back(* local_gen.begin());
        int mc_index = 0;
        local_gen.erase(* local_gen.begin());

	int connected = 0;
	// the rest of this is redundant if this is the last generation and
	// we are including the parent cells.
	if (gen == gen_cells.size() - 1 && !parents.empty())
	{
	    connected = gen_cells[gen].size();
	    mesh_check += connected - 1;
	}
	// find all of the cells at this generation level that are connected 
	// to the current mesh. Keep count of connected cells in connected
	while (connected != gen_cells[gen].size())
	{
	    // if this is not the first time through the loop reset the cell 
	    // iterator to the first cell in this generation that was not 
	    // included in a previous mesh, increment the mesh number, resize 
	    // the mesh_cells vector, and insert this cell in the new mesh. 
	    if (connected != 0)
	    {
		mesh_check += mesh_cells[mesh_number].size();
		++mesh_number;
		mesh_cells.resize(mesh_number + 1);
		mesh_cells[mesh_number].insert(* local_gen.begin());
		cell_set.push_back(* local_gen.begin());
		local_gen.erase(* local_gen.begin());
		mc_index = 0;
	    }
	    while (mc_index < cell_set.size())
	    {
		int cell_number = cell_set[mc_index];
	        // Look for adjacent connecting faces in both the negative
		// and positive x,y, (and possibly z) directions that are
		// at the same generation level. Due to the latter, use the
		// default of the get_adjCell function of only considering
		// the first index of the last element of adjCell. The need
		// to check for negative cell faces is imposed by the need
		// for the mesh cells to be able to "turn around negative 
		// corners" (e.g., top_hat). This is also the reason that
		// the cell_set vector (which is not sorted) is used.
		for (int face = 0; face < nsides_max; face++)
		{
		    // check for connected cells in all directions at this 
		    // generation level.
		    int adjCell = rttMesh->get_adjCell(cell_number,face);
		    // these adjustments only apply if we are trying to set 
		    // up a parent mesh
		    if (!parents.empty())
		    {
		        // if this is a parent cell and we are looking in the
		        // positive directions then we have to skip over a cell
		        // in the same direction to get the right adjacent cell
		        if (adjCell >= 0 && parents.count(cell_number) != 0 && 
			        gen != gen_cells.size() - 1 && 
			        gen_cells[gen + 1].count(cell_number) != 0)
		            adjCell = rttMesh->get_adjCell(adjCell,face);

		        // May have to skip over cells if the adjacent cell is
		        // really a parent cell, as indicated by the fact that 
		        // the indicated cell appears in a lower generation 
			// level and not in the current generation.
		        while (adjCell >= 0 && 
			       gen_cells[gen].count(adjCell) == 0 
			       && (gen == 0 ||  
			       gen_cells[gen - 1].count(adjCell) == 0))
			    adjCell = rttMesh->get_adjCell(adjCell,face);
		    }
		    // make sure that this cell has not already been added and
		    // is part of the current generation
		    if (adjCell >= 0 && 
			mesh_cells[mesh_number].count(adjCell) == 0 &&
			gen_cells[gen].count(adjCell) != 0)
		    {
		        mesh_cells[mesh_number].insert(adjCell);
			cell_set.push_back(adjCell);
		        local_gen.erase(adjCell);
		    }
		}
	        ++mc_index;
	    }
	    // found all the connected cells in this mesh. Increment connected
	    // to test if all the cells in this generation have been assigned.
	    // Echo the mesh structure to the screen.
	    connected += mesh_cells[mesh_number].size();
	    cell_set.resize(0);
	    cout << "cells for mesh " << mesh_number << " at generation " 
		 << gen << " = " << mesh_cells[mesh_number].size() << endl;
	}
	// done with this generation level. Increment mesh check. The mesh 
	// number increments with the generation level by definition.
	mesh_check += mesh_cells[mesh_number].size();
	++mesh_number;
    }
    Insist(mesh_check == ncells + parents.size(), 
	   "Error in cell mesh assignment!");

    return mesh_cells;
}

//---------------------------------------------------------------------------//
// Renormalize the node ordering to include only those nodes actually used in
// this mesh. This eliminates superfluous nodes in the ICEM RTT Format file 
// and is also "required" (okay, desireable) for a regressive mesh type.
//---------------------------------------------------------------------------//

vector<int> CAR_CU_Builder::map_Nodes(const SP<RTT_Format> & rttMesh)
{
    // Renormalize node ordering to include only those nodes actually used in
    // this mesh.

    // Set up some useful RTT Format variables
    // problem dimension (2D or 3D).
    int ndim = rttMesh->get_dims_ndim();
    // total number of cells.
    int ncells = rttMesh->get_dims_ncells();
    // maximum number of nodes per cell.
    int nnodes_max = rttMesh->get_dims_nnodes_max();
    // original node numbers
    set<int> org_node;
    // returned node map
    vector<int> node_map;

    for (int cell = 0; cell < ncells; cell++)
    {
	 for (int node = 0; node < nnodes_max; node++)
	     org_node.insert(rttMesh->get_cells_nodes(cell, node));
    }

    for (set<int>::const_iterator niter = org_node.begin();
	 niter != org_node.end(); niter++)
        node_map.push_back(* niter);

    return node_map;
}

//---------------------------------------------------------------------------//
// generate nodes that are centered on the cell faces, add them to the 
// cell_nodes vector, calculate their coordinates, and add these to the 
// node_coords vector.
//---------------------------------------------------------------------------//

/*!
 * \brief Generate nodes that are centered on the cell faces, adds them to the
 *        cell_nodes vector, calculates their coordinates, and adds the 
 *        coordinates to the node_coords vector. The numbering of the 
 *        face-centererd nodes begins after all of the cell corner nodes.
 * \param nnodes Total number of nodes (corner nodes only on input and both
 *               corner and face-centered nodes on output).
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 */
void CAR_CU_Builder::FC_Nodes(int nnodes, const SP<RTT_Format> & rttMesh)
{
    // Set up some useful RTT Format variables
    // total number of cells.
    int ncells = rttMesh->get_dims_ncells();
    // problem dimension (2D or 3D).
    int ndim = rttMesh->get_dims_ndim();
    // maximum number of nodes per cell.
    int nnodes_max = rttMesh->get_dims_nnodes_max();
    // maximum number of sides per cell.
    int nsides_max = rttMesh->get_dims_nsides_max();

    // determine the number of faces on the boundary of the problem in the 
    // positive x, y, and z directions. This information will be used to 
    // resize the node_coords vector
    int bndryCellFaces = 0;
    for (int face = ndim; face < 2 * ndim; face++)
    {
        set<int> bndryFace = rttMesh->get_bndryCells(face);
	for (set<int>::const_iterator fiter = bndryFace.begin();
	     fiter !=  bndryFace.end(); fiter++)
	{
	    if (rttMesh->get_adjCell(* fiter, face) < 0)
	        ++bndryCellFaces;
	}
    }
    // resize the node_coords vector to accomodate the face-centered nodes
    // (before the value of nnodes is changed).
    for (int d = 0 ; d < ndim; d++)
	 node_coords[d].resize(nnodes + ndim * ncells + bndryCellFaces);

    // We will assume that all of our cells have the same type, as defined
    // in the RTT Format file, so just consider cell 0 in the function call.
    int cellType = rttMesh->get_cells_type(0);
    // create face-centered nodes, calulate their coordinate values and save 
    // in the node_coords vector, and add the face centered values to the end
    // of the cell definitions in the cell_nodes vector.
    for (int cell = 0; cell < ncells; cell++)
    {
        // resize the cell pair vector to accomodate both the cell-corner nodes
        // that are defined in the RTT_Format file and the new face-centered 
        // nodes to be created.
        cell_nodes[cell].resize(nnodes_max + nsides_max);        
	for (int face = 0; face < nsides_max; face++)
	{
	    if (face < nsides_max/2 || rttMesh->get_adjCell(cell, face) < 0)
	    {
	        // get the sorted nodes that define this cell face. The first
	        // and last nodes in the set are on opposite corners of the 
	        // face. Note that the node numbers stored in cell_nodes are
	        // referenced to a start at one while the RTT Format node 
	        // numbers start at zero.
	        const set<int> faceNodes = 
	            rttMesh->get_cell_defs_side(cellType,face);

	        int lo_node = 
		    rttMesh->get_cells_nodes(cell, * faceNodes.begin());
		vector<double> lo_coord = rttMesh->get_nodes_coords(lo_node);
		int hi_node = 
		    rttMesh->get_cells_nodes(cell,* (--faceNodes.end()));
		vector<double> hi_coord = rttMesh->get_nodes_coords(hi_node);

		for (int d = 0 ; d < ndim; d++)
		    node_coords[d][nnodes] = (lo_coord[d] + hi_coord[d])/2.0;
		// numbering for the face-centered nodes to be created will 
		// start after all of the cell-corner nodes and correspond to 
		// faces -z, -y, -x, x, y, z. Note that unique face nodes are 
		// NOT assigned to each cell (saves a lot of memory).
		cell_nodes[cell][nnodes_max + face] = ++nnodes;
		if (face < nsides_max/2)
		    cell_nodes[cell][nnodes_max + face + ndim] = nnodes + ndim;
	    }
	}
    }
    return;
}

} // end namespace rtt_amr

//---------------------------------------------------------------------------//
//                              end of Builder.cc
//---------------------------------------------------------------------------//
