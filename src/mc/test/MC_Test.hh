//----------------------------------*-C++-*----------------------------------//
// MC_Test.hh
// Thomas M. Evans
// Wed Apr 21 19:33:36 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Some free functions that we will use throughout the MC test suite
//---------------------------------------------------------------------------//

#ifndef __mc_test_MC_Test_hh__
#define __mc_test_MC_Test_hh__

#include "../OS_Mesh.hh"
#include "../Topology.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include <vector>
#include <iostream>
#include <string>

namespace rtt_mc_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
	      << std::endl;
    return false;
}

#define FAILURE fail(__LINE__, __FILE__);

//===========================================================================//
// 2D Objects Interface
//===========================================================================//
// make an interface for a 6 cell mesh

class MC_Interface
{
    typedef std::vector<std::vector<double> > vec_double;
    typedef std::vector<std::string> vec_string;

  private:
    std::string coord;
    vec_double fine_edge;
    vec_string bnd;

  public:
    // constructor
    inline MC_Interface();
    
    // public copy functions for Mesh
    std::string get_coordinates() const {return coord;}
    vec_double get_fine_edge() const {return fine_edge;} 
    vec_string get_boundaries() const {return bnd;}
};

// constructor
MC_Interface::MC_Interface()
    : coord("xy"), fine_edge(2), bnd(4)
{
    // calculate the fine edges
    fine_edge[0].resize(4);
    fine_edge[1].resize(3);

    fine_edge[0][0] = -1;
    fine_edge[1][0] = -1;

    for (int i = 1; i < fine_edge[0].size(); i++)
	fine_edge[0][i] = fine_edge[0][i-1] + 1;
    for (int i = 1; i < fine_edge[1].size(); i++)
	fine_edge[1][i] = fine_edge[1][i-1] + 2;

    // calculate the boundaries
    for (int i = 1; i < 4; i++)
	bnd[i] = "vacuum";
    bnd[0] = "reflect";
}

//===========================================================================//
// Topology class tests
//===========================================================================//
// Full replication topology test --> uses base class reference to derived
// class type

bool topology_replication_test(rtt_dsxx::SP<rtt_mc::OS_Mesh> mesh, 
			       const rtt_mc::Topology &top)
{
    using std::vector;

    // passing condition
    bool p = true;

    // test num_cells
    if (top.num_cells() != mesh->num_cells()) p = FAILURE;
    if (top.num_cells(C4::node()) != mesh->num_cells()) p = FAILURE;
    
    // test num procs
    for (int cell = 1; cell <= top.num_cells(); cell++)
	if (top.num_procs(cell) != C4::nodes()) p = FAILURE;

    // verify the parallel scheme
    if (top.get_parallel_scheme() != "replication") p = FAILURE;

    // test global_cell functions
    for (int lc = 1; lc <= top.num_cells(C4::node()); lc++)
	if (top.global_cell(lc) != lc) p = FAILURE;

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int lc = 1; lc <= top.num_cells(proc); lc++)
	    if (top.global_cell(lc, proc) != lc) p = FAILURE;

    // test local cell functions
    for (int gc = 1; gc <= top.num_cells(); gc++)
	if (top.local_cell(gc) != gc) p = FAILURE;

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int gc = 1; gc <= top.num_cells(); gc++)
	    if (top.local_cell(gc, proc) != gc) p = FAILURE;

    // test boundary cell functions
    for (int proc = 0; proc < C4::nodes(); proc++)
    {
	if (top.get_boundary_cells(proc) != 0) p = FAILURE; 
	for (int gc = 1; gc <= top.num_cells(); gc++)
	    if (top.global_to_boundary(gc, proc) != 0) p = FAILURE;
    }

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int bc = 1; bc <= top.get_boundary_cells(proc); bc++)
	{ 
	    // shouldn't get here because there are no boundary
	    // cells
	    p = FAILURE; 
	}

    // test get_cells function
    {
	vector<int> cell_list(mesh->num_cells());
	for (int i = 1; i <= cell_list.size(); i++)
	    cell_list[i-1] = i;

	// compare
	for (int np = 0; np < C4::nodes(); np++)
	{
	    vector<int> got_cells = top.get_cells(np);
	    if (got_cells != cell_list) p = FAILURE;
	}
    }

    // test get_procs function
    {
	vector<int> proc_list(C4::nodes());
	for (int i = 0; i < C4::nodes(); i++)
	    proc_list[i] = i;

	// compare
	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    vector<int> got_procs = top.get_procs(cell);
	    if (got_procs != proc_list) p = FAILURE;
	}
    }

    // return passing condition
    return p;
}

//---------------------------------------------------------------------------//
// Full DD topology test --> uses base class reference to derived
// class type

bool topology_DD_test(rtt_dsxx::SP<rtt_mc::OS_Mesh> mesh,
		      const rtt_mc::Topology &top)
{
    using std::vector;
    using std::cout;
    using std::endl;

    // passing condition
    bool p = true;

    // verify the parallel scheme
    if (top.get_parallel_scheme() != "DD") p = FAILURE;

    // check the number of cells, this is a 2 processor/3 cell mesh
    // decomposition
    if (top.num_cells() != mesh->num_cells()) p = FAILURE;
    if (top.num_cells(C4::node()) != 3) p = FAILURE;

    // check the number of procs--each cell should be on one processor
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	if (top.num_procs(cell) != 1) p = FAILURE;

    // check global cell indices
    if (C4::node() == 0)
    {
	if (top.global_cell(1) != 1) p = FAILURE;
	if (top.global_cell(2) != 2) p = FAILURE;
	if (top.global_cell(3) != 3) p = FAILURE;
    }

    if (C4::node() == 1)
    {
	if (top.global_cell(1) != 4) p = FAILURE;
	if (top.global_cell(2) != 5) p = FAILURE;
	if (top.global_cell(3) != 6) p = FAILURE;
    }

    for (int np = 0; np < C4::nodes(); np++)
    {
	int offset = np * 3;
	for (int lc = 1; lc <= 3; lc++)
	    if (top.global_cell(lc, np) != lc + offset) p = FAILURE;
    }
	

    // check local cell indices
    if (C4::node() == 0)
    {
	if (top.local_cell(1) != 1) p = FAILURE;
	if (top.local_cell(2) != 2) p = FAILURE;
	if (top.local_cell(3) != 3) p = FAILURE;
	if (top.local_cell(4) != 0) p = FAILURE;
	if (top.local_cell(5) != 0) p = FAILURE;
	if (top.local_cell(6) != 0) p = FAILURE;
    }

    if (C4::node() == 1)
    {
	if (top.local_cell(1) != 0) p = FAILURE;
	if (top.local_cell(2) != 0) p = FAILURE;
	if (top.local_cell(3) != 0) p = FAILURE;
	if (top.local_cell(4) != 1) p = FAILURE;
	if (top.local_cell(5) != 2) p = FAILURE;
	if (top.local_cell(6) != 3) p = FAILURE;
    }

    {
	if (top.local_cell(1, 0) != 1) p = FAILURE;
	if (top.local_cell(2, 0) != 2) p = FAILURE;
	if (top.local_cell(3, 0) != 3) p = FAILURE;
	if (top.local_cell(4, 0) != 0) p = FAILURE;
	if (top.local_cell(5, 0) != 0) p = FAILURE;
	if (top.local_cell(6, 0) != 0) p = FAILURE;
	if (top.local_cell(1, 1) != 0) p = FAILURE;
	if (top.local_cell(2, 1) != 0) p = FAILURE;
	if (top.local_cell(3, 1) != 0) p = FAILURE;
	if (top.local_cell(4, 1) != 1) p = FAILURE;
	if (top.local_cell(5, 1) != 2) p = FAILURE;
	if (top.local_cell(6, 1) != 3) p = FAILURE;
    }

    // check boundary cells
    if (top.get_boundary_cells(0) != 3) p = FAILURE;
    if (top.get_boundary_cells(1) != 3) p = FAILURE;
    {
	if (top.global_to_boundary(1, 0) != 0) p = FAILURE;
	if (top.global_to_boundary(2, 0) != 0) p = FAILURE;
	if (top.global_to_boundary(3, 0) != 0) p = FAILURE;
	if (top.global_to_boundary(4, 0) != 1) p = FAILURE;
	if (top.global_to_boundary(5, 0) != 2) p = FAILURE;
	if (top.global_to_boundary(6, 0) != 3) p = FAILURE;
	if (top.global_to_boundary(1, 1) != 1) p = FAILURE;
	if (top.global_to_boundary(2, 1) != 2) p = FAILURE;
	if (top.global_to_boundary(3, 1) != 3) p = FAILURE;
	if (top.global_to_boundary(4, 1) != 0) p = FAILURE;
	if (top.global_to_boundary(5, 1) != 0) p = FAILURE;
	if (top.global_to_boundary(6, 1) != 0) p = FAILURE;
    }
    
    for (int np = 0; np < C4::nodes(); np++)
    {
	int offset;
	if (np == 0) offset = 3;
	if (np == 1) offset = 0;
	for (int bc = 1; bc <= top.get_boundary_cells(np); bc++)
	    if (top.boundary_to_global(bc, np) != bc + offset) p = FAILURE;
    }

    if (top.boundary_to_global(1, 0) != 4) p = FAILURE;
    if (top.boundary_to_global(2, 0) != 5) p = FAILURE;
    if (top.boundary_to_global(3, 0) != 6) p = FAILURE;
    if (top.boundary_to_global(1, 1) != 1) p = FAILURE;
    if (top.boundary_to_global(2, 1) != 2) p = FAILURE;
    if (top.boundary_to_global(3, 1) != 3) p = FAILURE;

    // check cell lists on processor
    vector<int> cell_list(3);
    if (C4::node() == 0) 
    {
	cell_list[0] = 1;
	cell_list[1] = 2;
	cell_list[2] = 3;
	if (top.get_cells(0) != cell_list) p = FAILURE;
    }
    
    if (C4::node() == 1) 
    {
	cell_list[0] = 4;
	cell_list[1] = 5;
	cell_list[2] = 6;
	if (top.get_cells(1) != cell_list) p = FAILURE;
    }

    // check proc lists for global cells
    vector<int> proc_list(1);
    proc_list[0] = 0;
    if (top.get_procs(1) != proc_list) p = FAILURE;
    if (top.get_procs(2) != proc_list) p = FAILURE;
    if (top.get_procs(3) != proc_list) p = FAILURE;
    proc_list[0] = 1;
    if (top.get_procs(4) != proc_list) p = FAILURE;
    if (top.get_procs(5) != proc_list) p = FAILURE;
    if (top.get_procs(6) != proc_list) p = FAILURE;

    // return passing condition
    return p;
}

} // end namespace rtt_mc_test

#endif                          // __mc_test_MC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/MC_Test.hh
//---------------------------------------------------------------------------//
