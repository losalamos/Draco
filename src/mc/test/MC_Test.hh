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

bool topology_replication_test(dsxx::SP<rtt_mc::OS_Mesh> mesh, 
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

    // if we haven't returned than we pass
    return p;
}

} // end namespace rtt_mc_test

#endif                          // __mc_test_MC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/MC_Test.hh
//---------------------------------------------------------------------------//
