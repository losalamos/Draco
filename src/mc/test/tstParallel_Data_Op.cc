//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstParallel_Data_Op.cc
 * \author Thomas M. Evans
 * \date   Fri Dec 10 13:12:10 1999
 * \brief  Parallel_Data_Operator test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../Rep_Topology.hh"
#include "../General_Topology.hh"
#include "../Parallel_Data_Operator.hh"
#include "../OS_Builder.hh"
#include "../OS_Mesh.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <string>

using namespace std;

using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::Parallel_Data_Operator;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc_test::Parser;
using rtt_dsxx::SP;

typedef OS_Mesh::CCSF<int>           ccsf_int;
typedef OS_Mesh::CCSF<int>::iterator field_itor;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__, __FILE__);

//---------------------------------------------------------------------------//
// TEST Parallel Data Operations on MT Fields
//---------------------------------------------------------------------------//

void test_MT_Fields()
{
    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();
    
    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a Parallel_Data_Operator
    Parallel_Data_Operator pdop(topology);

    // build some mesh fields
    ccsf_int field(mesh);
    field_itor begin = field.begin();
    field_itor end   = field.end();

    for (field_itor itr = begin; itr != end; itr++)
	*itr = 1;

    // now do a global sum of the field
    pdop.global_sum(begin, end);

    // check
    for (field_itor itr = begin; itr != end; itr++)
	if (*itr != C4::nodes()) ITFAILS;
}

template<class BT, class MT>
void test_MT_Fields()
{
    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    BT builder(interface);
    SP<MT> mesh = builder.build_Mesh();
    
    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a Parallel_Data_Operator
    Parallel_Data_Operator pdop(topology);

    // build some fields
    typename MT::CCSF<double> field(mesh);
    typename MT::CCSF<double>::iterator begin = field.begin();
    typename MT::CCSF<double>::iterator end   = field.end();
    typename MT::CCSF<double>::iterator itr;
    
    for (itr = begin; itr != end; itr++)
	*itr = .5;

    // now do a global sum of the field
    pdop.global_sum(begin, end);

    for (itr = begin; itr != end; itr++)
	if (*itr != static_cast<double>(C4::nodes()) * .5) ITFAILS;
}

//---------------------------------------------------------------------------//
// TEST Parallel Data Operations on STL Types and C-style arrays
//---------------------------------------------------------------------------//

template<class FT, class T>
void test_STL_Fields()
{
    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();
    
    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));
    
    // build a Parallel_Data_Operator
    Parallel_Data_Operator pdop(topology);

    // Build and fill a field type
    FT field(mesh->num_cells());
    typename FT::iterator begin = field.begin();
    typename FT::iterator end   = field.end();
    typename FT::iterator itr;

    // fill the field
    for (itr = begin; itr != end; itr++)
	*itr = 1;

    // now do a global sum of the field
    pdop.global_sum(begin, end);

    for (itr = begin; itr != end; itr++)
	if (*itr != C4::nodes()) ITFAILS;

    // now do a C-style "array" field
    T *c_data = new T[mesh->num_cells()];
    
    // fill the c_style array
    for (int i = 0; i < mesh->num_cells(); i++)
	c_data[i] = 1;
    
    // do a global reduction of the c-style array
    pdop.global_sum(c_data, c_data+mesh->num_cells());

    // check
    for (int i = 0; i < mesh->num_cells(); i++)
	if (c_data[i] != C4::nodes()) ITFAILS;
}

//---------------------------------------------------------------------------//
// TESTING GLOBAL EQUIVALENCES CHECK
//---------------------------------------------------------------------------//

template<class T>
void test_equivalence(const T value, const T mod_value)
{
    // set a local_value
    T local_value = value;

    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();
    
    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));
    
    // build a Parallel_Data_Operator
    Parallel_Data_Operator pdop(topology);

    // check if more than 1 node
    if (C4::nodes() > 1)
    {
	// at this point all processors should have the same value
	if (!pdop.check_global_equiv(local_value)) ITFAILS;
	

	// now change the first processor's value
	if (C4::node() == 0)
	    local_value = mod_value;

	if (C4::node() > 0)
	{
	    if (!pdop.check_global_equiv(local_value)) ITFAILS;
	}
	else
	{
	    if (pdop.check_global_equiv(local_value)) ITFAILS;
	}

	// reset all to the same value
	local_value = value;
	if (!pdop.check_global_equiv(local_value)) ITFAILS;
	
	// now change the last processor's value
	if (C4::node() == C4::nodes() - 1)
	    local_value = mod_value;

	if (C4::node() == C4::nodes() - 2)
	{
	    if (pdop.check_global_equiv(local_value)) ITFAILS;
	}
	else
	{
	    if (!pdop.check_global_equiv(local_value)) ITFAILS;
	}
    }
	 
    // reset all to the same value
    local_value = value;
    if (!pdop.check_global_equiv(local_value)) ITFAILS;

    // check if more than 2 nodes
    if (C4::nodes() > 2)
    {
	// now change a middle value
	if (C4::node() == C4::nodes()/2)
	    local_value = mod_value;
	
	if (C4::node() == C4::nodes()/2 - 1)
	{
	    if (pdop.check_global_equiv(local_value)) ITFAILS;
	}
	else if (C4::node() == C4::nodes()/2)
	{
	    if (pdop.check_global_equiv(local_value)) ITFAILS; 
	}
	else
	{
	    if (!pdop.check_global_equiv(local_value)) ITFAILS;
	}
    }
	 
    // reset all to the same value
    local_value = value;
    if (!pdop.check_global_equiv(local_value)) ITFAILS;

    // check if 1 node (serial)
    if (C4::nodes() == 1)
    {
	// do a serial problem test-->this is trivial but we want to check it
	// anyway
	local_value = mod_value;
	if (!pdop.check_global_equiv(local_value)) ITFAILS;
    }
}

// simple timing test
void timing()
{
    // set a value
    int value = 10;

    double begin;
    double end;
	
    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();
    
    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));
    
    // build a Parallel_Data_Operator
    Parallel_Data_Operator pdop(topology);

    if (C4::node() == 0)
	begin = C4::Wtime();

    for (int i = 0; i < 10; i++)
	pdop.check_global_equiv(value);

    if (C4::node() == 0)
	end = C4::Wtime();

    if (C4::node() == 0)
	cout << "Ran for " << end-begin << " seconds on " << C4::nodes()
	     << endl;
}

//---------------------------------------------------------------------------//
// TEST LOCAL-GLOBAL MAPPING
//---------------------------------------------------------------------------//

void test_mapping_replication()
{
    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();
    
    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));
    
    // build a Parallel_Data_Operator
    Parallel_Data_Operator pdop(topology);

    // make a MT::CCSF field
    ccsf_int local(mesh);
    ccsf_int global(mesh);
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	local(cell) = 10;

    // Do a local_global mapping with data-replicated data
    pdop.local_to_global(local, global,
			 Parallel_Data_Operator::Data_Replicated()); 

    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	if (global(cell) != local(cell)) ITFAILS;

    // Do a local_global mapping with data-decomposed data
    pdop.local_to_global(local, global,
			 Parallel_Data_Operator::Data_Decomposed()); 

    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	if (global(cell) != local(cell) * C4::nodes()) ITFAILS;

    // check local data
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	if (local(cell) != 10) ITFAILS;
}

void test_mapping_DD()
{
    // only perform this test on two processors
    if (C4::nodes() != 2) 
	return;

    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();

    // build a topology based on full DD
    Topology::vf_int cpp(C4::nodes());
    Topology::vf_int ppc(mesh->num_cells());
    Topology::vf_int bc(C4::nodes());

    // fill up cells per processor array
    cpp[0].resize(3);
    cpp[1].resize(3);

    cpp[0][0] = 1;
    cpp[0][1] = 2;
    cpp[0][2] = 3;
    cpp[1][0] = 4;
    cpp[1][1] = 5;
    cpp[1][2] = 6;

    // fill up processor per cells array
    for (int i = 0; i < ppc.size(); i++)
	ppc[i].resize(1);

    ppc[0][0] = 0;
    ppc[1][0] = 0;
    ppc[2][0] = 0;
    ppc[3][0] = 1;
    ppc[4][0] = 1;
    ppc[5][0] = 1;

    // fill up boundary cells array
    bc[0].resize(3);
    bc[1].resize(3);
    
    bc[0][0] = 4;
    bc[0][1] = 5;
    bc[0][2] = 6;   
    bc[1][0] = 1;
    bc[1][1] = 2;
    bc[1][2] = 3;

    // build topology
    SP<Topology> topology(new General_Topology(cpp, ppc, bc, "DD")); 
    
    // build a Parallel_Data_Operator
    Parallel_Data_Operator pdop(topology);

    // make fields
    vector<int> local(3, C4::node() + 1);
    vector<int> global_vect(6);
    
    pdop.local_to_global(local, global_vect,
			 Parallel_Data_Operator::Data_Distributed()); 

    // check value of global data
    if (global_vect[0] != 1) ITFAILS;
    if (global_vect[1] != 1) ITFAILS;
    if (global_vect[2] != 1) ITFAILS;
    if (global_vect[3] != 2) ITFAILS;
    if (global_vect[4] != 2) ITFAILS;
    if (global_vect[5] != 2) ITFAILS;

    // make sure local data unchanged
    for (int i = 0; i < local.size(); i++)
	if (local[i] != C4::node() + 1) ITFAILS;

    // now check mapping with different field types (vector and ccsf)
    ccsf_int global(mesh);
    
    pdop.local_to_global(local, global,
			 Parallel_Data_Operator::Data_Distributed()); 

    // check value of global data
    if (global[0] != 1) ITFAILS;
    if (global[1] != 1) ITFAILS;
    if (global[2] != 1) ITFAILS;
    if (global[3] != 2) ITFAILS;
    if (global[4] != 2) ITFAILS;
    if (global[5] != 2) ITFAILS;

    // make sure local data unchanged
    for (int i = 0; i < local.size(); i++)
	if (local[i] != C4::node() + 1) ITFAILS;

}

//---------------------------------------------------------------------------//
// MAIN
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);
    
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (!C4::node())
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    // Fields test
    test_MT_Fields();
    test_MT_Fields<OS_Builder,OS_Mesh>();

    // test with STL and c-style fields
    test_STL_Fields<vector<int>, int>();
    test_STL_Fields<deque<int>, int>();

    // test global equivalences
    test_equivalence(10, 11);           // int
    test_equivalence(10.0001, 11.0001); // double
    test_equivalence(10.0001, 10.0002); // double

    // test local-global cell mapping
    test_mapping_replication();
    test_mapping_DD();

    // status of test
    cout << endl;
    cout <<     "**************************************************" 
	 << endl;
    if (passed) 
    {
        cout << "**** Parallel_Data_Operator Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "**************************************************" 
	 << endl;
    cout << endl;

    cout << "Done testing Parallel_Data_Operator on node " << C4::node() 
	 << endl;

    C4::Finalize();
}
//---------------------------------------------------------------------------//
//                              end of tstParallel_Data_Op.cc
//---------------------------------------------------------------------------//
