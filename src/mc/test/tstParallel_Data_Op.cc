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
using rtt_mc::Parallel_Data_Operator;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc_test::MC_Interface;
using dsxx::SP;

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
    SP<MC_Interface> interface(new MC_Interface());
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
    SP<MC_Interface> interface(new MC_Interface());
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
    SP<MC_Interface> interface(new MC_Interface());
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
