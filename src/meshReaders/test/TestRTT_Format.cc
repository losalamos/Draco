//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestRTT_Format.cc
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Format class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestRTT_Format.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::string;

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_meshReaders_test::TestRTT_Format;
    
    return SP<TestApp>(new TestRTT_Format(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_meshReaders_test
{

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::set;
using std::ostream;
using std::fill;
using rtt_meshReaders::RTT_Format;

TestRTT_Format::TestRTT_Format(int argc, char * argv[], ostream & os_in)

    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestRTT_Format" << endl;
}

string TestRTT_Format::version() const
{
    return rtt_meshReaders::release();
}

/*!
 * \brief Tests the RTT_format mesh reader.
 *
 */
string TestRTT_Format::runTest()
{

    // Read and test a simplistic mesh file.
    string filename = "rttdef.mesh";
    RTT_Format mesh(filename);
    pass(" Construct ") << 
	"Read mesh without coreing in or firing an assertion." << endl;

    check_header(mesh);
    check_dims(mesh);
    check_node_flags(mesh);
    check_virtual(mesh);

    // Report results of test.
    if (passed())
    {
	return "All tests passed.";
    }
    return "Some tests failed.";
}

bool TestRTT_Format::check_virtual(const rtt_meshReaders::RTT_Format & mesh)
{
    // Exercize the virtual accessor functions for this mesh.
    bool all_passed = true;

    // Check node coords
    vector<vector<double> > coords = mesh.get_node_coords();
    vector<vector<double> > test_coords(4, vector<double> (3, 0.0));
    for (int i = 1; i < 4; i++)
        test_coords[i][i - 1] = static_cast<double>(i);

    if (test_coords != coords)
    {
        fail(" Node Coordinates ") << "Node coordinates not obtained." << endl;
	all_passed = false;
    }

    // Check coordinate units.
    os() << "Coordinate Units= " << mesh.get_node_coord_units() << endl;

    if (mesh.get_node_coord_units() != "cm")
    {
        fail(" Coordinates Units ") << "Coordinate units not obtained." 
				    << endl;
 	all_passed = false;
    }
     
    // Check element nodes.
    vector<vector<int> > element_nodes = mesh.get_element_nodes();

    const int side_nodes = 3;
    int test_elem_0[side_nodes] = {1, 2, 3};
    int test_elem_1[side_nodes] = {0, 3, 2};
    int test_elem_2[side_nodes] = {0, 1, 3};
    int test_elem_3[side_nodes] = {0, 2, 1};

    const int cell_nodes = 4;
    int test_elem_4[cell_nodes] = {0, 1, 2, 3};

    if (element_nodes[0] != vector<int>(test_elem_0,test_elem_0+side_nodes) ||
	element_nodes[1] != vector<int>(test_elem_1,test_elem_1+side_nodes) ||
	element_nodes[2] != vector<int>(test_elem_2,test_elem_2+side_nodes) ||
	element_nodes[3] != vector<int>(test_elem_3,test_elem_3+side_nodes) ||
	element_nodes[4] != vector<int>(test_elem_4,test_elem_4+cell_nodes))
    {
        fail(" Element Nodes ") << "Element nodes not obtained." << endl;
 	all_passed = false;
    }

    // Check Element Types.
    vector<rtt_meshReaders::Element_Definition::Element_Type>
	element_types = mesh.get_element_types();
    int test_element_types[8] = {0, 1, 3, 5, 10, 8, 12, 15};
    bool got_element_types = true;

    for (int i = 0; i < 8; i++)
        if (element_types[i] != test_element_types[i])
	    got_element_types = false;

    if (!got_element_types)
    {
	fail(" Element Types ") << "Element Types not obtained." << endl;
 	all_passed = false;
    }

    // Check node sets.
    map<string, set<int> > node_sets = mesh.get_node_sets();
    string node_flag_name[7] = { "node_type/interior", "node_type/dudded",
			         "node_type/parent", "boundary/reflective",
				 "boundary/vacuum", "source/no_source",
				 "source/rad_source"};
    vector<set<int > > flag_nodes(7);
    flag_nodes[0].insert(0);
    flag_nodes[1].insert(1); flag_nodes[1].insert(2);
    flag_nodes[2].insert(3);
    flag_nodes[3].insert(0); flag_nodes[3].insert(1);
    flag_nodes[4].insert(2); flag_nodes[4].insert(3);
    flag_nodes[5].insert(0); flag_nodes[5].insert(1); flag_nodes[5].insert(2);
    flag_nodes[6].insert(3);
    bool got_node_sets = true;

    for (int i = 0; i < 7; i++)
        if (node_sets.find(node_flag_name[i])->second != flag_nodes[i])
	    got_node_sets = false;

    if (!got_node_sets)
    {
        fail(" Node Sets ") << "Node sets not obtained." << endl;
 	all_passed = false;
    }

    // Check Element sets.
    map<string, set<int> > element_sets = mesh.get_element_sets();
    string element_flag_name[7] = {"boundary/reflective", "boundary/vacuum",
				   "material/control_rod", "material/shield",
				   "source/src_name1","source/src_name2"};
    vector<set<int > > flag_elements(6);
    flag_elements[0].insert(1); flag_elements[0].insert(2);
        flag_elements[0].insert(3);
    flag_elements[1].insert(0);
    flag_elements[2].insert(4);
    flag_elements[5].insert(4);
    bool got_element_sets = true;

    for (int i = 0; i < 6; i++)
        if (element_sets.find(element_flag_name[i])->second != 
	                                                  flag_elements[i])
	    got_element_sets = false;

    if (!got_element_sets)
    {
        fail(" Element Sets ") << "Element sets not obtained." << endl;
 	all_passed = false;
    }

    // Check title.
    string title = mesh.get_title();
    string test_title = "RTT_format mesh file definition, version 7.";
    os() << "Mesh title = " << title  << endl;

    if (title != test_title)
    {
        fail(" Title ") << "Title not obtained." << endl;
 	all_passed = false;
    }

    // Check invariant.
    bool test_invariant = mesh.invariant();

    if (!test_invariant)
    {
	fail(" Invariant ") << "Invariant not satisfied." << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" Virtual Accessors " ) << "Got all virtual accessors." << endl;
    else
	fail(" Virtual Accessors ") << "Errors in some virtual accessors." 
				    << endl;

    return all_passed;
}

bool TestRTT_Format::check_header(const rtt_meshReaders::RTT_Format & mesh)
{
    // Exercize the header accessor functions for this mesh.
    bool all_passed = true;

    // Check version.
    if (mesh.get_header_version() != "v1.0.0")
    {
        fail(" Header Version ") << "Header version not obtained." << endl;
 	all_passed = false;
    }

    // Check title.
    if (mesh.get_header_title() != 
	"RTT_format mesh file definition, version 7.")
    {
        fail(" Header Title ") << "Header title not obtained." << endl;
 	all_passed = false;
    }

    // Check date.
    if (mesh.get_header_date() != "24 Jul 97")
    {
        fail(" Header Date ") << "Header date not obtained." << endl;
 	all_passed = false;
    }

    // Check cycle.
    if (mesh.get_header_cycle() != 1)
    {
        fail(" Header Cycle ") << "Header cycle not obtained." << endl;
 	all_passed = false;
    }

    // Check time.
    if (mesh.get_header_time() != 0.0)
    {
        fail(" Header Time ") << "Header time not obtained." << endl;
 	all_passed = false;
    }

    // Check ncomments.
    if (mesh.get_header_ncomments() != 3)
    {
        fail(" Header Ncomments ") << "Header ncomments not obtained." << endl;
 	all_passed = false;
    }

    // Check comments.
    string test_comments[3]  = {"One tet mesh in an RTT mesh file format.",
				"Date     : 24 Jul 97",
				"Author(s): H. Trease, J.McGhee"};
    bool got_comments = true;

    for (int i= 0; i < 3; i++)
        if (test_comments[i] != mesh.get_header_comments(i))
	    got_comments = false;

    if (!got_comments)
    {
        fail(" Header Comments ") << "Header comments not obtained." << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" Header Accessors " ) << "Got all Header accessors." << endl;
    else
	fail(" Header Accessors ") << "Errors in some Header accessors." 
				   << endl;

    return all_passed;
}

bool TestRTT_Format::check_dims(const rtt_meshReaders::RTT_Format & mesh)
{
    // Exercize the dims accessor functions for this mesh.
    bool all_passed = true;

    // Check coordinate units.
    if (mesh.get_dims_coor_units() != "cm")
    {
        fail(" Dims coor_units ") << 
	     "Dims coor_units not obtained." << endl;
 	all_passed = false;
    }

    // Check problem time units.
    if (mesh.get_dims_prob_time_units() != "s")
    {
        fail(" Dims prob_time_units ") << 
	     "Dims prob_time_units not obtained." << endl;
 	all_passed = false;
    }

    // Check number of cell definitions.
    if (mesh.get_dims_ncell_defs() != 8)
    {
        fail(" Dims ncell_defs ") << 
	     "Dims ncell_defs not obtained." << endl;
 	all_passed = false;
    }

    // Check .
    if (mesh.get_dims_nnodes_max() != 8)
    {
        fail(" Dims nnodes_max ") << 
	     "Dims nnodes_max not obtained." << endl;
 	all_passed = false;
    }

     if (mesh.get_dims_nsides_max() != 6)
    {
        fail(" Dims nsides_max ") << 
	     "Dims nsides_max not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_nnodes_side_max() != 4)
    {
        fail(" Dims nnodes_side_max ") << 
	     "Dims nnodes_side_max not obtained." << endl;
 	all_passed = false;
    }

     if (mesh.get_dims_ndim() != 3)
    {
        fail(" Dims ndim ") << 
	     "Dims ndim not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_ndim_topo() != 3)
    {
        fail(" Dims ndim_topo ") << 
	     "Dims ndim_topo not obtained." << endl;
 	all_passed = false;
    }

     if (mesh.get_dims_nnodes() != 4)
    {
        fail(" Dims nnodes ") << 
	     "Dims nnodes not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_nnode_flag_types() != 3)
    {
        fail(" Dims nnode_flag_types ") << 
	     "Dims nnode_flag_types not obtained." << endl;
 	all_passed = false;
    }

    int nnode_flags[3] = {3, 2, 2};
    bool got_nnode_flags = true;
    for (int f = 0; f < 3; f++)
        if (mesh.get_dims_nnode_flags(f) != nnode_flags[f])
	    got_nnode_flags = false;
    if (!got_nnode_flags)
    {
        fail(" Dims nnode_flags ") << "Dims nnode_flags not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_nnode_data() != 3)
    {
        fail(" Dims nnode_data ") << 
	     "Dims nnode_data not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_nsides() != 4)
    {
        fail(" Dims nsides ") << 
	     "Dims nsides not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_nside_types() != 1)
    {
        fail(" Dims nside_types ") << 
	     "Dims nside_types not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_side_types(0) != 2)
    {
        fail(" Dims side_types ") << 
	     "Dims side_types not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_nside_flag_types() != 1)
    {
        fail(" Dims nside_flag_types ") << 
	     "Dims nside_flag_types not obtained." << endl;
 	all_passed = false;
    }

    int nside_flags[1] = {2};
    bool got_nside_flags = true;
    for (int f = 0; f < 1; f++)
        if (mesh.get_dims_nside_flags(f) != nside_flags[f])
	    got_nside_flags = false;
    if (!got_nside_flags)
    {
        fail(" Dims nside_flags ") << "Dims nside_flags not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_nside_data() != 2)
    {
        fail(" Dims nside_data ") << 
	     "Dims nside_data not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_ncells() != 1)
    {
        fail(" Dims ncells ") << 
	     "Dims ncells not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_ncell_types() != 1)
    {
        fail(" Dims ncell_types ") << 
	     "Dims ncell_types not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_cell_types(0) != 5)
    {
        fail(" Dims cell_types ") << 
	     "Dims cell_types not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_ncell_flag_types() != 2)
    {
        fail(" Dims ncell_flag_types ") << 
	     "Dims ncell_flag_types not obtained." << endl;
 	all_passed = false;
    }

    int ncell_flags[2] = {2, 2};
    bool got_ncell_flags = true;
    for (int f = 0; f < 2; f++)
        if (mesh.get_dims_ncell_flags(f) != ncell_flags[f])
	    got_ncell_flags = false;
    if (!got_ncell_flags)
    {
        fail(" Dims ncell_flags ") << "Dims ncell_flags not obtained." << endl;
 	all_passed = false;
    }

    if (mesh.get_dims_ncell_data() != 1)
    {
        fail(" Dims ncell_data ") << 
	     "Dims ncell_data not obtained." << endl;
 	all_passed = false;
    }

 if (all_passed)
        pass(" Dims Accessors " ) << "Got all Dims accessors." << endl;
    else
	fail(" Dims Accessors ") << "Errors in some Dims accessors." 
				 << endl;

    return all_passed;
}

bool TestRTT_Format::check_node_flags(const rtt_meshReaders::RTT_Format & mesh)
{
    // Exercize the node_flags accessor functions for this mesh.
    bool all_passed = true;

    // Check node flag types.
    string node_flag_type[3] = { "node_type", "boundary", "source"};
    bool got_node_flag_types = true;

    for (int i = 0; i < 3; i++)
        if (node_flag_type[i] != mesh.get_node_flags_flag_type(i))
	    got_node_flag_types = false;
    if (!got_node_flag_types)
    {
        fail(" NodeFlags flag_type ") << "Node Flags flag_types not obtained."
				      << endl;
 	all_passed = false;
    }

    // Check node flag numbers for each of the flag types.
    int node_flag_number[7] = {11, 21, 6, 1, 4, 101, 22};
    int node_j = 0;
    bool got_node_flag_numbers = true;

    for (int i = 0; i < 3; i++)
    {
         for (int j = 0; j < mesh.get_dims_nnode_flags(i); j++) 
	     if (node_flag_number[j + node_j] != 
		 mesh.get_node_flags_flag_number(i, j))
	         got_node_flag_numbers = false;
	 node_j += mesh.get_dims_nnode_flags(i);
    }
    if (!got_node_flag_numbers)
    {
        fail(" NodeFlags flag_number ") << 
	    "Node Flags flag_numbers not obtained." << endl;
 	all_passed = false;
    }

    // Check number of flags for each node flag type.
    int node_flag_size[3] = { 3, 2, 2};
    bool got_node_flag_size = true;

    for (int i = 0; i < 3; i++)
        if (node_flag_size[i] != mesh.get_node_flags_flag_size(i))
	    got_node_flag_size = false;
    if (!got_node_flag_size)
    {
        fail(" NodeFlags flag_size ") << "Node Flags flag_size not obtained."
				      << endl;
 	all_passed = false;
    }

    // Check node flag names for each of the flag types.
    string node_flag_name[7] = {"interior", "dudded", "parent", 
				"reflective", "vacuum", "no_source", 
				"rad_source"};
    node_j = 0;
    bool got_node_flag_name = true;

    for (int i = 0; i < 3; i++)
    {
         for (int j = 0; j < mesh.get_dims_nnode_flags(i); j++) 
	     if (node_flag_name[j + node_j] != 
		 mesh.get_node_flags_flag_name(i, j))
	         got_node_flag_name = false;
	 node_j += mesh.get_dims_nnode_flags(i);
    }
    if (!got_node_flag_name)
    {
        fail(" NodeFlags flag_name ") << 
	    "Node Flags flag_name not obtained." << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" NodeFlags Accessors " ) << "Got all NodeFlags accessors." 
				       << endl;
    else
	fail(" NodeFlags Accessors ") << "Errors in some NodeFlags accessors." 
				      << endl;

    return all_passed;
}

} // end namespace rtt_meshReaders_test


//---------------------------------------------------------------------------//
//                              end of TestRTT_Format.cc
//---------------------------------------------------------------------------//
