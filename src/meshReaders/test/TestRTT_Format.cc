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
    using rtt_format_test::TestRTT_Format;
    
    return SP<TestApp>(new TestRTT_Format(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_format_test
{

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::set;
using std::ostream;
using std::fill;
using rtt_format::RTT_Format;

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
    check_virtual(mesh);
    check_header(mesh);

    // Report results of test.
    if (passed())
    {
	return "All tests passed.";
    }
    return "Some tests failed.";
}

bool TestRTT_Format::check_virtual(const rtt_format::RTT_Format & mesh)
{
    // Exercize the virtual accessor functions for this mesh.

    // Check node coords
    vector<vector<double> > coords = mesh.get_node_coords();
    vector<vector<double> > test_coords(4, vector<double> (3, 0.0));
    for (int i = 1; i < 4; i++)
        test_coords[i][i - 1] = static_cast<double>(i);

    if (test_coords == coords)
        pass(" Node Coordinates ") << "Got node coordinates." << endl;
    else
        fail(" Node Coordinates ") << "Node coordinates not obtained." << endl;

    // Check coordinate units.
    string coord_units = mesh.get_node_coord_units();
    string test_coord_units = "cm";
    os() << "Coordinate Units= " << coord_units << endl;

    if (test_coord_units == coord_units)
        pass(" Coordinate Units ") << "Got coordinate units." << endl;
    else
        fail(" Coordinates Units ") << "Coordinate units not obtained." 
				    << endl;
      
    // Check element nodes.
    vector<vector<int> > element_nodes = mesh.get_element_nodes();

    const int side_nodes = 3;
    int test_elem_0[side_nodes] = {1, 2, 3};
    int test_elem_1[side_nodes] = {0, 3, 2};
    int test_elem_2[side_nodes] = {0, 1, 3};
    int test_elem_3[side_nodes] = {0, 2, 1};

    const int cell_nodes = 4;
    int test_elem_4[cell_nodes] = {0, 1, 2, 3};

    if (element_nodes[0] == vector<int>(test_elem_0,test_elem_0+side_nodes) &&
	element_nodes[1] == vector<int>(test_elem_1,test_elem_1+side_nodes) &&
	element_nodes[2] == vector<int>(test_elem_2,test_elem_2+side_nodes) &&
	element_nodes[3] == vector<int>(test_elem_3,test_elem_3+side_nodes) &&
	element_nodes[4] == vector<int>(test_elem_4,test_elem_4+cell_nodes))
        pass(" Element Nodes ") << "Got element nodes." << endl;
    else
        fail(" Element Nodes ") << "Element nodes not obtained." << endl;

    // Check Element Types.
    vector<rtt_meshReaders::Element_Definition::Element_Type>
	element_types = mesh.get_element_types();
    int test_element_types[8] = {0, 1, 3, 5, 10, 8, 12, 15};
    bool got_element_types = true;

    for (int i = 0; i < 8; i++)
        if (element_types[i] != test_element_types[i])
	    got_element_types = false;

    if (got_element_types)
	pass(" Element Types ") << "Read Element Types." << endl;
    else
        fail(" Element Types ") << "Element Types not obtained." << endl;

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

    if (got_node_sets)
        pass(" Node Sets ") << "Got node sets." << endl;
    else
        fail(" Node Sets ") << "Node sets not obtained." << endl;

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

    if (got_element_sets)
        pass(" Element Sets ") << "Got element sets." << endl;
    else
        fail(" Element Sets ") << "Element sets not obtained." << endl;

    // Check title.
    string title = mesh.get_title();
    string test_title = "RTT_format mesh file definition, version 7.";
    os() << "Mesh title = " << title  << endl;

    if (title == test_title)
        pass(" Title ") << "Got title." << endl;
    else
        fail(" Title ") << "Title not obtained." << endl;

    // Check invariant.
    bool test_invariant = mesh.invariant();
    if (test_invariant)
	pass(" Invariant ") << "Invoked invariant." << endl;
    else
	pass(" Invariant ") << "Invariant not satisfied." << endl;

     return true;
}

bool TestRTT_Format::check_header(const rtt_format::RTT_Format & mesh)
{
    // Check version.
    if (mesh.get_header_version() == "v1.0.0")
        pass(" Header Version ") << "Got header version." << endl;
    else
        fail(" Header Version ") << "Header version not obtained." << endl;

    // Check title.
    if (mesh.get_header_title() == 
	"RTT_format mesh file definition, version 7.")
        pass(" Header Title ") << "Got header title." << endl;
    else
        fail(" Header Title ") << "Header title not obtained." << endl;

    // Check date.
    if (mesh.get_header_date() == "24 Jul 97")
        pass(" Header Date ") << "Got header date." << endl;
    else
        fail(" Header Date ") << "Header date not obtained." << endl;

    // Check cycle.
    if (mesh.get_header_cycle() == 1)
        pass(" Header Cycle ") << "Got header cycle." << endl;
    else
        fail(" Header Cycle ") << "Header cycle not obtained." << endl;

    // Check time.
    if (mesh.get_header_time() == 0.0)
        pass(" Header Time ") << "Got header time." << endl;
    else
        fail(" Header Time ") << "Header time not obtained." << endl;

    // Check ncomments.
    if (mesh.get_header_ncomments() == 3)
        pass(" Header Ncomments ") << "Got header ncomments." << endl;
    else
        fail(" Header Ncomments ") << "Header ncomments not obtained." << endl;

    // Check comments.
    string test_comments[3]  = {"One tet mesh in an RTT mesh file format.",
				"Date     : 24 Jul 97",
				"Author(s): H. Trease, J.McGhee"};
    bool got_comments = true;

    for (int i= 0; i < 3; i++)
        if (test_comments[i] != mesh.get_header_comments(i))
	    got_comments = false;

    if (got_comments)
        pass(" Header Comments ") << "Got header comments." << endl;
    else
        fail(" Header Comments ") << "Header comments not obtained." << endl;

    return true;
}

} // end namespace rtt_format_test


//---------------------------------------------------------------------------//
//                              end of TestRTT_Format.cc
//---------------------------------------------------------------------------//
