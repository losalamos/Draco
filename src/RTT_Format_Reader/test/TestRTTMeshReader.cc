//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTMeshReader.cc
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Mesh_Reader class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestRTTMeshReader.hh"
#include "meshReaders/Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_RTT_Mesh_Reader_test::TestRTT_Mesh_Reader;
    
    return SP<TestApp>(new TestRTT_Mesh_Reader(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_RTT_Mesh_Reader_test
{

TestRTT_Mesh_Reader::TestRTT_Mesh_Reader(int argc, char * argv[], 
					     std::ostream & os_in)
    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestRTT_Mesh_Reader" << std::endl;
}

std::string TestRTT_Mesh_Reader::version() const
{
    return rtt_meshReaders::release();
}

/*!
 * \brief Tests the RTT_Mesh_Reader mesh reader.
 *
 */
std::string TestRTT_Mesh_Reader::runTest()
{
    // New meshes added to this test will have to be added to the enumeration
    // Meshes in the header file.
    const int MAX_MESHES = 1;
    std::string filename[MAX_MESHES] = {"rttdef.mesh"};
    Meshes mesh_type;

    for (int mesh_number = 0; mesh_number < MAX_MESHES; mesh_number++)
    {
        // Construct an RTT_Mesh_Reader class object from the data in the 
	// specified mesh file. 
        RTT_Mesh_Reader mesh(filename[mesh_number]);
	pass(" Construct ") << "Read " << filename[mesh_number] 
			    << " without coreing in or firing an assertion." 
			    << std::endl;
	bool all_passed = true;
	// The following switch allows addition of other meshes for testing,
	// with the "DEFINED" mesh providing an example.
        switch (mesh_number)
	{
	// Test all nested class accessor functions for a very simplistic 
	// mesh file (enum DEFINED).
	case (0):
	    mesh_type = DEFINED;
	    all_passed = all_passed && check_virtual(mesh, mesh_type);
	    break;

	default:
	    fail("runTest") << "Invalid mesh type encountered." << std::endl;
	    all_passed = false;
	    break;
	}
	if (!all_passed)
	    fail(filename[mesh_number]) << "Errors occured testing mesh " 
					<< "number " << mesh_type << std::endl;
    }

    // Report results of test.
    if (passed())
    {
	return "All tests passed.";
    }
    return "Some tests failed.";
}
bool TestRTT_Mesh_Reader::check_virtual(const RTT_Mesh_Reader & mesh, 
					const Meshes & meshtype)
{
    // Exercise the virtual accessor functions for this mesh.
    bool all_passed = true;
    std::vector<std::vector<double> > node_coords;
    std::string node_coord_units;
    std::vector<std::vector<int> > element_nodes;
    std::vector<rtt_meshReaders::Element_Definition::Element_Type> 
        element_types;
    std::map<std::string, std::set<int> > node_sets;
    std::map<std::string, std::set<int> > element_sets;
    std::string title;
    std::vector<double> coords(3,0.0);
    std::vector<int> side_nodes;
    std::set<int> flag_nodes;
    std::set<int> flag_elements;

    switch (meshtype)
    {
    case DEFINED:
        // set node coords per the input deck.
	node_coords.push_back(coords);
	coords[2] = 3.0;
	node_coords.push_back(coords);
	coords[1] = 2.0; coords[2] = 0.0;
	node_coords.push_back(coords);
	coords[0] = 1.0; coords[1] = 0.0;
	node_coords.push_back(coords);
	// set the coordinate units used for the nodes.
	node_coord_units = "cm";
	// load the node numbers for the single tet cell defined in the input
	// file (note that the node numbers are zero indexed).
	side_nodes.push_back(1); side_nodes.push_back(2);
	    side_nodes.push_back(3);
	element_nodes.push_back(side_nodes);
	side_nodes.resize(0);
	side_nodes.push_back(0); side_nodes.push_back(3);
	    side_nodes.push_back(2);
	element_nodes.push_back(side_nodes);
	side_nodes.resize(0);
	side_nodes.push_back(0); side_nodes.push_back(1);
	    side_nodes.push_back(3);
	element_nodes.push_back(side_nodes);
	side_nodes.resize(0);
	side_nodes.push_back(0); side_nodes.push_back(2);
	    side_nodes.push_back(1);
	element_nodes.push_back(side_nodes);
	side_nodes.resize(0);
	side_nodes.push_back(0); side_nodes.push_back(1); 
	    side_nodes.push_back(2); side_nodes.push_back(3);
	element_nodes.push_back(side_nodes);
	side_nodes.resize(0);
	// load the element types defined for RTT_Format according to the
	// corresponding Element_Definition::Element_Type.
	element_types.push_back(rtt_meshReaders::Element_Definition::TRI_3); 
	element_types.push_back(rtt_meshReaders::Element_Definition::TRI_3); 
	element_types.push_back(rtt_meshReaders::Element_Definition::TRI_3); 
	element_types.push_back(rtt_meshReaders::Element_Definition::TRI_3); 
	element_types.push_back(rtt_meshReaders::Element_Definition::TETRA_4);
	// load the node sets
	flag_nodes.insert(0);
	node_sets.insert(std::make_pair(std::string("node_type/interior"), 
					flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(1); flag_nodes.insert(2);
	node_sets.insert(std::make_pair(std::string("node_type/dudded"), 
					flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(3);
	node_sets.insert(std::make_pair(std::string("node_type/parent"), 
					flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(0); flag_nodes.insert(1);
	node_sets.insert(std::make_pair(std::string("boundary/reflective"), 
					flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(2); flag_nodes.insert(3);
	node_sets.insert(std::make_pair(std::string("boundary/vacuum"), 
					flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(0); flag_nodes.insert(1); flag_nodes.insert(2);
	node_sets.insert(std::make_pair(std::string("source/no_source"),
					flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(3);
	node_sets.insert(std::make_pair(std::string("source/rad_source"), 
					flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	// load the element (i.e., sides + cell) sets
	flag_elements.insert(1); flag_elements.insert(2); 
	    flag_elements.insert(3);
	element_sets.insert(std::make_pair(std::string("boundary/reflective"),
					   flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	flag_elements.insert(0);
	element_sets.insert(std::make_pair(std::string("boundary/vacuum"), 
					   flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	flag_elements.insert(4);
	element_sets.insert(std::make_pair(std::string("material/control_rod"),
					   flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	element_sets.insert(std::make_pair(std::string("material/shield"), 
					   flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	element_sets.insert(std::make_pair(std::string("rad_source/src_name1"),
					   flag_elements));
	flag_elements.insert(4);
	element_sets.insert(std::make_pair(std::string("rad_source/src_name2"),
					   flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	// set the mesh title
	title = "RTT_format mesh file definition, version 7.";
	break;

    default:
        fail("check_virtual") << "Invalid mesh type encountered." << std::endl;
	all_passed = false;
	return all_passed;
    }
    // Check node coords
    if (node_coords != mesh.get_node_coords())
    {
        fail(" Node Coordinates ") << "Node coordinates not obtained." 
				   << std::endl;
	all_passed = false;
    }
    // Check coordinate units.
    if (node_coord_units != mesh.get_node_coord_units())
    {
        fail(" Coordinates Units ") << "Coordinate units not obtained." 
				    << std::endl;
 	all_passed = false;
    }
    if (element_nodes != mesh.get_element_nodes())
    {
        fail(" Element Nodes ") << "Element nodes not obtained." << std::endl;
 	all_passed = false;
    }
    // Check Element Types.
    if (element_types != mesh.get_element_types())
    {
	fail(" Element Types ") << "Element Types not obtained." << std::endl;
 	all_passed = false;
    }
    // Check node sets.
    if (node_sets != mesh.get_node_sets())
    {
        fail(" Node Sets ") << "Node sets not obtained." << std::endl;
 	all_passed = false;
    }
    // Check Element sets.
    if (element_sets != mesh.get_element_sets())
    {
        fail(" Element Sets ") << "Element sets not obtained." << std::endl;
 	all_passed = false;
    }
    // Check title.
    if (title != mesh.get_title())
    {
        fail(" Title ") << "Title not obtained." << std::endl;
 	all_passed = false;
    }
    // Check invariant.
    if (!mesh.invariant())
    {
	fail(" Invariant ") << "Invariant not satisfied." << std::endl;
 	all_passed = false;
    }
    if (all_passed)
        pass(" Virtual Accessors " ) << "Got all virtual accessors." 
				     << std::endl;
    else
	fail(" Virtual Accessors ") << "Errors in some virtual accessors." 
				    << std::endl;

    return all_passed;
}

} // end namespace rtt_RTT_Mesh_Reader_test


//---------------------------------------------------------------------------//
//                end of RTT_Format_Reader/test/TestRTTMeshReader.cc
//---------------------------------------------------------------------------//
