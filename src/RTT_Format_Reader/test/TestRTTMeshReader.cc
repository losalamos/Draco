//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTMeshReader.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 27 10:41:12 2002
 * \brief  RTT_Mesh_Reader test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "RTT_Format_Reader_test.hh"
#include "TestRTTMeshReader.hh"
#include "../Release.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

using rtt_RTT_Format_Reader::RTT_Mesh_Reader;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void runTest()
{
    // New meshes added to this test will have to be added to the enumeration
    // Meshes in the header file.
    const int MAX_MESHES = 1;
    std::string filename[MAX_MESHES] = {"rttdef.mesh"};
    rtt_RTT_Mesh_Reader_test::Meshes mesh_type;

    for (int mesh_number = 0; mesh_number < MAX_MESHES; mesh_number++)
    {
        // Construct an RTT_Mesh_Reader class object from the data in the 
	// specified mesh file. 
        RTT_Mesh_Reader mesh(filename[mesh_number]);
	{
	    ostringstream m;
	    m << "Read " << filename[mesh_number] 
	      << " without coreing in or firing an assertion." 
	      << std::endl;
	    PASSMSG(m.str());
	}
	bool all_passed = true;
	// The following switch allows addition of other meshes for testing,
	// with the "DEFINED" mesh providing an example.
        switch (mesh_number)
	{
	    // Test all nested class accessor functions for a very simplistic 
	    // mesh file (enum DEFINED).
	case (0):
	    mesh_type = rtt_RTT_Mesh_Reader_test::DEFINED;
	    all_passed = all_passed && rtt_RTT_Mesh_Reader_test::check_virtual(mesh, mesh_type);
	    break;

	default:
	    FAILMSG("Invalid mesh type encountered.");
	    all_passed = false;
	    break;
	}
	if (!all_passed)
	{
	    ostringstream m;
	    m << "Errors occured testing mesh " 
	      << "number " << mesh_type << std::endl;
	    FAILMSG(m.str());
	}
    }

    // Report results of test.
    if (rtt_RTT_Format_Reader_test::passed)
    {
	PASSMSG("All tests passed.");
    }
    else
    {
	FAILMSG("Some tests failed.");
    }
}

//---------------------------------------------------------------------------//

namespace rtt_RTT_Mesh_Reader_test
{

bool check_virtual(const RTT_Mesh_Reader & mesh, const Meshes & meshtype)
{
    // Exercise the virtual accessor functions for this mesh.
    bool all_passed = true;
    std::vector<std::vector<double> > node_coords;
    std::string node_coord_units;
    std::vector<std::vector<int> > element_nodes;
    std::vector<rtt_meshReaders::Element_Definition::Element_Type> 
        element_types;
    std::vector<rtt_meshReaders::Element_Definition::Element_Type> 
        unique_element_types;
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
	// load the unique element types defined for RTT_Format according to 
	// the corresponding Element_Definition::Element_Type.
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::NODE); 
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::BAR_2); 
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::TRI_3); 
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::QUAD_4); 
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::PYRA_5); 
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::TETRA_4); 
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::PENTA_6); 
	unique_element_types.push_back(
	    rtt_meshReaders::Element_Definition::HEXA_8); 
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
	FAILMSG("Invalid mesh type encountered.");
	all_passed = false;
	return all_passed;
    }
    // Check node coords
    if (node_coords != mesh.get_node_coords())
    {
	FAILMSG("Node coordinates not obtained.");
	all_passed = false;
    }
    // Check coordinate units.
    if (node_coord_units != mesh.get_node_coord_units())
    {
	FAILMSG("Coordinate units not obtained.");
 	all_passed = false;
    }
    if (element_nodes != mesh.get_element_nodes())
    {
	FAILMSG("Element nodes not obtained.");
 	all_passed = false;
    }
    // Check Element Types.
    if (element_types != mesh.get_element_types())
    {
	FAILMSG("Element Types not obtained.");
 	all_passed = false;
    }
    // Check Unique Element Types.
    if (unique_element_types != mesh.get_unique_element_types())
    {
	FAILMSG("Unique Element Types not obtained.");
 	all_passed = false;
    }
    // Check node sets.
    if (node_sets != mesh.get_node_sets())
    {
	FAILMSG("Node sets not obtained.");
 	all_passed = false;
    }
    // Check Element sets.
    if (element_sets != mesh.get_element_sets())
    {
	FAILMSG("Element sets not obtained.");
 	all_passed = false;
    }
    // Check title.
    if (title != mesh.get_title())
    {
	FAILMSG("Title not obtained.");
 	all_passed = false;
    }
    // Check invariant.
    if (!mesh.invariant())
    {
	FAILMSG("Invariant not satisfied.");
 	all_passed = false;
    }
    if (all_passed)
    {
        PASSMSG("Got all virtual accessors.");
    }
    else
    {
	FAILMSG("Errors in some virtual accessors.");
    }

    return all_passed;
}

} // end of rtt_RTT_Format_Reader_test

//---------------------------------------------------------------------------//

void reportVersion()
{
    using std::cout;
    using std::endl;
    cout << "\nThis is RTT_Format_Reader\n" 
	 << "Version "<< rtt_RTT_Format_Reader::release()
	 << "\n==================================================\n" << endl;
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_RTT_Format_Reader::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	reportVersion();
	runTest();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing TestRTTMeshReader, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_RTT_Format_Reader_test::passed) 
    {
        cout << "**** TestRTTMeshReader Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing TestRTTMeshReader." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of TestRTTMeshReader.cc
//---------------------------------------------------------------------------//
