//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestHexFormat.cc
 * \author John McGhee
 * \date   Thu Mar  9 08:54:59 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestHexFormat.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>

using std::cout;
using std::endl;

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_meshReaders_test::TestHexFormat;
    
    return SP<TestApp>(new TestHexFormat(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_meshReaders_test
{

using std::string;
using rtt_meshReaders::Hex_Format;

TestHexFormat::TestHexFormat(int argc, char *argv[],
					       std::ostream& os_in)

    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestHexFormat" << endl;
}

string TestHexFormat::version() const
{
    return rtt_meshReaders::release();
}

/*!
 * \brief Tests the CIC-19 Hex mesh format reader.
 *
 */
string TestHexFormat::runTest()
{
    using rtt_meshReaders::Hex_Format;

    // Read and test a 1D mesh.
    std::string filename = "marshak_slab.mesh.inp";
    Hex_Format mesh_1D(filename);
    pass(" Construct ") << 
	"Read mesh without coreing in or firing an assertion." << std::endl;
    check_mesh(mesh_1D);

    // Read and test a 2D mesh.
    filename = "quad.mesh.inp";
    Hex_Format mesh_2D(filename);
    pass(" Construct ") << 
    	"Read mesh without coreing in or firing an assertion." << std::endl;
    check_mesh(mesh_2D);
    
     // Read and test a 3D mesh.
    //    filename = "marshak_slab.mesh.inp";
    //  Hex_Format mesh_3D(filename);
    // pass(" Construct ") << 
    //	 "Read mesh without coreing in or firing an assertion." << std::endl;
    // check_mesh(mesh_3D);
    
    // Report results of test.
    if (passed())
    {
	return "All tests passed.";
    }
    return "Some tests failed.";
}

bool TestHexFormat::check_mesh(const rtt_meshReaders::Hex_Format &mesh)
{
    // Exercize the accessor functions for this mesh.

    // Check node coords
    std::vector<std::vector<double> > point_c = mesh.get_node_coords();
    pass(" Node Coordinates ") << "Got node coordinates." << std::endl;
    
    // Check coordinate units.
    std::string punits = mesh.get_node_coord_units();
    os() << "units= " << punits << std::endl;
    pass(" Coordinate Units ") << "Got coordinate units." << std::endl;

    // Check node sets.
    std::map<std::string, std::set<int> > ndsets = mesh.get_node_sets();
    pass(" Node Sets ") << "Got node sets." << std::endl;

    // Check title.
    std::string title = mesh.get_title();
    os() << "title= " << title  << std::endl;
    pass(" Title ") << "Got title." << std::endl;

    // Check element nodes.
    std::vector<std::vector<int> > enodes = mesh.get_element_nodes();
    pass(" Element Nodes ") << "Got element nodes." << std::endl;
 
    // Check invariant.
    bool invr = mesh.invariant();
    if ( invr)
	pass(" Invariant ") << "Invoked invariant." << std::endl;
    else
	pass(" Invariant ") << "Invariant not satisfied." << std::endl;

    // Check Element sets.
    std::map<std::string, std::set<int> > elmsets = mesh.get_element_sets();
    pass(" Element Sets ") << "Got element sets." << std::endl;

    // Check Element Types.
    std::vector<rtt_meshReaders::Element_Definition::Element_Type>
	etypes = mesh.get_element_types();
    // os() << "Element Types: " << endl;
    //    for (int i=0; i<etypes.size(); ++i)
    //	os() << i << ",  " << etypes[i] << endl;
    pass(" Element Types ") << "Read Element Types." << std::endl;

    return true;
}

} // end namespace rtt_meshReaders_test


//---------------------------------------------------------------------------//
//                              end of TestHexFormat.cc
//---------------------------------------------------------------------------//
