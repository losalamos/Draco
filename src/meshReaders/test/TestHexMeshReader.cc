//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestHexMeshReader.cc
 * \author John McGhee
 * \date   Thu Mar  9 08:54:59 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestHexMeshReader.hh"
#include "../Release.hh"
#include "../Element_Definition.hh"
#include "../Hex_Mesh_Reader.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_meshReaders_test::TestHexMeshReader;
    
    return SP<TestApp>(new TestHexMeshReader(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_meshReaders_test
{

using std::string;
using rtt_meshReaders::Hex_Mesh_Reader;

TestHexMeshReader::TestHexMeshReader(int argc, char *argv[],
					       std::ostream& os_in)

    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestHexMeshReader" << endl;
}

string TestHexMeshReader::version() const
{
    return rtt_meshReaders::release();
}

/*!
 * \brief Tests the CIC-19 Hex mesh format reader.
 *
 */
string TestHexMeshReader::runTest()
{
    using rtt_meshReaders::Hex_Mesh_Reader;

    os() << std::endl << "******* CIC-19 Hex Mesh Reader Tests *******" 
	 << std::endl;
    // Read and test a 1D mesh.
    std::string filename = "slab.mesh";
    Hex_Mesh_Reader mesh_1D(filename);
    pass(" Construct 1D") << 
	"Read a 1D mesh without coreing in or firing an assertion." 
			  << std::endl;
    check_mesh(mesh_1D, "slab");

    // Read and test a 2D mesh.
    filename = "quad.mesh";
    Hex_Mesh_Reader mesh_2D(filename);
    pass(" Construct 2D ") << 
    	"Read a 2D mesh without coreing in or firing an assertion." 
			   << std::endl;
    check_mesh(mesh_2D, "quad");

    // Read and test a 3D mesh.
    filename = "cube.mesh";
    Hex_Mesh_Reader mesh_3D(filename);
    pass(" Construct 3D ") << 
	"Read a 3D mesh without coreing in or firing an assertion." 
			   << std::endl;
    check_mesh(mesh_3D, "cube");
    os() << endl;
   
    // Report results of test.
    if (passed())
    {
	return "All tests passed.";
    }
    return "Some tests failed.";
}

bool TestHexMeshReader::check_mesh(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh, const std::string &testid)
{

    // Exercize the accessor functions for this mesh and spot check
    // the results.

    bool pass_nc = check_nodes(mesh, testid);
    bool pass_cu = check_node_units(mesh);
    bool pass_ns = check_node_sets(mesh);
    bool pass_ti = check_title(mesh);
    bool pass_en = check_element_nodes(mesh, testid);
    bool pass_in = check_invariant(mesh);
    bool pass_es = check_element_sets(mesh, testid);
    bool pass_et = check_element_types(mesh, testid);

    return pass_nc && pass_cu && pass_ns && pass_ti && pass_en
	&& pass_in && pass_es && pass_et;
}

bool TestHexMeshReader::check_nodes(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh, const std::string &testid)
{
    // Check node coords -- Need to do a "fuzzy" check here, these are doubles!
    bool pass_nc = true;
    std::vector<std::vector<double> > point_c = mesh.get_node_coords();
     if (testid == "slab") 
 	pass_nc = 
 	    compare_double(point_c[0][0]  , 0. ) && 
	    compare_double(point_c[10][0] , 2.5) && 
	    compare_double(point_c[100][0], 25.);
     else if (testid == "quad")
 	pass_nc = 
 	    compare_double(point_c[0][0]  , 0.)   && 
	    compare_double(point_c[10][0] , 12.5) && 
	    compare_double(point_c[440][0], 25.)  &&
 	    compare_double(point_c[0][1]  , 0.)   && 
	    compare_double(point_c[10][1] , 0.)   && 
	    compare_double(point_c[440][1], 25.);
     else if (testid == "cube") 
     {
	 bool pass1 =
	     compare_double(point_c[0][0]  , 0.)  &&
	     compare_double(point_c[0][1]  , 0.)  && 
	     compare_double(point_c[0][2]  , 0.)  ;
	 bool pass2=  
	     compare_double(point_c[10][0] , 0.8) && 
	     compare_double(point_c[10][1] , 0.2) && 
	     compare_double(point_c[10][2] , 0.)  ; 
	 bool pass3=
	     compare_double(point_c[215][0], 1.) &&
	     compare_double(point_c[215][1], 1.) &&
	     compare_double(point_c[215][2], 1.);
	 pass_nc = pass1 && pass2 && pass3;
     }
     else
 	Insist(false,"Unrecognized test id string!");
     if (pass_nc) 
 	pass(" Node Coordinates ") << "Got node coordinates." << std::endl;
     else
 	fail(" Node Coordinates ") << "Error in node coordinates." << 
	    std::endl;	    
     return pass_nc;
}

bool TestHexMeshReader::check_node_units(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh)
{     
    // Check coordinate units.
    std::string punits = mesh.get_node_coord_units();
    bool pass_cu = punits == "unknown";
    if (pass_cu)
	pass(" Coordinate Units ") << "Got coordinate units." << std::endl;
    else
	fail(" Coordinate Units ") << "Error in coordinate units." <<
	    std::endl;
    return pass_cu;
}

bool TestHexMeshReader::check_node_sets(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh, const std::string &testid)
{
    // Check node sets.
    std::map<std::string, std::set<int> > ndsets = mesh.get_node_sets();
    bool pass_ns = ndsets.size() == 1;
    if (testid == "slab")
    {
	pass_ns = pass_ns && check_map(ndsets,"Interior",0,100);
    }
    else if (testid == "quad")
    {
	pass_ns = pass_ns && check_map(ndsets,"Interior",0,440);
    }
    else if (testid == "cube")
    {
	pass_ns = pass_ns && check_map(ndsets,"Interior",0,215);
    }
    else
	Insist(false,"Unrecognized test id string!");
    if (pass_ns)
	pass(" Node Sets ") << "Got node sets." << std::endl;
    else
	fail(" Node Sets ") << "Error in node sets." << std::endl;
    return pass_ns;
}

bool TestHexMeshReader::check_title(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh)
{
    // Check title.
    std::string title = mesh.get_title();
    bool pass_ti = title == "Untitled -- CIC-19 Hex Mesh";
    if (pass_ti)
	pass(" Title ") << "Got title." << std::endl;
    else
	fail(" Title ") << "Error in title." << std::endl;
    return pass_ti;
}

bool TestHexMeshReader::check_element_nodes(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh, const std::string &testid)
{
    // Check element nodes.
    bool pass_en = true;
    std::vector<std::vector<int> > enodes = mesh.get_element_nodes();
    if (testid == "slab")
	pass_en = 
	    enodes[0][0]   == 0   && enodes[0][1]   == 1  &&
	    enodes[10][0]  == 10  && enodes[10][1]  == 11 &&
	    enodes[99][0]  == 99  && enodes[99][1]  == 100 &&
	    enodes[100][0] == 0   && enodes[101][0] == 100 ;
    else if (testid == "quad")
	pass_en = 
	    enodes[0][0]   == 0   && enodes[0][1]   == 1   &&
	    enodes[0][2]   == 22  && enodes[0][3]   == 21  &&
	    enodes[10][0]  == 10  && enodes[10][1]  == 11  &&
	    enodes[10][2]  == 32  && enodes[10][3]  == 31  &&
	    enodes[399][0] == 418 && enodes[399][1] == 419 &&
	    enodes[399][2] == 440 && enodes[399][3] == 439 &&
	    enodes[400][0] == 421 && enodes[400][1] == 420 &&
	    enodes[419][0] == 440 && enodes[419][1] == 439 &&
 	    enodes[420][0] == 21  && enodes[420][1] == 0   && 
 	    enodes[479][0] == 19  && enodes[479][1] == 20; 
    else if (testid == "cube") 
    {
	bool pass1 =
	    enodes[0][0]   == 0   && enodes[0][1]   == 1   && 
	    enodes[0][2]   == 7   && enodes[0][3]   == 6   && 
	    enodes[0][4]   == 36  && enodes[0][5]   == 37  &&
	    enodes[0][6]   == 43  && enodes[0][7]   == 42  ;
	bool pass2 =
	    enodes[10][0]  == 12  && enodes[10][1]  == 13  &&
	    enodes[10][2]  == 19  && enodes[10][3]  == 18  &&
	    enodes[10][4]  == 48  && enodes[10][5]  == 49  &&
	    enodes[10][6]  == 55  && enodes[10][7]  == 54  ;
	bool pass3 =
	    enodes[124][0] == 172 && enodes[124][1] == 173 &&
	    enodes[124][2] == 179 && enodes[124][3] == 178 &&
	    enodes[124][4] == 208 && enodes[124][5] == 209 &&
	    enodes[124][6] == 215 && enodes[124][7] == 214 ;
	bool pass4 =
	    enodes[125][0] == 5   && enodes[125][1] == 11  &&
	    enodes[125][2] == 47  && enodes[125][3] == 41  ;
	bool pass5 =
	    enodes[249][0] == 208 && enodes[249][1] == 209 &&
	    enodes[249][2] == 215 && enodes[249][3] == 214 ;
	bool pass6 =
	    enodes[250][0] == 0   && enodes[250][1] == 36  &&
	    enodes[250][2] == 42  && enodes[250][3] == 6   ;
	bool pass7 =
	    enodes[274][0] == 168 && enodes[274][1] == 204 &&
	    enodes[274][2] == 210 && enodes[274][3] == 174 ;
	pass_en = pass1 && pass2 && pass3 && pass4 && pass5 && 
	    pass6 && pass7; 
    }
    else
	Insist(false,"Unrecognized test id string!");
    if (pass_en) 
	pass(" Element Nodes ") << "Got element nodes." << std::endl;
    else
	fail(" Element Nodes ") << "Error in element nodes." << std::endl;
    return pass_en;
}

bool TestHexMeshReader::check_invariant(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh)
{ 
    // Check invariant.
    bool invr = mesh.invariant();
    if ( invr)
	pass(" Invariant ") << "Invoked invariant." << std::endl;
    else
	fail(" Invariant ") << "Invariant not satisfied." << std::endl;
    return invr;
}

bool TestHexMeshReader::check_element_sets(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh, const std::string &testid)
{
    typedef std::map<std::string, std::set<int> > mt;
    bool pass_es = true;
    const mt elmsets = mesh.get_element_sets();
    if (testid == "slab") 
    {
	pass_es == pass_es && elmsets.size() == 4;
	pass_es = pass_es && check_map(elmsets,"Interior",0,100);
	pass_es = pass_es && check_map(elmsets,"Interior_Region_1",0,100);
	pass_es = pass_es && check_map(elmsets,"Vacuum_Boundary",100,102);
	pass_es = pass_es && check_map(elmsets,"Vacuum_Boundary_Region_1",
				       100,102);
    }
    else if (testid == "quad")
    {
	pass_es = pass_es && elmsets.size() == 6;
	pass_es = pass_es && check_map(elmsets,"Interior",0,400);
	pass_es = pass_es && check_map(elmsets,"Interior_Region_1",20,400);
	pass_es = pass_es && check_map(elmsets,"Interior_Region_2",0,20);
	pass_es = pass_es && check_map(elmsets,"Vacuum_Boundary",400,420);
	pass_es = pass_es && check_map(elmsets,"Vacuum_Boundary_Region_1",
				       400,420);
	pass_es = pass_es && check_map(elmsets,"Reflective_Boundary",
				       420,480);

    }
    else if (testid == "cube") 
    {
	bool pass0 = elmsets.size() == 10;
	bool pass1 = check_map(elmsets,"Interior",0,125);
	bool pass2 = check_map(elmsets,"Interior_Region_10",0,2);
	bool pass3 = check_map(elmsets,"Interior_Region_14",2,3);
	bool pass4 = check_map(elmsets,"Interior_Region_1",3,122);
	bool pass5 = check_map(elmsets,"Interior_Region_2",122,125);
	bool pass6 = check_map(elmsets,"Vacuum_Boundary",125,250);
	bool pass7 = check_map(elmsets,"Vacuum_Boundary_Region_3",125,130);
	bool pass8 = check_map(elmsets,"Vacuum_Boundary_Region_4",130,131);
	bool pass9 = check_map(elmsets,"Vacuum_Boundary_Region_1",131,250);
	bool pass10= check_map(elmsets,"Reflective_Boundary",250,275);
	pass_es = pass0 && pass1 && pass2 && pass3 && pass4 && pass5
	    && pass6 && pass7 && pass8 && pass9 && pass10;
    }
    else
	Insist(false,"Unrecognized test id string!");
    if (pass_es)
	pass(" Element Sets ") << "Got element sets." << std::endl;
    else
	fail(" Element Sets ") << "Error in element sets." << std::endl;

    return pass_es;
}

bool TestHexMeshReader::check_element_types(
    const rtt_meshReaders::Hex_Mesh_Reader &mesh, const std::string &testid)
{
    // Check Element Types.
    typedef rtt_meshReaders::Element_Definition et;
    bool pass_et = true;
    std::vector<rtt_meshReaders::Element_Definition::Element_Type>
	etypes = mesh.get_element_types();
    if (testid == "slab")
    { 
	for (int i=0; i<100; ++i)
	    pass_et = pass_et && etypes[i] == et::BAR_2;
	for (int i=100; i<102; ++i)
	    etypes[i] == et::NODE;
    }
    else if (testid == "quad")
    {
	for (int i=0; i<400; ++i)
	    pass_et = pass_et && etypes[i] == et::QUAD_4;
	for (int i=400; i<480; ++i)
	    etypes[i] == et::BAR_2;
    }
    else if (testid == "cube") 
    {
	for (int i=0; i<125; ++i)
	    pass_et = pass_et && etypes[i] == et::HEXA_8;
	for (int i=125; i<275; ++i)
	    etypes[i] == et::QUAD_4;
    }
    else
	Insist(false,"Unrecognized test id string!");
    if (pass_et) 
	pass(" Element Types ") << "Read Element Types." << std::endl;
    else
	fail(" Element Types ") << "Error in Element Types." << std::endl;
    return pass_et;
}

bool TestHexMeshReader::compare_double(const double &lhs, const double &rhs)
{
    // Note that this is only good for doubles close to one.
    return std::fabs(lhs-rhs) <= 0.00001;
}

bool TestHexMeshReader::check_map(const std::map<std::string, std::set<int> >
			      &elmsets, const std::string &name, 
			      const int &begin, const int &end) 
{
    bool pass = true;
    typedef std::map<std::string, std::set<int> > mt;
    mt::const_iterator iter = elmsets.find(name);
    if (iter != elmsets.end())
    {
	const std::set<int> &elem_subset = (*iter).second;
	if (elem_subset.size() == end-begin) 
	{
	    for (int i=begin; i<end; ++i)
	    {
		std::set<int>::const_iterator siter = elem_subset.find(i);
		pass = pass && (siter != elem_subset.end());
	    }
	}
	else
	    pass = false;
    }
    else
	pass = false;
    return pass;
}

} // end namespace rtt_meshReaders_test


//---------------------------------------------------------------------------//
//                              end of TestHexMeshReader.cc
//---------------------------------------------------------------------------//
