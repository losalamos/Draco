//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders_Services/test/TestmeshReadersServices.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 27 14:52:05 2002
 * \brief  Test meshReaders connectivity services.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "meshReaders_Services_test.hh"
#include "TestmeshReadersServices.hh"
#include "../Release.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

using rtt_RTT_Format_Reader::RTT_Mesh_Reader;
using rtt_meshReaders_Services::Connect;

using rtt_meshReaders_Services_test::bndry_flags;
using rtt_meshReaders_Services_test::flag_numbs;
using rtt_meshReaders_Services_test::mesh;
using rtt_meshReaders_Services_test::connect;
using rtt_meshReaders_Services_test::compareXYZ;

typedef vector<int> vector_int;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void runTest()
{
    // New meshes added to this test will have to be added to the enumeration
    // Meshes in the header file.
    const int MAX_MESHES = 2;
    std::string filename[MAX_MESHES] = {"rttdef.mesh", "rttamr.mesh"};
    rtt_meshReaders_Services_test::Meshes mesh_type;

    for (int mesh_number = 0; mesh_number < MAX_MESHES; mesh_number++)
    {
	bndry_flags.resize(0);
	flag_numbs.resize(0);
	bool all_passed = true;
        switch (mesh_number)
	{
	    // Test all nested class accessor functions for a very simplistic 
	    // mesh file (enum DEFINED).
	case (0):
	    mesh_type = rtt_meshReaders_Services_test::DEFINED;
	    bndry_flags.push_back("boundary/reflective");
	    flag_numbs.push_back(1);
	    bndry_flags.push_back("boundary/vacuum");
	    flag_numbs.push_back(2);
	    break;

	    // Test a fairly simple AMR mesh (enum AMR). Just read it in and
	    // verify some data. That's enough!
	case (1):
	    mesh_type = rtt_meshReaders_Services_test::AMR;
	    bndry_flags.push_back("BNDRY/VACUU"); 
	    flag_numbs.push_back(1);
	    bndry_flags.push_back("BNDRY/REFLE"); 
	    flag_numbs.push_back(2);
	    break;

	default:
	    FAILMSG("Invalid mesh type encountered.");
	    all_passed = false;
	    break;
	}

        // Construct an Reader class object from the data in the 
	// specified mesh file. 
        mesh = new RTT_Mesh_Reader(filename[mesh_number]);
	{
	    ostringstream m;
	    m << "Read " << filename[mesh_number] 
	      << " without coreing in or firing an"
	      << " assertion." << std::endl;
	    PASSMSG(m.str());
	}

	connect = new Connect(mesh, bndry_flags, flag_numbs,
			      compareXYZ<std::vector<double> >());
	{
	    ostringstream m;
	    m << "Connected mesh without coreing in or"
	      << " firing an assertion." << std::endl;
	    PASSMSG(m.str());
	}

	all_passed = all_passed && check_connect(mesh, connect, mesh_type);

	if (!all_passed)
	{
	    ostringstream m;
	    m << "Errors occured testing mesh " 
	      << "number " << mesh_type << std::endl;
	    FAILMSG(m.str());
	}
    }

    // Report results of test.
    if (rtt_meshReaders_Services_test::passed)
    {
	PASSMSG("All tests passed.");
    }
    else
    {
	FAILMSG("Some tests failed.");
    }
}

//---------------------------------------------------------------------------//

namespace rtt_meshReaders_Services_test
{

std::vector<std::string>      bndry_flags;
std::vector<int>              flag_numbs;
rtt_dsxx::SP<RTT_Mesh_Reader> mesh;
rtt_dsxx::SP<Connect>         connect;

bool check_connect(rtt_dsxx::SP<RTT_Mesh_Reader> mesh, 
		   rtt_dsxx::SP<Connect> connect, 
		   const Meshes & meshtype)
{
    // Exercise the connectivity accessor functions for this mesh.
    bool all_passed = true;
    std::vector<std::vector<std::vector<int> > > adjacent_cell;
    std::multimap<int, int> bndryFaces;
    std::multimap<int, vector_int > Side_to_Cell_Face;
    // Rest of these variables are for the user to specify the tests to be
    // performed (i.e., the user has to initialize them in a case).
    vector_int cells_to_check; // Probably don't want to check every cell.
    vector_int bndryFacesCount;
    vector_int cell_face(2,0);
    
    switch (meshtype)
    {
    case DEFINED:
        cells_to_check.push_back(0); // there is only one cell in the mesh.
	adjacent_cell.resize(cells_to_check.size());
	// Resize the adjacent cell vector for the number of sides per cell.
	for (int c = 0; c < cells_to_check.size(); c++)
	    adjacent_cell[c].resize(4);

	// All of the "adjacent cells" are boundaries for this case, and the
	// adjacent cell is returned as the negative of the boundary flag #. 
	adjacent_cell[0][0].push_back(-1);
	adjacent_cell[0][1].push_back(-1); adjacent_cell[0][2].push_back(-2);
	adjacent_cell[0][3].push_back(-1);
	// bndryFaces pairs are keyed on face number (i.e., 0 - 5 for a hex) to
	// the cell number as a value and contain all the faces that are either
	// on the outer geometric boundary of the problem or a junction between
	// cells with differing refinement levels in an AMR mesh. Only the
	// subset of the multimap that is entered here will be checked. 
	bndryFaces.insert(std::make_pair(0,0)); 
	bndryFaces.insert(std::make_pair(1,0));
	bndryFaces.insert(std::make_pair(2,0)); 
	bndryFaces.insert(std::make_pair(3,0));
	// bndryFacesCount is input to check the number of cells that have 
	// faces that occur on the boundaries. The int values to be entered
	// here correspond to the entire mesh. 
	bndryFacesCount.resize(4); // Four faces for a tet.
	bndryFacesCount[0] = bndryFacesCount[1] = bndryFacesCount[2] =
	    bndryFacesCount[3] = 1;
	// Side_to_Cell_Face correlates the sides data in the RTT_Format file
	// with the equivalent cell face. The side is the key and the cell and
	// face number form a vector that is the data value. Only the  subset
	// of the multimap that is entered here will be checked.
	cell_face[1] = 2; 
	Side_to_Cell_Face.insert(std::make_pair(0,cell_face));
	cell_face[1] = 3; 
	Side_to_Cell_Face.insert(std::make_pair(1,cell_face));
	cell_face[1] = 1; 
	Side_to_Cell_Face.insert(std::make_pair(2,cell_face));
	cell_face[1] = 0; 
	Side_to_Cell_Face.insert(std::make_pair(3,cell_face));

	break;

    case AMR:
        cells_to_check.push_back(0); cells_to_check.push_back(9);
	cells_to_check.push_back(19); cells_to_check.push_back(28);
	cells_to_check.push_back(29); cells_to_check.push_back(39);
	cells_to_check.push_back(49); cells_to_check.push_back(56);
	adjacent_cell.resize(cells_to_check.size());
	// Resize the adjacent cell vector for the number of sides per cell.
	for (int c = 0; c < cells_to_check.size(); c++)
	    adjacent_cell[c].resize(6);

	// Specify the "adjacent cells". The adjacent cell is returned as the 
	// negative of the boundary flag # for cell on the outer boundary.
	// Cell 0.
	adjacent_cell[0][0].push_back(-2); adjacent_cell[0][1].push_back(-2);
	adjacent_cell[0][2].push_back(1); adjacent_cell[0][3].push_back(4); 
	adjacent_cell[0][4].push_back(-2); adjacent_cell[0][5].push_back(16);
	// Cell 9.
	adjacent_cell[1][0].push_back(-2); adjacent_cell[1][1].push_back(8);
	adjacent_cell[1][2].push_back(11); adjacent_cell[1][3].push_back(14); 
	adjacent_cell[1][4].push_back(7); adjacent_cell[1][5].push_back(28);
	// Cell 19.
	adjacent_cell[2][0].push_back(18); adjacent_cell[2][1].push_back(-2);
	adjacent_cell[2][2].push_back(21); adjacent_cell[2][3].push_back(28);
 	adjacent_cell[2][4].push_back(17); adjacent_cell[2][5].push_back(42);
	// Cell 28. This cell is the junction between generations on every 
	// face (i.e, four adjacent cells per cell face).
	adjacent_cell[3][0].push_back(6); adjacent_cell[3][0].push_back(8); 
	adjacent_cell[3][0].push_back(7); adjacent_cell[3][0].push_back(9);
	adjacent_cell[3][1].push_back(18); adjacent_cell[3][1].push_back(20);
	adjacent_cell[3][1].push_back(19); adjacent_cell[3][1].push_back(21);
	adjacent_cell[3][2].push_back(29); adjacent_cell[3][2].push_back(31);
	adjacent_cell[3][2].push_back(30); adjacent_cell[3][2].push_back(32);
	adjacent_cell[3][3].push_back(35); adjacent_cell[3][3].push_back(37);
	adjacent_cell[3][3].push_back(36); adjacent_cell[3][3].push_back(38);
 	adjacent_cell[3][4].push_back(24); adjacent_cell[3][4].push_back(26);
	adjacent_cell[3][4].push_back(25); adjacent_cell[3][4].push_back(27);
	adjacent_cell[3][5].push_back(47); adjacent_cell[3][5].push_back(49);
	adjacent_cell[3][5].push_back(48); adjacent_cell[3][5].push_back(50);
	// Cell 29.
	adjacent_cell[4][0].push_back(10); adjacent_cell[4][1].push_back(22);
	adjacent_cell[4][2].push_back(-2); adjacent_cell[4][3].push_back(31);
 	adjacent_cell[4][4].push_back(28); adjacent_cell[4][5].push_back(30);
	// Cell 39.
	adjacent_cell[5][0].push_back(15); adjacent_cell[5][1].push_back(31);
	adjacent_cell[5][2].push_back(-2); adjacent_cell[5][3].push_back(-2);
 	adjacent_cell[5][4].push_back(37); adjacent_cell[5][5].push_back(40);
	// Cell 49.
	adjacent_cell[6][0].push_back(28); adjacent_cell[6][1].push_back(43);
	adjacent_cell[6][2].push_back(51); adjacent_cell[6][3].push_back(50);
 	adjacent_cell[6][4].push_back(47); adjacent_cell[6][5].push_back(-2);
	// Cell 56.
	adjacent_cell[7][0].push_back(40); adjacent_cell[7][1].push_back(52);
	adjacent_cell[7][2].push_back(-2); adjacent_cell[7][3].push_back(-2);
 	adjacent_cell[7][4].push_back(55); adjacent_cell[7][5].push_back(-2);
	// bndryFaces pairs are keyed on face number (i.e., 0 - 5 for a hex) to
	// the cell number as a value and contain all the faces that are either
	// on the outer geometric boundary of the problem or a junction between
	// cells with differing refinement levels in an AMR mesh. Only the
	// subset of the multimap that is entered here will be checked. 
	bndryFaces.insert(std::make_pair(0,0)); 
	bndryFaces.insert(std::make_pair(1,0));
	bndryFaces.insert(std::make_pair(4,0));
	bndryFaces.insert(std::make_pair(0,9)); 
	bndryFaces.insert(std::make_pair(1,19)); 
	bndryFaces.insert(std::make_pair(0,28)); 
	bndryFaces.insert(std::make_pair(1,28));
	bndryFaces.insert(std::make_pair(2,28)); 
	bndryFaces.insert(std::make_pair(3,28));
	bndryFaces.insert(std::make_pair(4,28)); 
	bndryFaces.insert(std::make_pair(5,28));
	bndryFaces.insert(std::make_pair(2,29)); 
	bndryFaces.insert(std::make_pair(2,39));
	bndryFaces.insert(std::make_pair(3,39));
	bndryFaces.insert(std::make_pair(5,49));
	bndryFaces.insert(std::make_pair(2,56));
	bndryFaces.insert(std::make_pair(3,56));
	bndryFaces.insert(std::make_pair(5,56));
	// bndryFacesCount is input to check the number of cells that have 
	// faces that occur on the boundaries. The int values to be entered
	// here correspond to the entire mesh. 
	bndryFacesCount.resize(6); // Six faces for a Hexahedron.
	bndryFacesCount[0] = bndryFacesCount[1] = bndryFacesCount[2] =
	    bndryFacesCount[3] = bndryFacesCount[4] = bndryFacesCount[5] = 21;
	// Side_to_Cell_Face correlates the sides data in the RTT_Format file
	// with the equivalent cell face. The side is the key and the cell and
	// face number form a vector that is the data value. Only the  subset
	// of the multimap that is entered here will be checked.
	cell_face[0] = 0; cell_face[1] = 0; 
	Side_to_Cell_Face.insert(std::make_pair(93,cell_face));
	cell_face[1] = 1;
	Side_to_Cell_Face.insert(std::make_pair(94,cell_face));
	cell_face[1] = 4; 
	Side_to_Cell_Face.insert(std::make_pair(95,cell_face));
	cell_face[0] = 9; cell_face[1] = 0; 
	Side_to_Cell_Face.insert(std::make_pair(78,cell_face));
	cell_face[0] = 19; cell_face[1] = 1; 
	Side_to_Cell_Face.insert(std::make_pair(58,cell_face));
	cell_face[0] = 29; cell_face[1] = 2; 
	Side_to_Cell_Face.insert(std::make_pair(47,cell_face));
	cell_face[0] = 39; cell_face[1] = 2; 
	Side_to_Cell_Face.insert(std::make_pair(33,cell_face));
	cell_face[1] = 3; 
	Side_to_Cell_Face.insert(std::make_pair(35,cell_face));
	cell_face[0] = 49; cell_face[1] = 5; 
	Side_to_Cell_Face.insert(std::make_pair(15,cell_face));
	cell_face[0] = 56; cell_face[1] = 2; 
	Side_to_Cell_Face.insert(std::make_pair(0,cell_face));
	cell_face[1] = 3;
	Side_to_Cell_Face.insert(std::make_pair(1,cell_face));
	cell_face[1] = 5;
	Side_to_Cell_Face.insert(std::make_pair(2,cell_face));
	
	break;

    default:
	FAILMSG("Invalid mesh type encountered.");
	all_passed = false;
	return all_passed;
    }
    // Check the number of cells that are adjacent to the cells that are
    // specified in cells_to_check.
    bool got_adj_cell_size = true;
    for (int c = 0; c < cells_to_check.size(); c++)
    {
        int cell = cells_to_check[c];
	for (int f = 0; f < adjacent_cell[c].size(); f++)
	    if (adjacent_cell[c][f].size() != 
		connect->get_adjCell_size(cell,f))
	        got_adj_cell_size = false;
    }
    if (!got_adj_cell_size)
    {
        FAILMSG("Connectivity adjacent cell size not obtained.");
 	all_passed = false;
    }

    // Check the cell faces that are adjacent to the cells that are
    // specified in cells_to_check.
    bool got_adj_cell = true;
    for (int c = 0; c < cells_to_check.size(); c++)
    {
        int cell = cells_to_check[c];
	for (int f = 0; f < adjacent_cell[c].size(); f++)
	    for (int a = 0; a < adjacent_cell[c][f].size(); a++) 
	        if (adjacent_cell[c][f][a] != connect->get_adjCell(cell,f,a))
		    got_adj_cell = false;
    }
    if (!got_adj_cell)
    {
        FAILMSG("Connectivity adjacent cell not obtained.");
 	all_passed = false;
    }

    // Check the count of boundary cell faces.
    bool got_bndry_face_count = true;
    for (int f = 0; f < bndryFacesCount.size(); f++)
        if (bndryFacesCount[f] != connect->get_bndryFaces_count(f))
	    got_bndry_face_count = false;
    if (!got_bndry_face_count)
    {
        FAILMSG("Connectivity boundary face count not obtained.");
 	all_passed = false;
    }

    // Check the get_bndryCells accessor function for those face/cell pairs
    // that were specified for this case (in bndryFaces for the switch).
    bool got_bndry_cells = true;
    for (std::multimap<int, int>::iterator f = bndryFaces.begin();
	 f != bndryFaces.end(); f++)
	if (connect->get_bndryCells(f->first).count(f->second) == 0)
	    got_bndry_cells = false;
    if (!got_bndry_cells)
    {
        FAILMSG("Connectivity boundary cells not obtained.");
 	all_passed = false;
    }

    // Check the get_bndryFace accessor function for those face/cell pairs
    // that were specified for this case (in bndryFaces for the switch).
    bool got_bndry_face = true;
    for (std::multimap<int, int>::iterator f = bndryFaces.begin();
	 f != bndryFaces.end(); f++)
	if (!connect->check_bndryFace(f->second, f->first))
	    got_bndry_face = false;
    if (!got_bndry_face)
    {
        FAILMSG("Connectivity boundary face not obtained.");
 	all_passed = false;
    }

    // Check the get_Cell_from_Side accessor function for those sides
    // that were specified for this case (Side_to_Cell_Face in the switch).
    bool got_cell_from_side = true;
    for (std::multimap<int, vector_int >::iterator s = 
	     Side_to_Cell_Face.begin();
	 s != Side_to_Cell_Face.end(); s++)
	if (connect->get_Cell_from_Side(s->first) != s->second[0])
	    got_cell_from_side = false;
    if (!got_cell_from_side)
    {
        FAILMSG("Connectivity cell not obtained from side.");
 	all_passed = false;
    }

    // Check the get_Cell_Face_from_Side accessor function for those sides
    // that were specified for this case (Side_to_Cell_Face in the switch).
    bool got_cell_face_from_side = true;
    for (std::multimap<int, vector_int >::iterator s = 
	     Side_to_Cell_Face.begin();
	 s != Side_to_Cell_Face.end(); s++)
	if (connect->get_Cell_Face_from_Side(s->first) != s->second[1])
	    got_cell_face_from_side = false;
    if (!got_cell_face_from_side)
    {
        FAILMSG("Connectivity cell face not obtained from side.");
 	all_passed = false;
    }

    if (all_passed)
    {
        PASSMSG("Got all connectivity accessors.");
    }
    else
    {
	FAILMSG("Errors in some connectivity accessors.");
    }

    return all_passed;
}

} // end of rtt_meshReaders_Services_test

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_meshReaders_Services::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	runTest();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing TestmeshReadersServices, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_meshReaders_Services_test::passed) 
    {
        cout << "**** TestmeshReadersServices Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing TestmeshReadersServices." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of TestmeshReadersServices.cc
//---------------------------------------------------------------------------//
