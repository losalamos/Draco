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
#include "../Element_Definition.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>
#include <vector>

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
using std::multimap;
using std::make_pair;
using std::set;
using std::pair;
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
    // New meshes added to this test will have to be added to the enumeration
    // Meshes in the header file.
    const int MAX_MESHES = 2;
    string filename[MAX_MESHES] = {"rttdef.mesh", "rttdef.mesh"};
    Meshes mesh_type;
    bool renumber;

    for (int mesh_number = 0; mesh_number < MAX_MESHES; mesh_number++)
    {
        // Set the option to renumber the RTT_Format nodes, sides, and cells
        // in ascending order for those mesh-types that require this option 
        // (e.g., AMR mesh). Default is no renumbering.
        switch (mesh_number)
	{
	case(1):
	    renumber = true;
	    break;

	default:
	    renumber = false;
	}
        // Construct an RTT_Format class object from the data in the specified
        //  mesh file. 
        RTT_Format mesh(filename[mesh_number], renumber);
	pass(" Construct ") << "Read and connected " << filename[mesh_number] 
			    << " without coreing in or firing an assertion." 
			    << endl;
	bool all_passed = true;
	// The following switch allows addition of other meshes for testing,
	// with the "DEFINED" mesh providing an example. Only the check_dims
	// tests is required and it will be automatically called by the other
	// tests (with the exception of check header) if not invoked herein.
	// The comparison data must also be provided for additional meshes
	// within the switch structure residing in the test functions. 
        switch (mesh_number)
	{
	// Test all nested class accessor functions for a very simplistic 
	// mesh file (enum DEFINED).
	case (0):
	    mesh_type = DEFINED;
	    all_passed = all_passed && check_header(mesh, mesh_type);
	    all_passed = all_passed && check_dims(mesh, mesh_type);
	    all_passed = all_passed && check_node_flags(mesh, mesh_type);
	    all_passed = all_passed && check_side_flags(mesh, mesh_type);
	    all_passed = all_passed && check_cell_flags(mesh, mesh_type);
	    all_passed = all_passed && check_node_data_ids(mesh, mesh_type);
	    all_passed = all_passed && check_side_data_ids(mesh, mesh_type);
	    all_passed = all_passed && check_cell_data_ids(mesh, mesh_type);
	    all_passed = all_passed && check_cell_defs(mesh, mesh_type);
	    all_passed = all_passed && check_nodes(mesh, mesh_type);
	    all_passed = all_passed && check_sides(mesh, mesh_type);
	    all_passed = all_passed && check_cells(mesh, mesh_type);
	    all_passed = all_passed && check_node_data(mesh, mesh_type);
	    all_passed = all_passed && check_side_data(mesh, mesh_type);
	    all_passed = all_passed && check_cell_data(mesh, mesh_type);
	    all_passed = all_passed && check_virtual(mesh, mesh_type);
	    all_passed = all_passed && check_connectivity(mesh, mesh_type);
	    break;

	// Test all nested class accessor functions for the same simplistic 
	// mesh file with node, side, and cell renumbering implemented 
	// (enum SORTED). Note that nothing should change from the previous
	// case until the check_cell_defs accessor function call, but the 
	// preceding accessor function calls are included to verify that this
	// fact remains true.
	case (1):
	    mesh_type = SORTED;
	    all_passed = all_passed && check_header(mesh, mesh_type);
	    all_passed = all_passed && check_dims(mesh, mesh_type);
	    all_passed = all_passed && check_node_flags(mesh, mesh_type);
	    all_passed = all_passed && check_side_flags(mesh, mesh_type);
	    all_passed = all_passed && check_cell_flags(mesh, mesh_type);
	    all_passed = all_passed && check_node_data_ids(mesh, mesh_type);
	    all_passed = all_passed && check_side_data_ids(mesh, mesh_type);
	    all_passed = all_passed && check_cell_data_ids(mesh, mesh_type);
	    all_passed = all_passed && check_cell_defs(mesh, mesh_type);
	    all_passed = all_passed && check_renumber(mesh, mesh_type);
	    all_passed = all_passed && check_nodes(mesh, mesh_type);
	    all_passed = all_passed && check_sides(mesh, mesh_type);
	    all_passed = all_passed && check_cells(mesh, mesh_type);
	    all_passed = all_passed && check_node_data(mesh, mesh_type);
	    all_passed = all_passed && check_side_data(mesh, mesh_type);
	    all_passed = all_passed && check_cell_data(mesh, mesh_type);
	    break;

	default:
	    fail("runTest") << "Invalid mesh type encountered." << endl;
	    all_passed = false;
	    break;
	}
	if (!all_passed)
	    fail(filename[mesh_number]) << "Errors occured testing mesh " 
					<< "number " << mesh_type << endl;
    }

    // Report results of test.
    if (passed())
    {
	return "All tests passed.";
    }
    return "Some tests failed.";
}

bool TestRTT_Format::check_header(const rtt_meshReaders::RTT_Format & mesh,
				  const Meshes & meshtype)
{
    // Exercise the header accessor functions for this mesh.
    bool all_passed = true;
    string version, title, date;
    int cycle, ncomments;
    vector<string> comments;
    double time;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        version = "v1.0.0";
	title = "RTT_format mesh file definition, version 7.";
	date = "24 Jul 97";
	cycle = 1;
	time = 0.0;
	ncomments = 3;
	comments.push_back("One tet mesh in an RTT mesh file format.");
	comments.push_back("Date     : 24 Jul 97");
	comments.push_back("Author(s): H. Trease, J.McGhee");
	break;

    default:
        fail("check_header") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }
    // Check version.
    if (mesh.get_header_version() != version)
    {
        fail(" Header Version ") << "Header version not obtained." << endl;
 	all_passed = false;
    }
    // Check title.
    if (mesh.get_header_title() != title)
    {
        fail(" Header Title ") << "Header title not obtained." << endl;
 	all_passed = false;
    }
    // Check date.
    if (mesh.get_header_date() != date)
    {
        fail(" Header Date ") << "Header date not obtained." << endl;
 	all_passed = false;
    }
    // Check cycle.
    if (mesh.get_header_cycle() != cycle)
    {
        fail(" Header Cycle ") << "Header cycle not obtained." << endl;
 	all_passed = false;
    }
    // Check time.
    if (mesh.get_header_time() != time)
    {
        fail(" Header Time ") << "Header time not obtained." << endl;
 	all_passed = false;
    }
    // Check ncomments.
    if (mesh.get_header_ncomments() != ncomments)
    {
        fail(" Header Ncomments ") << "Header ncomments not obtained." << endl;
 	all_passed = false;
    }
    // Check comments.
    bool got_comments = true;
    for (int i= 0; i < ncomments; i++)
        if (comments[i] != mesh.get_header_comments(i))
	    got_comments = false;
    if (!got_comments)
    {
        fail(" Header Comments ") << "Header comments not obtained." << endl;
 	all_passed = false;
    }
    // Check that all Header class accessors passed their tests.
    if (all_passed)
        pass(" Header Accessors " ) << "Got all Header accessors." << endl;
    else
	fail(" Header Accessors ") << "Errors in some Header accessors." 
				   << endl;

    return all_passed;
}

bool TestRTT_Format::check_dims(const rtt_meshReaders::RTT_Format & mesh,
				const Meshes & meshtype)
{
    // Exercise the dims accessor functions for this mesh.
    bool all_passed = true;
    string coor_units, prob_time_units;
    int ncell_defs;
    int nnodes_max;
    int nsides_max;
    int nnodes_side_max;
    int ndim;
    int ndim_topo;
    int nnodes;
    int nnode_flag_types;
    vector<int> nnode_flags;
    int nnode_data;
    int nsides;
    int nside_types;
    vector<int> side_types;
    int nside_flag_types;
    vector<int> nside_flags;
    int nside_data;
    int ncells;
    int ncell_types;
    vector<int> cell_types;
    int ncell_flag_types;
    vector<int> ncell_flags;
    int ncell_data;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
	coor_units = "cm";
	prob_time_units= "s";
        ncell_defs = 8;
	nnodes_max = 8;
	nsides_max = 6;
	nnodes_side_max = 4;
	ndim = 3;
	ndim_topo = 3;
	nnodes = 4;
	nnode_flag_types = 3;
	    nnode_flags.push_back(3); 
	    nnode_flags.push_back(2);
	    nnode_flags.push_back(2);
	nnode_data = 3;
	nsides = 4;
	nside_types = 1;
	    // All side types are decremented relative to the value in the
	    // input file for zero indexing.
	    side_types.push_back(2);
	nside_flag_types = 1;
	    nside_flags.push_back(2);
	nside_data =2;
	ncells = 1;
	ncell_types = 1;
	    // All cell types are decremented relative to the value in the
	    // input file for zero indexing.
	    cell_types.push_back(5);
	ncell_flag_types = 2;
	    ncell_flags.push_back(2);
	    ncell_flags.push_back(2);
	ncell_data = 1;
	break;

    default:
        fail("check_dims") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }
    // Check coordinate units.
    if (mesh.get_dims_coor_units() != coor_units)
    {
        fail(" Dims coor_units ") << "Dims coor_units not obtained." << endl;
 	all_passed = false;
    }
    // Check problem time units.
    if (mesh.get_dims_prob_time_units() != prob_time_units)
    {
        fail(" Dims prob_time_units ") << "Dims prob_time_units not obtained."
				       << endl;
 	all_passed = false;
    }
    // Check number of cell definitions.
    if (mesh.get_dims_ncell_defs()!= ncell_defs)
    {
        fail(" Dims ncell_defs ") << "Dims ncell_defs not obtained." << endl;
 	all_passed = false;
    }
    // Check maximum number of nodes for cells in the "cell_defs" block.
    if (mesh.get_dims_nnodes_max() != nnodes_max)
    {
        fail(" Dims nnodes_max ") << "Dims nnodes_max not obtained." << endl;
 	all_passed = false;
    }
    // Check maximum number of sides for cells in the "cell_defs" block.
    if (mesh.get_dims_nsides_max() != nsides_max)
    {
        fail(" Dims nsides_max ") << "Dims nsides_max not obtained." << endl;
 	all_passed = false;
    }
    // Check maximum number of nodes/side for cells in the "cell_defs" block.
    if (mesh.get_dims_nnodes_side_max() != nnodes_side_max)
    {
        fail(" Dims nnodes_side_max ") << "Dims nnodes_side_max not obtained."
				       << endl;
 	all_passed = false;
    }
    // Check number of spatial dimensions.
    if (mesh.get_dims_ndim() != ndim)
    {
        fail(" Dims ndim ") << "Dims ndim not obtained." << endl;
 	all_passed = false;
    }
    // Check number of topological dimensions.
    if (mesh.get_dims_ndim_topo() != ndim_topo)
    {
        fail(" Dims ndim_topo ") << "Dims ndim_topo not obtained." << endl;
 	all_passed = false;
    }
    // Check total number of nodes in the mesh.
    if (mesh.get_dims_nnodes() != nnodes)
    {
        fail(" Dims nnodes ") << "Dims nnodes not obtained." << endl;
 	all_passed = false;
    }
    // Check number of node flag types.
    if (mesh.get_dims_nnode_flag_types() != nnode_flag_types)
    {
        fail(" Dims nnode_flag_types ") << 
	     "Dims nnode_flag_types not obtained." << endl;
 	all_passed = false;
    }
    // Check number of flags/node flag type.
    bool got_nnode_flags = true;
    for (int f = 0; f < nnode_flag_types; f++)
        if (mesh.get_dims_nnode_flags(f) != nnode_flags[f])
	    got_nnode_flags = false;
    if (!got_nnode_flags)
    {
        fail(" Dims nnode_flags ") << "Dims nnode_flags not obtained." << endl;
 	all_passed = false;
    }
    // Check number of node data fields.
    if (mesh.get_dims_nnode_data() != nnode_data)
    {
        fail(" Dims nnode_data ") << "Dims nnode_data not obtained." << endl;
 	all_passed = false;
    }
    // Check number of sides in the mesh.
    if (mesh.get_dims_nsides() != nsides)
    {
        fail(" Dims nsides ") << "Dims nsides not obtained." << endl;
 	all_passed = false;
    }
    // Check number of side types actually present in "side" block.
    if (mesh.get_dims_nside_types() != nside_types)
    {
        fail(" Dims nside_types ") << "Dims nside_types not obtained." << endl;
 	all_passed = false;
    }
    // Check side type indexes used in "side" block.
    bool got_side_types = true;
    for (int s = 0; s < nside_types; s++)
        if (mesh.get_dims_side_types(s) != side_types[s])
	    got_side_types = false;
    if (!got_side_types)
    {
        fail(" Dims side_types ") << "Dims side_types not obtained." << endl;
 	all_passed = false;
    }
    // Check number of side flag types.
    if (mesh.get_dims_nside_flag_types() != nside_flag_types)
    {
        fail(" Dims nside_flag_types ") << 
	     "Dims nside_flag_types not obtained." << endl;
 	all_passed = false;
    }
    // Check number of side flags/side flag type.
    bool got_nside_flags = true;
    for (int f = 0; f < nside_flag_types; f++)
        if (mesh.get_dims_nside_flags(f) != nside_flags[f])
	    got_nside_flags = false;
    if (!got_nside_flags)
    {
        fail(" Dims nside_flags ") << "Dims nside_flags not obtained." << endl;
 	all_passed = false;
    }
    // Check number of side data fields.
    if (mesh.get_dims_nside_data() != nside_data)
    {
        fail(" Dims nside_data ") << "Dims nside_data not obtained." << endl;
 	all_passed = false;
    }
    // Check total number of cells in the mesh.
    if (mesh.get_dims_ncells() != ncells)
    {
        fail(" Dims ncells ") << "Dims ncells not obtained." << endl;
 	all_passed = false;
    }
    // Check number of cell types actually present in "cells" block.
    if (mesh.get_dims_ncell_types() != ncell_types)
    {
        fail(" Dims ncell_types ") << "Dims ncell_types not obtained." << endl;
 	all_passed = false;
    }
    // Check cell type indexes used in "cells" block.
    bool got_ncell_types = true;
    for (int f = 0; f < ncell_types; f++)
        if (mesh.get_dims_cell_types(f) != cell_types[f])
	    got_ncell_types = false;
    if (!got_ncell_types) 
    {
        fail(" Dims cell_types ") << "Dims cell_types not obtained." << endl;
 	all_passed = false;
    }
    // Check number of cell flag types.
    if (mesh.get_dims_ncell_flag_types() != ncell_flag_types)
    {
        fail(" Dims ncell_flag_types ") << 
	    "Dims ncell_flag_types not obtained." << endl;
 	all_passed = false;
    }
    // Check number of flags/cell flag type.
    bool got_ncell_flags = true;
    for (int f = 0; f < ncell_flag_types; f++)
        if (mesh.get_dims_ncell_flags(f) != ncell_flags[f])
	    got_ncell_flags = false;
    if (!got_ncell_flags)
    {
        fail(" Dims ncell_flags ") << "Dims ncell_flags not obtained." << endl;
 	all_passed = false;
    }
    // Check number of cell data fields.
    if (mesh.get_dims_ncell_data() != ncell_data)
    {
        fail(" Dims ncell_data ") << "Dims ncell_data not obtained." << endl;
 	all_passed = false;
    }
    // Check that all Dims class accessors passed their tests.
    if (all_passed)
        pass(" Dims Accessors " ) << "Got all Dims accessors." << endl;
    else
	fail(" Dims Accessors ") << "Errors in some Dims accessors." 
				 << endl;

    // Retain the result of testing the Dims integrity for this mesh type.
    Dims_validated.insert(make_pair(meshtype, all_passed));

    return all_passed;
}
bool TestRTT_Format::check_node_flags(const rtt_meshReaders::RTT_Format & mesh,
				      const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the node_flags accessor functions for this mesh.
    bool all_passed = true;
    vector<string> flagTypes;
    vector<vector<pair<int, string> > > flag_num_name;
    vector<pair<int, string> > num_name;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        flagTypes.push_back("node_type");
	    num_name.push_back(make_pair(11,"interior"));
	    num_name.push_back(make_pair(21,"dudded"));
	    num_name.push_back(make_pair(6,"parent"));
	    flag_num_name.push_back(num_name);
	    num_name.resize(0);
	flagTypes.push_back("boundary");
	    num_name.push_back(make_pair(1,"reflective"));
	    num_name.push_back(make_pair(4,"vacuum"));
	    flag_num_name.push_back(num_name);
	    num_name.resize(0);
	flagTypes.push_back("source");
	    num_name.push_back(make_pair(101,"no_source"));
	    num_name.push_back(make_pair(22,"rad_source"));
	    flag_num_name.push_back(num_name);
	    num_name.resize(0);
	break;

    default:
        fail("check_node_flags") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check node flag types.
    bool got_node_flag_types = true;
    for (int i = 0; i < mesh.get_dims_nnode_flag_types(); i++)
        if (flagTypes[i] != mesh.get_node_flags_flag_type(i))
	    got_node_flag_types = false;
    if (!got_node_flag_types)
    {
        fail(" NodeFlags flag_type ") << "Node Flags flag_types not obtained."
				      << endl;
 	all_passed = false;
    }
    // Check node flag numbers for each of the flag types.
    bool got_node_flag_numbers = true;
    for (int i = 0; i < mesh.get_dims_nnode_flag_types(); i++)
    {
         num_name = flag_num_name[i];
         for (int j = 0; j < mesh.get_dims_nnode_flags(i); j++) 
	     if (num_name[j].first != 
		 mesh.get_node_flags_flag_number(i,j))
	         got_node_flag_numbers = false;
    }
    if (!got_node_flag_numbers)
    {
        fail(" NodeFlags flag_number ") << 
	     "Node Flags flag_numbers not obtained." << endl;
 	all_passed = false;
    }
    // Check number of flags for each node flag type.
    bool got_node_flag_size = true;
    for (int i = 0; i < mesh.get_dims_nnode_flag_types(); i++)
    {
        if (flag_num_name[i].size() != mesh.get_node_flags_flag_size(i))
	    got_node_flag_size = false;
    }
    if (!got_node_flag_size)
    {
        fail(" NodeFlags flag_size ") << "Node Flags flag_size not obtained."
				      << endl;
 	all_passed = false;
    }
    // Check node flag names for each of the flag types.
    bool got_node_flag_name = true;
    for (int i = 0; i < mesh.get_dims_nnode_flag_types(); i++)
    {
        num_name = flag_num_name[i];
        for (int j = 0; j < mesh.get_dims_nnode_flags(i); j++) 
	     if (num_name[j].second != 
		 mesh.get_node_flags_flag_name(i,j))
	         got_node_flag_name = false;
    }
    if (!got_node_flag_name)
    {
        fail(" NodeFlags flag_name ") << "Node Flags flag_name not obtained." 
				      << endl;
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
bool TestRTT_Format::check_side_flags(const rtt_meshReaders::RTT_Format & mesh,
				      const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the side_flags accessor functions for this mesh.
    bool all_passed = true;
    vector<string> flagTypes;
    vector<vector<pair<int, string> > > flag_num_name;
    vector<pair<int, string> > num_name;
    int bndry, src;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
	flagTypes.push_back("boundary");
	    num_name.push_back(make_pair(1,"reflective"));
	    num_name.push_back(make_pair(2,"vacuum"));
	    flag_num_name.push_back(num_name);
	    num_name.resize(0);
	    bndry = 0;
	    src = -1;
	break;

    default:
        fail("check_side_flags") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check side flag types.
    bool got_side_flag_types = true;
    for (int i = 0; i < mesh.get_dims_nside_flag_types(); i++)
        if (flagTypes[i] != mesh.get_side_flags_flag_type(i))
	    got_side_flag_types = false;
    if (!got_side_flag_types)
    {
        fail(" SideFlags flag_type ") << "Side Flags flag_types not obtained."
				      << endl;
 	all_passed = false;
    }
    // Check side flag numbers for each of the flag types.
    bool got_side_flag_numbers = true;
    for (int i = 0; i < mesh.get_dims_nside_flag_types(); i++)
    {
         num_name = flag_num_name[i];
         for (int j = 0; j < mesh.get_dims_nside_flags(i); j++) 
	     if (num_name[j].first != 
		 mesh.get_side_flags_flag_number(i,j))
	         got_side_flag_numbers = false;
    }
    if (!got_side_flag_numbers)
    {
        fail(" SideFlags flag_number ") << 
	     "Side Flags flag_numbers not obtained." << endl;
 	all_passed = false;
    }
    // Check number of flags for each side flag type.
    bool got_side_flag_size = true;
    for (int i = 0; i < mesh.get_dims_nside_flag_types(); i++)
    {
        if (flag_num_name[i].size() != mesh.get_side_flags_flag_size(i))
	    got_side_flag_size = false;
    }
    if (!got_side_flag_size)
    {
        fail(" SideFlags flag_size ") << "Side Flags flag_size not obtained."
				      << endl;
 	all_passed = false;
    }
    // Check side flag names for each of the flag types.
    bool got_side_flag_name = true;
    for (int i = 0; i < mesh.get_dims_nside_flag_types(); i++)
    {
        num_name = flag_num_name[i];
        for (int j = 0; j < mesh.get_dims_nside_flags(i); j++) 
	     if (num_name[j].second != 
		 mesh.get_side_flags_flag_name(i,j))
	         got_side_flag_name = false;
    }
    if (!got_side_flag_name)
    {
        fail(" SideFlags flag_name ") << "Side Flags flag_name not obtained."
				      << endl;
 	all_passed = false;
    }
    // Check side flag required boundary flag number.
    if (bndry != mesh.get_side_flags_boundary_flag_number())
    {
        fail(" SideFlags bndry_flag ") << 
	     "Side Flags boundary flag not obtained." << endl;
 	all_passed = false;
    }
    // Check side flag optional surface source flag number.
    if (src != mesh.get_side_flags_surface_src_flag_number())
    {
        fail(" SideFlags src_flag ") << 
	     "Side Flags surface source flag not obtained." << endl;
 	all_passed = false;
    }
    if (all_passed)
        pass(" SideFlags Accessors " ) << "Got all SideFlags accessors." 
				       << endl;
    else
	fail(" SideFlags Accessors ") << "Errors in some SideFlags accessors." 
				      << endl;

    return all_passed;
}
bool TestRTT_Format::check_cell_flags(const rtt_meshReaders::RTT_Format & mesh,
				      const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the cell_flags accessor functions for this mesh.
    bool all_passed = true;
    vector<string> flagTypes;
    vector<vector<pair<int, string> > > flag_num_name;
    vector<pair<int, string> > num_name;
    int matl, vsrc, rsrc;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        flagTypes.push_back("material");
	    num_name.push_back(make_pair(1,"control_rod"));
	    num_name.push_back(make_pair(2,"shield"));
	    flag_num_name.push_back(num_name);
	    num_name.resize(0);
	    matl = 0;
	flagTypes.push_back("rad_source");
	    num_name.push_back(make_pair(1,"src_name1"));
	    num_name.push_back(make_pair(2,"src_name2"));
	    flag_num_name.push_back(num_name);
	    num_name.resize(0);
	    rsrc = 1;
	vsrc = -1;
	break;

    default:
        fail("check_cell_flags") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check cell flag types.
    bool got_cell_flag_types = true;
    for (int i = 0; i < mesh.get_dims_ncell_flag_types(); i++)
        if (flagTypes[i] != mesh.get_cell_flags_flag_type(i))
	    got_cell_flag_types = false;
    if (!got_cell_flag_types)
    {
        fail(" CellFlags flag_type ") << "Cell Flags flag_types not obtained."
				      << endl;
 	all_passed = false;
    }
    // Check cell flag numbers for each of the flag types.
    bool got_cell_flag_numbers = true;
    for (int i = 0; i < mesh.get_dims_ncell_flag_types(); i++)
    {
         num_name = flag_num_name[i];
         for (int j = 0; j < mesh.get_dims_ncell_flags(i); j++) 
	     if (num_name[j].first != 
		 mesh.get_cell_flags_flag_number(i,j))
	         got_cell_flag_numbers = false;
    }
   if (!got_cell_flag_numbers)
    {
        fail(" CellFlags flag_number ") << 
	     "Cell Flags flag_numbers not obtained." << endl;
 	all_passed = false;
    }
    // Check number of flags for each cell flag type.
    bool got_cell_flag_size = true;
    for (int i = 0; i < mesh.get_dims_ncell_flag_types(); i++)
    {
        if (flag_num_name[i].size() != mesh.get_cell_flags_flag_size(i))
	    got_cell_flag_size = false;
    }
    if (!got_cell_flag_size)
    {
        fail(" CellFlags flag_size ") << "Cell Flags flag_size not obtained."
				      << endl;
 	all_passed = false;
    }
    // Check cell flag names for each of the flag types.
    bool got_cell_flag_name = true;
    for (int i = 0; i < mesh.get_dims_ncell_flag_types(); i++)
    {
        num_name = flag_num_name[i];
        for (int j = 0; j < mesh.get_dims_ncell_flags(i); j++) 
	     if (num_name[j].second != 
		 mesh.get_cell_flags_flag_name(i,j))
	         got_cell_flag_name = false;
    }
    if (!got_cell_flag_name)
    {
        fail(" CellFlags flag_name ") << "Cell Flags flag_name not obtained." 
				      << endl;
 	all_passed = false;
    }
    // Check cell flag required material flag number.
    if (matl != mesh.get_cell_flags_material_flag_number())
    {
        fail(" CellFlags matl_flag ") << 
	     "Cell Flags material flag not obtained." << endl;
 	all_passed = false;
    }
    // Check cell flag optional volume source flag number.
    if (vsrc != mesh.get_cell_flags_volume_src_flag_number())
    {
        fail(" CellFlags vsrc_flag ") << 
	     "Cell Flags volume source flag not obtained." << endl;
 	all_passed = false;
    }
     // Check cell flag optional radiation source flag number.
    if (rsrc != mesh.get_cell_flags_radiation_src_flag_number())
    {
        fail(" CellFlags rsrc_flag ") << 
	     "Cell Flags volume source flag not obtained." << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" CellFlags Accessors " ) << "Got all CellFlags accessors." 
				       << endl;
    else
	fail(" CellFlags Accessors ") << "Errors in some CellFlags accessors." 
				      << endl;

    return all_passed;
}

bool TestRTT_Format::check_node_data_ids(const rtt_meshReaders::RTT_Format & 
					 mesh, const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the node_data_ids accessor functions for this mesh.
    bool all_passed = true;
    vector<string> names;
    vector<string> units;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        names.push_back("density"); units.push_back("gm/cm**3");
        names.push_back("ion_temp"); units.push_back("keV");
        names.push_back("x_vel"); units.push_back("cm/sec");
	break;

    default:
        fail("check_node_data_ids") << "Invalid mesh type encountered." 
				    << endl;
	all_passed = false;
	return all_passed;
    }

    // Check node data id names.
    bool got_node_data_id_names = true;
    for (int i = 0; i < mesh.get_dims_nnode_data(); i++)
    {
        if (names[i] != mesh.get_node_data_id_name(i))
	    got_node_data_id_names = false;
    }
    if (!got_node_data_id_names)
    {
        fail(" NodeDataIDs name ") << "NodeDataIDs names not obtained." 
				   << endl;
 	all_passed = false;
    }

    // Check node data id units.
    bool got_node_data_id_units = true;
    for (int i = 0; i < mesh.get_dims_nnode_data(); i++)
    {
        if (units[i] != mesh.get_node_data_id_units(i))
	    got_node_data_id_units = false;
    }
    if (!got_node_data_id_units)
    {
        fail(" NodeDataIDs unit ") << "NodeDataIDs units not obtained." 
				   << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" NodeDataIDs Accessors " ) << "Got all NodeDataIDs accessors." 
					 << endl;
    else
	fail(" NodeDataIDs Accessors ") << 
             "Errors in some NodeDataIDs accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_side_data_ids(const rtt_meshReaders::RTT_Format & 
					 mesh, const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the side_data_ids accessor functions for this mesh.
    bool all_passed = true;
    vector<string> names;
    vector<string> units;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        names.push_back("density"); units.push_back("gm/cm**3");
        names.push_back("ion_temp"); units.push_back("keV");
	break;

    default:
        fail("check_side_data_ids") << "Invalid mesh type encountered." 
				    << endl;
	all_passed = false;
	return all_passed;
    }

    // Check side data id names.
    bool got_side_data_id_names = true;
    for (int i = 0; i < mesh.get_dims_nside_data(); i++)
    {
        if (names[i] != mesh.get_side_data_id_name(i))
	    got_side_data_id_names = false;
    }
    if (!got_side_data_id_names)
    {
        fail(" SideDataIDs name ") << "SideDataIDs names not obtained." 
				   << endl;
 	all_passed = false;
    }

    // Check side data id units.
    bool got_side_data_id_units = true;
    for (int i = 0; i < mesh.get_dims_nside_data(); i++)
    {
        if (units[i] != mesh.get_side_data_id_units(i))
	    got_side_data_id_units = false;
    }
    if (!got_side_data_id_units)
    {
        fail(" SideDataIDs unit ") << "SideDataIDs units not obtained." 
				   << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" SideDataIDs Accessors " ) << "Got all SideDataIDs accessors." 
					 << endl;
    else
	fail(" SideDataIDs Accessors ") << 
             "Errors in some SideDataIDs accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_cell_data_ids(const rtt_meshReaders::RTT_Format & 
					 mesh, const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the cell_data_ids functions for this mesh.
    bool all_passed = true;
    vector<string> names;
    vector<string> units;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        names.push_back("density"); units.push_back("gm/cm**3");
	break;

    default:
        fail("check_cell_data_ids") << "Invalid mesh type encountered." 
				    << endl;
	all_passed = false;
	return all_passed;
    }

    // Check cell data id names.
    bool got_cell_data_id_names = true;
    for (int i = 0; i < mesh.get_dims_ncell_data(); i++)
    {
        if (names[i] != mesh.get_cell_data_id_name(i))
	    got_cell_data_id_names = false;
    }
    if (!got_cell_data_id_names)
    {
        fail(" CellDataIDs name ") << "CellDataIDs names not obtained." 
				   << endl;
 	all_passed = false;
    }

    // Check cell data id units.
    bool got_cell_data_id_units = true;
    for (int i = 0; i < mesh.get_dims_ncell_data(); i++)
    {
        if (units[i] != mesh.get_cell_data_id_units(i))
	    got_cell_data_id_units = false;
    }
    if (!got_cell_data_id_units)
    {
        fail(" CellDataIDs unit ") << "CellDataIDs units not obtained." 
				   << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" CellDataIDs Accessors " ) << "Got all CellDataIDs accessors." 
					 << endl;
    else
	fail(" CellDataIDs Accessors ") << 
             "Errors in some CellDataIDs accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_cell_defs(const rtt_meshReaders::RTT_Format & mesh,
				     const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the cell_defs accessor functions for this mesh.
    bool all_passed = true;
    vector<string> names;
    vector<int> nnodes;
    vector<int> nsides;
    vector<vector<int> > side_types;
    vector<vector<set<int> > > sides;
    vector<vector<vector<int> > > ordered_sides;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        // Size vectors and load data consistent between both the standard
        //  and sorted cell definitions (names, nnodes, nsides).
        side_types.resize(8);
	sides.resize(8);
	ordered_sides.resize(8);
	// Size and load point data
        names.push_back("point");
	nnodes.push_back(1);
	nsides.push_back(0);
	side_types[0].resize(0);
	sides[0].resize(0);
	ordered_sides[0].resize(0);
	// Size and load line data.
        names.push_back("line");
	nnodes.push_back(2);
	nsides.push_back(2);
	side_types[1].resize(2,0);
	sides[1].resize(2);
	ordered_sides[1].resize(2);
	ordered_sides[1][0].push_back(0); 
	ordered_sides[1][1].push_back(1);
	// Size triangle data.
        names.push_back("triangle");
	nnodes.push_back(3);
	nsides.push_back(3);
	side_types[2].resize(3,1);
	sides[2].resize(3);
	ordered_sides[2].resize(3);
	// Size quad data.
        names.push_back("quad");
	nnodes.push_back(4);
	nsides.push_back(4);
        side_types[3].resize(4,1);
	sides[3].resize(4);
	ordered_sides[3].resize(4);
	// Size quad pyramid data.
        names.push_back("quad_pyr");
	nnodes.push_back(5);
	nsides.push_back(5);
	side_types[4].resize(5,2);
	side_types[4][0] = 3;
	sides[4].resize(5);
	ordered_sides[4].resize(5);
	// Size tetrahedron data.
        names.push_back("tetrahedron");
	nnodes.push_back(4);
	nsides.push_back(4);
	side_types[5].resize(4,2);
	sides[5].resize(4);
	ordered_sides[5].resize(4);
	// Size tri_prism data.
        names.push_back("tri_prism");
        nnodes.push_back(6);
        nsides.push_back(5);
        side_types[6].resize(5,3);
        sides[6].resize(5);
        ordered_sides[6].resize(5);
	// Size hexahedron data.
        names.push_back("hexahedron");
        nnodes.push_back(8);
        nsides.push_back(6);
        side_types[7].resize(6,3);
        sides[7].resize(6);
        ordered_sides[7].resize(6);
	switch (meshtype)
	{
	case DEFINED:
	    // triangle
	    ordered_sides[2][0].push_back(1); 
	        ordered_sides[2][0].push_back(2); 
	    ordered_sides[2][1].push_back(2); 
	        ordered_sides[2][1].push_back(0); 
	    ordered_sides[2][2].push_back(0); 
	        ordered_sides[2][2].push_back(1); 
	    // quad
	    ordered_sides[3][0].push_back(0); 
	        ordered_sides[3][0].push_back(1); 
	    ordered_sides[3][1].push_back(1); 
	        ordered_sides[3][1].push_back(2); 
	    ordered_sides[3][2].push_back(2); 
	        ordered_sides[3][2].push_back(3); 
 	    ordered_sides[3][3].push_back(3); 
	        ordered_sides[3][3].push_back(0); 
	    //quad_pyr
	    ordered_sides[4][0].push_back(0); 
	        ordered_sides[4][0].push_back(3); 
	        ordered_sides[4][0].push_back(2); 
	        ordered_sides[4][0].push_back(1); 
	    ordered_sides[4][1].push_back(0); 
	        ordered_sides[4][1].push_back(1); 
	        ordered_sides[4][1].push_back(4); 
	    ordered_sides[4][2].push_back(1); 
	        ordered_sides[4][2].push_back(2); 
	        ordered_sides[4][2].push_back(4); 
 	    ordered_sides[4][3].push_back(2); 
	        ordered_sides[4][3].push_back(3); 
	        ordered_sides[4][3].push_back(4); 
 	    ordered_sides[4][4].push_back(3); 
	        ordered_sides[4][4].push_back(0); 
	        ordered_sides[4][4].push_back(4); 
	    //tetrahedron
	    ordered_sides[5][0].push_back(1); 
	        ordered_sides[5][0].push_back(2); 
	        ordered_sides[5][0].push_back(3); 
	    ordered_sides[5][1].push_back(0); 
	        ordered_sides[5][1].push_back(3); 
	        ordered_sides[5][1].push_back(2); 
	    ordered_sides[5][2].push_back(0); 
	        ordered_sides[5][2].push_back(1); 
	        ordered_sides[5][2].push_back(3); 
 	    ordered_sides[5][3].push_back(0); 
	        ordered_sides[5][3].push_back(2); 
	        ordered_sides[5][3].push_back(1); 
	    // tri_prism
	    // the side_types vary between the cases for this cell def. 
	    side_types[6][0] = side_types[6][1] = 2;
	    ordered_sides[6][0].push_back(0); 
	        ordered_sides[6][0].push_back(2); 
	        ordered_sides[6][0].push_back(1); 
	    ordered_sides[6][1].push_back(3); 
	        ordered_sides[6][1].push_back(4); 
	        ordered_sides[6][1].push_back(5); 
	    ordered_sides[6][2].push_back(0); 
	        ordered_sides[6][2].push_back(1); 
	        ordered_sides[6][2].push_back(4); 
	        ordered_sides[6][2].push_back(3); 
 	    ordered_sides[6][3].push_back(0); 
	        ordered_sides[6][3].push_back(3); 
	        ordered_sides[6][3].push_back(5); 
	        ordered_sides[6][3].push_back(2); 
 	    ordered_sides[6][4].push_back(1); 
	        ordered_sides[6][4].push_back(2); 
	        ordered_sides[6][4].push_back(5); 
	        ordered_sides[6][4].push_back(4); 
	    // hexahedron
	    ordered_sides[7][0].push_back(0); 
	        ordered_sides[7][0].push_back(3); 
	        ordered_sides[7][0].push_back(2); 
	        ordered_sides[7][0].push_back(1); 
	    ordered_sides[7][1].push_back(4); 
	        ordered_sides[7][1].push_back(5); 
	        ordered_sides[7][1].push_back(6); 
	        ordered_sides[7][1].push_back(7); 
	    ordered_sides[7][2].push_back(0); 
	        ordered_sides[7][2].push_back(1); 
	        ordered_sides[7][2].push_back(5); 
	        ordered_sides[7][2].push_back(4); 
 	    ordered_sides[7][3].push_back(1); 
	        ordered_sides[7][3].push_back(2); 
	        ordered_sides[7][3].push_back(6); 
	        ordered_sides[7][3].push_back(5); 
 	    ordered_sides[7][4].push_back(2); 
	        ordered_sides[7][4].push_back(3); 
	        ordered_sides[7][4].push_back(7); 
	        ordered_sides[7][4].push_back(6); 
 	    ordered_sides[7][5].push_back(0); 
	        ordered_sides[7][5].push_back(4); 
	        ordered_sides[7][5].push_back(7); 
	        ordered_sides[7][5].push_back(3);
	    break;

	case SORTED:
	    // triangle
	    ordered_sides[2][0].push_back(1); 
	        ordered_sides[2][0].push_back(0); 
	    ordered_sides[2][1].push_back(0); 
	        ordered_sides[2][1].push_back(2); 
	    ordered_sides[2][2].push_back(2); 
	        ordered_sides[2][2].push_back(1); 
	    // quad
	    ordered_sides[3][0].push_back(1); 
	        ordered_sides[3][0].push_back(0); 
	    ordered_sides[3][1].push_back(0); 
	        ordered_sides[3][1].push_back(2); 
	    ordered_sides[3][2].push_back(3); 
	        ordered_sides[3][2].push_back(1); 
 	    ordered_sides[3][3].push_back(2); 
	        ordered_sides[3][3].push_back(3); 
	    //quad_pyr
	    ordered_sides[4][0].push_back(0); 
	        ordered_sides[4][0].push_back(1); 
	        ordered_sides[4][0].push_back(3); 
	        ordered_sides[4][0].push_back(2); 
	    ordered_sides[4][1].push_back(0); 
	        ordered_sides[4][1].push_back(4); 
	        ordered_sides[4][1].push_back(1); 
	    ordered_sides[4][2].push_back(0); 
	        ordered_sides[4][2].push_back(2); 
	        ordered_sides[4][2].push_back(4); 
 	    ordered_sides[4][3].push_back(1); 
	        ordered_sides[4][3].push_back(4); 
	        ordered_sides[4][3].push_back(3); 
 	    ordered_sides[4][4].push_back(2); 
	        ordered_sides[4][4].push_back(3); 
	        ordered_sides[4][4].push_back(4); 
	    //tetrahedron
	    ordered_sides[5][0].push_back(0); 
	        ordered_sides[5][0].push_back(1); 
	        ordered_sides[5][0].push_back(2); 
	    ordered_sides[5][1].push_back(0); 
	        ordered_sides[5][1].push_back(3); 
	        ordered_sides[5][1].push_back(1); 
	    ordered_sides[5][2].push_back(0); 
	        ordered_sides[5][2].push_back(2); 
	        ordered_sides[5][2].push_back(3); 
 	    ordered_sides[5][3].push_back(1); 
	        ordered_sides[5][3].push_back(3); 
	        ordered_sides[5][3].push_back(2); 
	    // tri_prism
	    // the side_types vary between the cases for this cell def. 
	    side_types[6][0] = side_types[6][4] = 2;
	    ordered_sides[6][0].push_back(0); 
	        ordered_sides[6][0].push_back(1); 
	        ordered_sides[6][0].push_back(2); 
	    ordered_sides[6][1].push_back(0); 
	        ordered_sides[6][1].push_back(3);
		ordered_sides[6][1].push_back(4);
	        ordered_sides[6][1].push_back(1); 
	    ordered_sides[6][2].push_back(0); 
	        ordered_sides[6][2].push_back(2); 
	        ordered_sides[6][2].push_back(5); 
	        ordered_sides[6][2].push_back(3); 
 	    ordered_sides[6][3].push_back(1); 
	        ordered_sides[6][3].push_back(4); 
	        ordered_sides[6][3].push_back(5); 
	        ordered_sides[6][3].push_back(2); 
 	    ordered_sides[6][4].push_back(3); 
	        ordered_sides[6][4].push_back(5); 
	        ordered_sides[6][4].push_back(4); 
	    // hexahedron
	    ordered_sides[7][0].push_back(0); 
	        ordered_sides[7][0].push_back(1); 
	        ordered_sides[7][0].push_back(3); 
	        ordered_sides[7][0].push_back(2); 
	    ordered_sides[7][1].push_back(0); 
	        ordered_sides[7][1].push_back(4); 
	        ordered_sides[7][1].push_back(5); 
	        ordered_sides[7][1].push_back(1); 
	    ordered_sides[7][2].push_back(0); 
	        ordered_sides[7][2].push_back(2); 
	        ordered_sides[7][2].push_back(6); 
	        ordered_sides[7][2].push_back(4); 
 	    ordered_sides[7][3].push_back(1); 
	        ordered_sides[7][3].push_back(5); 
	        ordered_sides[7][3].push_back(7); 
	        ordered_sides[7][3].push_back(3); 
 	    ordered_sides[7][4].push_back(2); 
	        ordered_sides[7][4].push_back(3); 
	        ordered_sides[7][4].push_back(7); 
	        ordered_sides[7][4].push_back(6); 
 	    ordered_sides[7][5].push_back(4); 
	        ordered_sides[7][5].push_back(6); 
	        ordered_sides[7][5].push_back(7); 
	        ordered_sides[7][5].push_back(5);
	    break;
	}
	// Load the (sorted) sides data from the ordered_sides data specified
	// above. The first element of ordered_sides and sides vectors 
	// corresponds to a point (no sides, and thus the second "array" 
	// size is zero) so start at cell definition one.
	for (int c = 1; c < ordered_sides.size(); c++)
	    for (int s = 0; s < ordered_sides[c].size(); s++)
	        for (int n = 0; n < ordered_sides[c][s].size(); n++)
		    sides[c][s].insert(ordered_sides[c][s][n]);
	break;

    default:
        fail("check_cell_defs") << "Invalid mesh type encountered." 
				    << endl;
	all_passed = false;
	return all_passed;
    }

    // Check cell definition names.
    bool got_cell_defs_names = true;
    for (int i = 0; i < mesh.get_dims_ncell_defs(); i++)
    {
        if (names[i] != mesh.get_cell_defs_name(i))
	    got_cell_defs_names = false;
    }
    if (!got_cell_defs_names)
    {
        fail(" CellDefs name ") << "CellDefs names not obtained." 
				   << endl;
 	all_passed = false;
    }

    // Check cell definition number of nodes.
    bool got_cell_defs_nnodes = true;
    for (int i = 0; i < mesh.get_dims_ncell_defs(); i++)
    {
        if (nnodes[i] != mesh.get_cell_defs_nnodes(i))
	    got_cell_defs_nnodes = false;
    }
    if (!got_cell_defs_nnodes)
    {
        fail(" CellDefs nnodes ") << "CellDefs nnodes not obtained." 
				   << endl;
 	all_passed = false;
    }
    // Check cell definition number of sides.
    bool got_cell_defs_nsides = true;
    for (int i = 0; i < mesh.get_dims_ncell_defs(); i++)
    {
        if (nsides[i] != mesh.get_cell_defs_nsides(i))
	    got_cell_defs_nsides = false;
    }
    if (!got_cell_defs_nsides)
    {
        fail(" CellDefs nsides ") << "CellDefs nsides not obtained." 
				   << endl;
 	all_passed = false;
    }
    // Check cell definition side types.
    bool got_cell_defs_side_types = true;
    for (int i = 0; i < mesh.get_dims_ncell_defs(); i++)
    {
        for (int s = 0; s < mesh.get_cell_defs_nsides(i); s++)
	    if (side_types[i][s] != mesh.get_cell_defs_side_types(i,s))
	        got_cell_defs_side_types = false;
    }
    if (!got_cell_defs_side_types)
    {
        fail(" CellDefs side_types ") << "CellDefs side_types not obtained." 
				   << endl;
 	all_passed = false;
    }

    // Check cell definition side sets.
    bool got_cell_defs_side = true;
    for (int i = 0; i < mesh.get_dims_ncell_defs(); i++)
    {
        for (int s = 0; s < mesh.get_cell_defs_nsides(i); s++)
	    if (sides[i][s] != mesh.get_cell_defs_side(i,s))
	        got_cell_defs_side = false;
    }
    if (!got_cell_defs_side)
    {
        fail(" CellDefs side ") << "CellDefs side not obtained." << endl;
 	all_passed = false;
    }
    // Check cell definition ordered_side sets.
    bool got_cell_defs_ordered_side = true;
    for (int i = 0; i < mesh.get_dims_ncell_defs(); i++)
    {
        for (int s = 0; s < mesh.get_cell_defs_nsides(i); s++)
	    if (ordered_sides[i][s] != mesh.get_cell_defs_ordered_side(i,s))
	        got_cell_defs_ordered_side = false;
    }
    if (!got_cell_defs_ordered_side)
    {
        fail(" CellDefs ordered_side ") << 
	     "CellDefs ordered_side not obtained." << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" CellDefs Accessors " ) << "Got all CellDefs accessors." << endl;
    else
	fail(" CellDefs Accessors ") << 
             "Errors in some CellDefs accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_renumber(const rtt_meshReaders::RTT_Format & mesh, 
				    const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;
    // Return if the renumbering flag is not set for the mesh constructor.
    if (!mesh.get_dims_renumber())
        return false;

    // Exercise the renumbering accessor functions for this mesh.
    bool all_passed = true;
    vector<int> node_map(mesh.get_dims_nnodes());
    vector<int> side_map(mesh.get_dims_nsides());
    vector<int> cell_map(mesh.get_dims_ncells());

    switch (meshtype)
    {
    case SORTED:
        node_map[0] = 0; node_map[1] = 3; node_map[2] = 2; node_map[3] = 1; 
        side_map[0] = 3; side_map[1] = 0; side_map[2] = 1; side_map[3] = 2; 
        cell_map[0] = 0; 
	break;

    default:
        fail("check_renumber") << "Invalid mesh type encountered." 
				    << endl;
	all_passed = false;
	return all_passed;
    }

    // Check node_map
    bool got_node_map = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
    {
        if (node_map[i] != mesh.get_nodes_map(i))
	    got_node_map = false;
    }
    if (!got_node_map)
    {
        fail(" Renumber node_map ") << "Renumber node_map not obtained." 
				    << endl;
 	all_passed = false;
    }
    // Check side_map
    bool got_side_map = true;
    for (int i = 0; i < mesh.get_dims_nsides(); i++)
    {
        if (side_map[i] != mesh.get_sides_map(i))
	    got_side_map = false;
    }
    if (!got_side_map)
    {
        fail(" Renumber side_map ") << "Renumber side_map not obtained." 
				    << endl;
 	all_passed = false;
    }
    // Check cell_map
    bool got_cell_map = true;
    for (int i = 0; i < mesh.get_dims_ncells(); i++)
    {
        if (cell_map[i] != mesh.get_cells_map(i))
	    got_cell_map = false;
    }
    if (!got_cell_map)
    {
        fail(" Renumber cell_map ") << "Renumber cell_map not obtained." 
				    << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" Renumber Accessors " ) << "Got all Renumber accessors." 
					 << endl;
    else
	fail(" Renumber Accessors ") << 
             "Errors in some Renumber accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_nodes(const rtt_meshReaders::RTT_Format & mesh, 
				 const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the nodes accessor functions for this mesh.
    bool all_passed = true;
    vector<vector<double> > coords(mesh.get_dims_nnodes(), 
				   vector<double>(mesh.get_dims_ndim(),0.0));
    vector<int> parents(mesh.get_dims_nnodes());
    vector<vector<int> > flags(mesh.get_dims_nnodes());

    switch (meshtype)
    {
    case DEFINED:
    case SORTED:
        // set node coords per the input deck.
	coords[2][1] = 2.0;
	// set node parents per the input deck.
	for (int i = 0; i < mesh.get_dims_nnodes(); i++)
	    parents[i] = i;
	// load the node flags.
	flags[0].push_back(11); flags[0].push_back(1); flags[0].push_back(101);
	flags[2].push_back(21); flags[2].push_back(4); flags[2].push_back(101);
	switch (meshtype)
	{
	case DEFINED:
            coords[1][2] = 3.0;
            coords[3][0] = 1.0;
            flags[1].push_back(21); flags[1].push_back(1); 
	        flags[1].push_back(101);
	    flags[3].push_back(6); flags[3].push_back(4); 
	        flags[3].push_back(22);
	    break;

	case SORTED: // Nodes 2 and 4 swapped.
            coords[1][0] = 1.0;
            coords[3][2] = 3.0;
            flags[1].push_back(6); flags[1].push_back(4); 
                flags[1].push_back(22);
            flags[3].push_back(21); flags[3].push_back(1); 
                flags[3].push_back(101);
            flags[1].push_back(6); flags[1].push_back(4); 
                flags[1].push_back(22);
            break;
	}
	break;

    default:
        fail("check_nodes") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check all of the node coords.
    if (coords != mesh.get_nodes_coords())
    {
        fail(" Nodes Coordinates ") << "Nodes coordinates not obtained." 
				    << endl;
	all_passed = false;
    }
    // Check all of the coordinate directions for a single node.
    bool got_node_coords = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
	if (coords[i] != mesh.get_nodes_coords(i))
	    got_node_coords = false;
    if (!got_node_coords)
    {
        fail(" Nodes Coordinates ") << "Node coordinates not obtained." 
				   << endl;
	all_passed = false;
    }
    // Check a single coordinate direction for a single node.
    bool got_node_coord = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
    {
        for (int d = 0; d < mesh.get_dims_ndim(); d++)
	    if (coords[i][d] != mesh.get_nodes_coords(i,d))
	        got_node_coord = false;
    }
    if (!got_node_coord)
    {
        fail(" Nodes Coordinate ") << "Node coordinate not obtained." << endl;
	all_passed = false;
    }
    // Retrieve node numbers from coordinates.
    bool got_nodes_node = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
	if (i != mesh.get_nodes_node(coords[i]))
	    got_nodes_node = false;
    if (!got_nodes_node)
    {
        fail(" Nodes node ") << "Nodes node not obtained." << endl;
	all_passed = false;
    }
    // Check the node parents.
    bool got_nodes_parents = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
	if (parents[i] != mesh.get_nodes_parents(i))
	    got_nodes_parents = false;
    if (!got_nodes_parents)
    {
        fail(" Nodes parents ") << "Nodes parents not obtained." << endl;
	all_passed = false;
    }
    // Check the node flags.
    bool got_nodes_flags = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
    {
        for (int f = 0; f < mesh.get_dims_nnode_flag_types(); f++)
	    if (flags[i][f] != mesh.get_nodes_flags(i,f))
	        got_nodes_flags = false;
    }
    if (!got_nodes_flags)
    {
        fail(" Nodes Flags ") << "Nodes flags not obtained." << endl;
	all_passed = false;
    }

    if (all_passed)
        pass(" Nodes Accessors " ) << "Got all Nodes accessors." << endl;
    else
	fail(" Nodes Accessors ") << "Errors in some Nodes accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_sides(const rtt_meshReaders::RTT_Format & mesh, 
				 const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the sides accessor functions for this mesh.
    bool all_passed = true;
    vector<int> sideType;
    vector<vector<int> > nodes;
    vector<vector<int> > flags;

    switch (meshtype)
    {
    case DEFINED:
	// set side sideType per the input deck.
	sideType.resize(mesh.get_dims_nsides(), 2);
	// load the side nodes.
	nodes.resize(mesh.get_dims_nsides());
	nodes[0].push_back(1); nodes[0].push_back(2); nodes[0].push_back(3); 
	nodes[1].push_back(0); nodes[1].push_back(3); nodes[1].push_back(2); 
	nodes[2].push_back(0); nodes[2].push_back(1); nodes[2].push_back(3); 
	nodes[3].push_back(0); nodes[3].push_back(2); nodes[3].push_back(1); 
	// load the side flags.
	flags.resize(mesh.get_dims_nsides());
	flags[0].push_back(2);
	flags[1].push_back(1);
	flags[2].push_back(1);
	flags[3].push_back(1);
	break;

    case SORTED: // (1->4, 2->1, 3->2, 4->3 is side translation from original)
	// set sorted side sideType per the node renumbering.
        // (1->4, 2->1, 3->2, 4->3 is side translation from original).
	sideType.resize(mesh.get_dims_nsides(), 2);
	// load the side nodess.
	nodes.resize(mesh.get_dims_nsides());
	nodes[0].push_back(0); nodes[0].push_back(1); nodes[0].push_back(2); 
	nodes[1].push_back(0); nodes[1].push_back(3); nodes[1].push_back(1); 
	nodes[2].push_back(0); nodes[2].push_back(2); nodes[2].push_back(3); 
	nodes[3].push_back(3); nodes[3].push_back(2); nodes[3].push_back(1); 
	// load the side flags.
	flags.resize(mesh.get_dims_nsides());
	flags[0].push_back(1);
	flags[1].push_back(1);
	flags[2].push_back(1);
	flags[3].push_back(2);
	break;

   default:
        fail("check_sides") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check the side sideType.
    bool got_sides_type = true;
    for (int i = 0; i < mesh.get_dims_nsides(); i++)
	if (sideType[i] != mesh.get_sides_type(i))
	    got_sides_type = false;
    if (!got_sides_type)
    {
        fail(" Side type ") << "Side type not obtained." << endl;
	all_passed = false;
    }
    // Check all of the side nodes.
    if (nodes != mesh.get_sides_nodes())
    {
        fail(" Sides nodes ") << "Sides nodes not obtained." << endl;
	all_passed = false;
    }
    // Check all of the nodes for a single side.
    bool got_side_nodes = true;
    for (int i = 0; i < mesh.get_dims_nsides(); i++)
	if (nodes[i] != mesh.get_sides_nodes(i))
	    got_side_nodes = false;
    if (!got_side_nodes)
    {
        fail(" Sides nodes ") << "Side nodes not obtained." << endl;
	all_passed = false;
    }
    // Check a single node for a single side.
    bool got_side_node = true;
    for (int i = 0; i < mesh.get_dims_nsides(); i++)
    {
        for (int n = 0; n < mesh.get_cell_defs_nnodes(mesh.get_sides_type(i));
	     n++)
	    if (nodes[i][n] != mesh.get_sides_nodes(i,n))
	        got_side_node = false;
    }
    if (!got_side_node)
    {
        fail(" Sides node ") << "Side node not obtained." << endl;
	all_passed = false;
    }
    // Check the side flags.
    bool got_sides_flags = true;
    for (int i = 0; i < mesh.get_dims_nsides(); i++)
    {
        for (int f = 0; f < mesh.get_dims_nside_flag_types(); f++)
	    if (flags[i][f] != mesh.get_sides_flags(i,f))
	        got_sides_flags = false;
    }
    if (!got_sides_flags)
    {
        fail(" Side Flags ") << "Side flags not obtained." << endl;
	all_passed = false;
    }

    if (all_passed)
        pass(" Sides Accessors " ) << "Got all Sides accessors." << endl;
    else
	fail(" Sides Accessors ") << "Errors in some Sides accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_cells(const rtt_meshReaders::RTT_Format & mesh, 
				 const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the cells functions for this mesh.
    bool all_passed = true;
    vector<int> cellType;
    vector<vector<int> > nodes;
    vector<vector<int> > flags;

    switch (meshtype)
    {
    case DEFINED:
    case SORTED: // Nothing changes for this case.
	// set cell cellType per the input deck.
	cellType.resize(mesh.get_dims_ncells(), 5);
	// load the cell flags.
	flags.resize(mesh.get_dims_ncells());
	flags[0].push_back(1); flags[0].push_back(2);
	// Load the cell nodes.
	nodes.resize(mesh.get_dims_ncells());
	nodes[0].push_back(0); nodes[0].push_back(1); nodes[0].push_back(2); 
	nodes[0].push_back(3);
	break;

    default:
        fail("check_cells") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check the cell cellType.
    bool got_cells_type = true;
    for (int i = 0; i < mesh.get_dims_ncells(); i++)
	if (cellType[i] != mesh.get_cells_type(i))
	    got_cells_type = false;
    if (!got_cells_type)
    {
        fail(" Cell type ") << "Cell type not obtained." << endl;
	all_passed = false;
    }
    // Check all of the cell nodes.
    if (nodes != mesh.get_cells_nodes())
    {
        fail(" Cells nodes ") << "Cells nodes not obtained." << endl;
	all_passed = false;
    }
    // Check all of the nodes for a single cell.
    bool got_cell_nodes = true;
    for (int i = 0; i < mesh.get_dims_ncells(); i++)
	if (nodes[i] != mesh.get_cells_nodes(i))
	    got_cell_nodes = false;
    if (!got_cell_nodes)
    {
        fail(" Cells nodes ") << "Cell nodes not obtained." << endl;
	all_passed = false;
    }
    // Check a single node for a single cell.
    bool got_cell_node = true;
    for (int i = 0; i < mesh.get_dims_ncells(); i++)
    {
        for (int n = 0; n < mesh.get_cell_defs_nnodes(mesh.get_cells_type(i));
	     n++)
	    if (nodes[i][n] != mesh.get_cells_nodes(i,n))
	        got_cell_node = false;
    }
    if (!got_cell_node)
    {
        fail(" Cells node ") << "Cell node not obtained." << endl;
	all_passed = false;
    }
    // Check the cell flags.
    bool got_cells_flags = true;
    for (int i = 0; i < mesh.get_dims_ncells(); i++)
    {
        for (int f = 0; f < mesh.get_dims_ncell_flag_types(); f++)
	    if (flags[i][f] != mesh.get_cells_flags(i,f))
	        got_cells_flags = false;
    }
    if (!got_cells_flags)
    {
        fail(" Cell Flags ") << "Cell flags not obtained." << endl;
	all_passed = false;
    }

    if (all_passed)
        pass(" Cells Accessors " ) << "Got all Cells accessors." << endl;
    else
	fail(" Cells Accessors ") << "Errors in some Cells accessors." << endl;

    return all_passed;
}
bool TestRTT_Format::check_node_data(const rtt_meshReaders::RTT_Format & mesh, 
				     const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the node_data functions for this mesh.
    bool all_passed = true;
    vector<vector<double> > data(mesh.get_dims_nnodes(), 
			    vector<double>(mesh.get_dims_nnode_data(),0.0));

    switch (meshtype)
    {
    case DEFINED:
        // set node data per the input deck.
	data[1][2] = 3.0;
	data[2][1] = 2.0;
	data[3][0] = 1.0;
	break;

    case SORTED:
        // set node data per the input deck.
	data[1][0] = 1.0;
	data[2][1] = 2.0;
	data[3][2] = 3.0;
	break;

    default:
        fail("check_node_data") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check all of the node data.
    if (data != mesh.get_node_data())
    {
        fail(" NodeData ") << "NodeData not obtained for all nodes/fields." 
			   << endl;
	all_passed = false;
    }
    // Check all of the data fields for a single node.
    bool got_node_data_fields = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
	if (data[i] != mesh.get_node_data(i))
	    got_node_data_fields = false;
    if (!got_node_data_fields)
    {
        fail(" NodeData ") << "NodeData fields not obtained for a node." 
			   << endl;
	all_passed = false;
    }
    // Check a single data field for a single node.
    bool got_node_data = true;
    for (int i = 0; i < mesh.get_dims_nnodes(); i++)
    {
        for (int d = 0; d < mesh.get_dims_nnode_data(); d++)
	    if (data[i][d] != mesh.get_node_data(i,d))
	        got_node_data = false;
    }
    if (!got_node_data)
    {
        fail(" NodeData ") << "NodeData value not obtained." << endl;
	all_passed = false;
    }

    if (all_passed)
        pass(" NodeData Accessors " ) << "Got all NodeData accessors." << endl;
    else
	fail(" NodeData Accessors ") << "Errors in some NodeData accessors."
				     << endl;

    return all_passed;
}
bool TestRTT_Format::check_side_data(const rtt_meshReaders::RTT_Format & mesh, 
				     const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the side_data functions for this mesh.
    bool all_passed = true;
    vector<vector<double> > data(mesh.get_dims_nsides(), 
			    vector<double>(mesh.get_dims_nside_data(),0.0));

    switch (meshtype)
    {
    case DEFINED:
        // set side data per the input deck.
	data[2][1] = 2.0;
	data[3][0] = 1.0;
	break;

    case SORTED:
        // set side data per the input deck.
	data[1][1] = 2.0;
	data[2][0] = 1.0;
	break;

    default:
        fail("check_side_data") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check all of the side data.
    if (data != mesh.get_side_data())
    {
        fail(" SideData ") << "SideData not obtained for all sides/fields." 
			   << endl;
	all_passed = false;
    }
    // Check all of the data fields for a single side.
    bool got_side_data_fields = true;
    for (int i = 0; i < mesh.get_dims_nsides(); i++)
	if (data[i] != mesh.get_side_data(i))
	    got_side_data_fields = false;
    if (!got_side_data_fields)
    {
        fail(" SideData ") << "SideData fields not obtained for a side." 
			   << endl;
	all_passed = false;
    }
    // Check a single data field for a single side.
    bool got_side_data = true;
    for (int i = 0; i < mesh.get_dims_nsides(); i++)
    {
        for (int d = 0; d < mesh.get_dims_nside_data(); d++)
	    if (data[i][d] != mesh.get_side_data(i,d))
	        got_side_data = false;
    }
    if (!got_side_data)
    {
        fail(" SideData ") << "SideData value not obtained." << endl;
	all_passed = false;
    }

    if (all_passed)
        pass(" SideData Accessors " ) << "Got all SideData accessors." << endl;
    else
	fail(" SideData Accessors ") << "Errors in some SideData accessors."
				     << endl;

    return all_passed;
}
bool TestRTT_Format::check_cell_data(const rtt_meshReaders::RTT_Format & mesh, 
				     const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the cell_data functions for this mesh.
    bool all_passed = true;
    vector<vector<double> > data;
    vector<double> cell_data(mesh.get_dims_ncell_data(),0.0);

    switch (meshtype)
    {
    case DEFINED:
    case SORTED: // Nothing changes for this case.
        // set cell data per the input deck.
	data.push_back(cell_data);
	break;

    default:
        fail("check_cell_data") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check all of the cell data.
    if (data != mesh.get_cell_data())
    {
        fail(" CellData ") << "CellData not obtained for all cells/fields." 
			   << endl;
	all_passed = false;
    }
    // Check all of the data fields for a single cell.
    bool got_cell_data_fields = true;
    for (int i = 0; i < mesh.get_dims_ncells(); i++)
	if (data[i] != mesh.get_cell_data(i))
	    got_cell_data_fields = false;
    if (!got_cell_data_fields)
    {
        fail(" CellData ") << "CellData fields not obtained for a cell." 
			   << endl;
	all_passed = false;
    }
    // Check a single data field for a single cell.
    bool got_cell_data = true;
    for (int i = 0; i < mesh.get_dims_ncells(); i++)
    {
        for (int d = 0; d < mesh.get_dims_ncell_data(); d++)
	    if (data[i][d] != mesh.get_cell_data(i,d))
	        got_cell_data = false;
    }
    if (!got_cell_data)
    {
        fail(" CellData ") << "CellData value not obtained." << endl;
	all_passed = false;
    }

    if (all_passed)
        pass(" CellData Accessors " ) << "Got all CellData accessors." << endl;
    else
	fail(" CellData Accessors ") << "Errors in some CellData accessors."
				     << endl;

    return all_passed;
}
bool TestRTT_Format::check_virtual(const rtt_meshReaders::RTT_Format & mesh, 
				   const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the virtual accessor functions for this mesh.
    bool all_passed = true;
    vector<vector<double> > node_coords;
    string node_coord_units;
    vector<vector<int> > element_nodes;
    vector<rtt_meshReaders::Element_Definition::Element_Type> element_types;
    map<string, set<int> > node_sets;
    map<string, set<int> > element_sets;
    string title;
    vector<double> coords(mesh.get_dims_ndim(),0.0);
    vector<int> side_nodes;
    set<int> flag_nodes;
    set<int> flag_elements;

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
	element_types.push_back(rtt_meshReaders::Element_Definition::NODE);
	element_types.push_back(rtt_meshReaders::Element_Definition::BAR_2);
	element_types.push_back(rtt_meshReaders::Element_Definition::TRI_3); 
	element_types.push_back(rtt_meshReaders::Element_Definition::QUAD_4);
	element_types.push_back(rtt_meshReaders::Element_Definition::PYRA_5);
	element_types.push_back(rtt_meshReaders::Element_Definition::TETRA_4);
	element_types.push_back(rtt_meshReaders::Element_Definition::PENTA_6);
	element_types.push_back(rtt_meshReaders::Element_Definition::HEXA_8); 
	// load the node sets
	flag_nodes.insert(0);
	node_sets.insert(make_pair("node_type/interior", flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(1); flag_nodes.insert(2);
	node_sets.insert(make_pair("node_type/dudded", flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(3);
	node_sets.insert(make_pair("node_type/parent", flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(0); flag_nodes.insert(1);
	node_sets.insert(make_pair("boundary/reflective", flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(2); flag_nodes.insert(3);
	node_sets.insert(make_pair("boundary/vacuum", flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(0); flag_nodes.insert(1); flag_nodes.insert(2);
	node_sets.insert(make_pair("source/no_source", flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	flag_nodes.insert(3);
	node_sets.insert(make_pair("source/rad_source", flag_nodes));
	flag_nodes.erase(flag_nodes.begin(),flag_nodes.end());
	// load the element (i.e., sides + cell) sets
	flag_elements.insert(1); flag_elements.insert(2); 
	    flag_elements.insert(3);
	element_sets.insert(make_pair("boundary/reflective", flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	flag_elements.insert(0);
	element_sets.insert(make_pair("boundary/vacuum", flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	flag_elements.insert(4);
	element_sets.insert(make_pair("material/control_rod", flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	element_sets.insert(make_pair("material/shield", flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	element_sets.insert(make_pair("rad_source/src_name1", flag_elements));
	flag_elements.insert(4);
	element_sets.insert(make_pair("rad_source/src_name2", flag_elements));
	flag_elements.erase(flag_elements.begin(),flag_elements.end());
	// set the mesh title
	title = "RTT_format mesh file definition, version 7.";
	break;

    default:
        fail("check_virtual") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }
    // Check node coords
    if (node_coords != mesh.get_node_coords())
    {
        fail(" Node Coordinates ") << "Node coordinates not obtained." << endl;
	all_passed = false;
    }
    // Check coordinate units.
    if (mesh.get_dims_coor_units() != mesh.get_node_coord_units())
    {
        fail(" Coordinates Units ") << "Coordinate units not obtained." 
				    << endl;
 	all_passed = false;
    }
    if (element_nodes != mesh.get_element_nodes())
    {
        fail(" Element Nodes ") << "Element nodes not obtained." << endl;
 	all_passed = false;
    }
    // Check Element Types.
    if (element_types != mesh.get_element_types())
    {
	fail(" Element Types ") << "Element Types not obtained." << endl;
 	all_passed = false;
    }
    // Check node sets.
    if (node_sets != mesh.get_node_sets())
    {
        fail(" Node Sets ") << "Node sets not obtained." << endl;
 	all_passed = false;
    }
    // Check Element sets.
    if (element_sets != mesh.get_element_sets())
    {
        fail(" Element Sets ") << "Element sets not obtained." << endl;
 	all_passed = false;
    }
    // Check title.
    if (title != mesh.get_title())
    {
        fail(" Title ") << "Title not obtained." << endl;
 	all_passed = false;
    }
    // Check invariant.
    if (!mesh.invariant())
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
bool TestRTT_Format::check_connectivity(const rtt_meshReaders::RTT_Format & 
					mesh, const Meshes & meshtype)
{
    // Return if the Dims data is corrupt.        
    if (!verify_Dims(mesh, meshtype))
        return false;

    // Exercise the connectivity accessor functions for this mesh.
    bool all_passed = true;
    vector<vector<vector<int> > > adjacent_cell;
    multimap<int, int> bndryFaces;
    multimap<int, vector<int> > Side_to_Cell_Face;
    // Rest of these variables are for the user to specify the tests to be
    // performed (i.e., the user has to initialize them in a case).
    vector<int> cells_to_check; // Probably don't want to check every cell.
    vector<int> sides_to_check; // Probably don't want to check every side.
    vector<int> bndryFacesCount;
    vector<int> cell_face(2,0);
    
    switch (meshtype)
    {
    case DEFINED:
        cells_to_check.push_back(0); // there is only one cell in the mesh.
	adjacent_cell.resize(cells_to_check.size());
	// Resize the adjacent cell vector for the number of sides per cell.
	for (int c = 0; c < cells_to_check.size(); c++)
	{
	    int cell = cells_to_check[c];
	    int cell_type = mesh.get_cells_type(cell);
	    int cell_sides = mesh.get_cell_defs_nsides(cell_type);
	    for (int s = 0; s < cell_sides; s++)
	        sides_to_check.push_back(s); // only four sides
	    adjacent_cell[c].resize(cell_sides);
	}
	// All of the "adjacent cells" are boundaries for this case, and the
	// adjacent cell is returned as the negative of the boundary flag #. 
	adjacent_cell[0][0].push_back(-2);
	adjacent_cell[0][1].push_back(-1); adjacent_cell[0][2].push_back(-1);
	adjacent_cell[0][3].push_back(-1);
	// bndryFaces pairs are keyed on face number (i.e., 0 - 5 for a hex) to
	// the cell number as a value and contain all the faces that are either
	// on the outer geometric boundary of the problem or a junction between
	// cells with differing refinement levels in an AMR mesh. Only the
	// subset of the multimap that is entered here will be checked. 
	bndryFaces.insert(make_pair(0,0)); bndryFaces.insert(make_pair(1,0));
	bndryFaces.insert(make_pair(2,0)); bndryFaces.insert(make_pair(3,0));
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
	Side_to_Cell_Face.insert(make_pair(0,cell_face));
	cell_face[1] = 1;
	Side_to_Cell_Face.insert(make_pair(1,cell_face));
	cell_face[1] = 2;
	Side_to_Cell_Face.insert(make_pair(2,cell_face));
	cell_face[1] = 3;
	Side_to_Cell_Face.insert(make_pair(3,cell_face));
	break;

    default:
        fail("check_connectivity") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }

    // Check the number of cells that are adjacent to the cells that are
    // specified in cells_to_check.
    bool got_adj_cell_size = true;
    for (int c = 0; c < cells_to_check.size(); c++)
    {
        int cell = cells_to_check[c];
	for (int f = 0; f < adjacent_cell[cell].size(); f++)
	    if (adjacent_cell[cell][f].size() != mesh.get_adjCell_size(cell,f))
	        got_adj_cell_size = false;
    }
    if (!got_adj_cell_size)
    {
        fail(" Connectivity ") << 
	     "Connectivity adjacent cell size not obtained." << endl;
 	all_passed = false;
    }

    // Check the cell faces that are adjacent to the cells that are
    //  specified in cells_to_check.
    bool got_adj_cell = true;
    for (int c = 0; c < cells_to_check.size(); c++)
    {
        int cell = cells_to_check[c];
	for (int f = 0; f < adjacent_cell[cell].size(); f++)
	    for (int a = 0; a < adjacent_cell[cell][f].size(); a++) 
	        if (adjacent_cell[cell][f][a] != mesh.get_adjCell(cell,f,a))
	        got_adj_cell = false;
    }
    if (!got_adj_cell)
    {
        fail(" Connectivity ") << "Connectivity adjacent cell not obtained." 
			       << endl;
 	all_passed = false;
    }

    // Check the count of boundary cell faces.
    bool got_bndry_face_count = true;
    for (int f = 0; f < bndryFacesCount.size(); f++)
        if (bndryFacesCount[f] != mesh.get_bndryFaces_count(f))
	    got_bndry_face_count = false;
    if (!got_bndry_face_count)
    {
        fail(" Connectivity ") << 
	     "Connectivity boundary face count not obtained." << endl;
 	all_passed = false;
    }

    // Check the get_bndryCells accessor function for those face/cell pairs
    // that were specified for this case (in bndryFaces for the switch).
    bool got_bndry_cells = true;
    for (multimap<int, int>::iterator f = bndryFaces.begin();
	 f != bndryFaces.end(); f++)
	    if (mesh.get_bndryCells(f->first).count(f->second) == 0)
	        got_bndry_cells = false;
    if (!got_bndry_cells)
    {
        fail(" Connectivity ") << "Connectivity boundary cells not obtained." 
			       << endl;
 	all_passed = false;
    }

    // Check the get_bndryFace accessor function for those face/cell pairs
    // that were specified for this case (in bndryFaces for the switch).
    bool got_bndry_face = true;
    for (multimap<int, int>::iterator f = bndryFaces.begin();
	 f != bndryFaces.end(); f++)
	    if (!mesh.check_bndryFace(f->second, f->first))
	        got_bndry_face = false;
    if (!got_bndry_face)
    {
        fail(" Connectivity ") << "Connectivity boundary face not obtained." 
			       << endl;
 	all_passed = false;
    }

    // Check the get_Cell_from_Side accessor function for those sides
    // that were specified for this case (Side_to_Cell_Face in the switch).
    bool got_cell_from_side = true;
    for (multimap<int, vector<int> >::iterator s = Side_to_Cell_Face.begin();
	 s != Side_to_Cell_Face.end(); s++)
	    if (mesh.get_Cell_from_Side(s->first) != s->second[0])
	        got_cell_from_side = false;
    if (!got_cell_from_side)
    {
        fail(" Connectivity ") << "Connectivity cell not obtained from side." 
			       << endl;
 	all_passed = false;
    }

    // Check the get_Cell_Face_from_Side accessor function for those sides
    // that were specified for this case (Side_to_Cell_Face in the switch).
    bool got_cell_face_from_side = true;
    for (multimap<int, vector<int> >::iterator s = Side_to_Cell_Face.begin();
	 s != Side_to_Cell_Face.end(); s++)
	    if (mesh.get_Cell_Face_from_Side(s->first) != s->second[1])
	        got_cell_face_from_side = false;
    if (!got_cell_face_from_side)
    {
        fail(" Connectivity ") << 
	     "Connectivity cell face not obtained from side." << endl;
 	all_passed = false;
    }

    if (all_passed)
        pass(" Connectivity Accessors " ) << 
	     "Got all connectivity accessors." << endl;
    else
	fail(" Connectivity Accessors ") << 
	     "Errors in some connectivity accessors." << endl;

    return all_passed;
}

} // end namespace rtt_meshReaders_test


//---------------------------------------------------------------------------//
//                              end of TestRTT_Format.cc
//---------------------------------------------------------------------------//
