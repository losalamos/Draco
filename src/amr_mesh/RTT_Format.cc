//----------------------------------*-C++-*--------------------------------//
// RTT_Format.cc
// Shawn Pautz (TychoMesh.cc original) / B.T. Adams (Extended to RTT_Format.cc)
//
// The node, side, and cell numbering that is implemented in ICEM appears to be
// completely arbitrary. A more structured arrangement is preferrable for some
// applications (notably mine), so an option was added to renumber the mesh so
// that the first node, side, and cell numbers correspond to the lowest x,y,and
// z coordinates and the last correspond to the highest x,y,z coordinates (with
// x varying fastest, followed by y, and finally z). The cell definitions are
// also modified so that, in a hexahedron, the lowest nodes correspond to side
// "-z", followed by "-y", "-x", "x", "y", and "z". This cell definition has
// the distinct advantage that, in Cartesian coordinates with hexahedrals, the
// ICEM cell definitions can be translated to this coordinate system simply by
// sorting the cell's nodes into ascending order, regardless of the particular 
// cell orientations that were used in ICEM. Sorting the sides based upon node
// numbers also yields these definitions. The constructor for the RTT_Format 
// class defaults to no renumbering, but this option was added to this class 
// because calculating the connectivity is CPU intensive and this addition 
// prevents having to repeat the work. Runtime testing indicates that there is
// a significant DECREASE in the time required to read and connect the mesh if
// a continuous adaptive refinement mesh is used with a few thousand cells or 
// more and sorting is implemented. This behaviour results from the fact that
// a bilinear search routine can be used connect the refined cells adjacent to
// an area of lesser refinement if sorting has been performed, while a linear
// search of all the nodes must be performed otherwise.
// 7 June 99
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "RTT_Format.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <utility>
#include <algorithm>
#include <iterator>

using std::string;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::getline;
using std::multimap;
using std::pair;
using std::map;
using std::set;
using std::make_pair;
using std::vector;
using std::find;
using std::sort;
using std::min;
using std::max;
using std::greater;

namespace rtt_format
{

RTT_Format::RTT_Format(const string & RTT_File, const bool & renumber)
{
    if (C4::node() == 0)
    {
	readMesh(RTT_File,renumber);
	calculateConnectivity();
    }
}

void RTT_Format::calculateConnectivity()
{
    spConnectivity = new Connectivity(dims, * spCellDefs, * spCells, 
				      * spSides, * spNodes);
}

void RTT_Format::readMesh(const string & RTT_File, const bool & renumber)
{
    if (C4::node() == 0)
    {
        const char * file = RTT_File.c_str();
	ifstream meshfile(file, ios::in);
	if (!meshfile)
	{
	    cout << "File could not be opened\n";
	}

	try
	{
	    readKeyword(meshfile);
	    header.readHeader(meshfile);
	    dims.readDims(meshfile,renumber);
	    createMembers();
	    readFlagBlocks(meshfile);
	    readDataIDs(meshfile);
	    spCellDefs->readCellDefs(meshfile);
	    spNodes->readNodes(meshfile);
	    spSides->readSides(meshfile);
	    spCells->readCells(meshfile);
	    spNodeData->readNodeData(meshfile);
	    spSideData->readSideData(meshfile);
	    spCellData->readCellData(meshfile);
	    readEndKeyword(meshfile);
	}
	catch (dsxx::assertion as)
	{
	    cout << "Assertion thrown: " << as.what() << endl;
	    Insist(false, as.what());
	}
    }
}

void RTT_Format::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "rtt_ascii", "Invalid mesh file: Not an RTT file");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::createMembers()
{
    spNodeFlags = new NodeFlags(dims);
    spSideFlags = new SideFlags(dims);
    spCellFlags = new CellFlags(dims);
    spNodeDataIds = new NodeDataIDs(dims);
    spSideDataIds = new SideDataIDs(dims);
    spCellDataIds = new CellDataIDs(dims);
    spCellDefs = new CellDefs(dims);
    spNodes = new Nodes(* spNodeFlags, dims);
    spSides = new Sides(* spSideFlags, dims, * spCellDefs, * spNodes);
    spCells = new Cells(* spCellFlags, dims, * spCellDefs, * spNodes);
    spNodeData = new NodeData(dims, * spNodes);
    spSideData = new SideData(dims, * spNodes);
    spCellData = new CellData(dims, * spNodes);
}

void RTT_Format::readFlagBlocks(ifstream & meshfile)
{
    spNodeFlags->readNodeFlags(meshfile);
    spSideFlags->readSideFlags(meshfile);
    spCellFlags->readCellFlags(meshfile);
}

void RTT_Format::readDataIDs(ifstream & meshfile)
{
    spNodeDataIds->readDataIDs(meshfile);
    spSideDataIds->readDataIDs(meshfile);
    spCellDataIds->readDataIDs(meshfile);
}

void RTT_Format::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_rtt_mesh",
	   "Invalid mesh file: RTT file missing end");
    getline(meshfile, dummyString);
}

void RTT_Format::Header::readHeader(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::Header::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "header", "Invalid mesh file: Header block missing");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::Header::readData(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString >> version;
    Insist(dummyString == "version",
	   "Invalid mesh file: Header block missing version");
    Insist(version == "v1.0.0", "Invalid mesh file: Wrong version");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> title;
    Insist(dummyString == "title",
	   "Invalid mesh file: Header block missing title");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> date;
    Insist(dummyString == "date",
	   "Invalid mesh file: Header block missing date");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> cycle;
    Insist(dummyString == "cycle",
	   "Invalid mesh file: Header block missing cycle");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> time;
    Insist(dummyString == "time",
	   "Invalid mesh file: Header block missing time");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> ncomments;
    Insist(dummyString == "ncomments",
	   "Invalid mesh file: Header block missing ncomments");
    comments.redim(ncomments);
    for (int i = 0; i < ncomments; ++i)
	getline(meshfile, comments(i));
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::Header::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_header",
	   "Invalid mesh file: Header block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::Dims::readDims(ifstream & meshfile, const bool & renumber_)
{
    renumber = renumber_;
    readKeyword(meshfile);
    readUnits(meshfile);
    readCellDefs(meshfile);
    readDimensions(meshfile);
    readNodes(meshfile);
    readSides(meshfile);
    readCells(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::Dims::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "dims", "Invalid mesh file: Dimension block missing");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::Dims::readUnits(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString >> coor_units;
    Insist(dummyString == "coor_units",
	   "Invalid mesh file: Dimension block missing coor_units");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> prob_time_units;
    Insist(dummyString == "prob_time_units",
	   "Invalid mesh file: Dimension block missing prob_time_units");
    getline(meshfile, dummyString);
}

void RTT_Format::Dims::readCellDefs(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString >> ncell_defs;
    Insist(dummyString == "ncell_defs",
	   "Invalid mesh file: Dimension block missing ncell_defs");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nnodes_max;
    Insist(dummyString == "nnodes_max",
	   "Invalid mesh file: Dimension block missing nnodes_max");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nsides_max;
    Insist(dummyString == "nsides_max",
	   "Invalid mesh file: Dimension block missing nsides_max");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nnodes_side_max;
    Insist(dummyString == "nnodes_side_max",
	   "Invalid mesh file: Dimension block missing nnodes_side_max");
    getline(meshfile, dummyString);
}

void RTT_Format::Dims::readDimensions(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString >> ndim;
    Insist(dummyString == "ndim",
	   "Invalid mesh file: Dimension block missing ndim");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> ndim_topo;
    Insist(dummyString == "ndim_topo",
	   "Invalid mesh file: Dimension block missing ndim_topo");
    getline(meshfile, dummyString);
}

void RTT_Format::Dims::readNodes(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString >> nnodes;
    Insist(dummyString == "nnodes",
	   "Invalid mesh file: Dimension block missing nnodes");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nnode_flag_types;
    Insist(dummyString == "nnode_flag_types",
	   "Invalid mesh file: Dimension block missing nnode_flag_types");
    getline(meshfile, dummyString);

    nnode_flags.redim(nnode_flag_types);
    meshfile >> dummyString;
    Insist(dummyString == "nnode_flags",
	   "Invalid mesh file: Dimension block missing nnode_flags");
    for (dsxx::Mat1<int>::iterator iter = nnode_flags.begin();
	 iter < nnode_flags.end(); ++iter)
	meshfile >> *iter;
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nnode_data;
    Insist(dummyString == "nnode_data",
	   "Invalid mesh file: Dimension block missing nnode_data");
    getline(meshfile, dummyString);
}

void RTT_Format::Dims::readSides(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString >> nsides;
    Insist(dummyString == "nsides",
	   "Invalid mesh file: Dimension block missing nsides");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nside_types;
    Insist(dummyString == "nside_types",
	   "Invalid mesh file: Dimension block missing nside_types");
    getline(meshfile, dummyString);

    side_types.redim(nside_types);
    meshfile >> dummyString;
    Insist(dummyString == "side_types",
	   "Invalid mesh file: Dimension block missing side_types");
    for (dsxx::Mat1<int>::iterator iter = side_types.begin();
	 iter < side_types.end(); ++iter)
    {
	meshfile >> *iter;
	--(*iter);
    }
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nside_flag_types;
    Insist(dummyString == "nside_flag_types",
	   "Invalid mesh file: Dimension block missing nside_flag_types");
    getline(meshfile, dummyString);

    nside_flags.redim(nside_flag_types);
    meshfile >> dummyString;
    Insist(dummyString == "nside_flags",
	   "Invalid mesh file: Dimension block missing nside_flags");
    for (dsxx::Mat1<int>::iterator iter = nside_flags.begin();
	 iter < nside_flags.end(); ++iter)
	meshfile >> *iter;
    getline(meshfile, dummyString);

    meshfile >> dummyString >> nside_data;
    Insist(dummyString == "nside_data",
	   "Invalid mesh file: Dimension block missing nside_data");
    getline(meshfile, dummyString);
}

void RTT_Format::Dims::readCells(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString >> ncells;
    Insist(dummyString == "ncells",
	   "Invalid mesh file: Dimension block missing ncells");
    getline(meshfile, dummyString);

    meshfile >> dummyString >> ncell_types;
    Insist(dummyString == "ncell_types",
	   "Invalid mesh file: Dimension block missing ncell_types");
    getline(meshfile, dummyString);

    cell_types.redim(ncell_types);
    meshfile >> dummyString;
    Insist(dummyString == "cell_types",
	   "Invalid mesh file: Dimension block missing cell_types");
    for (dsxx::Mat1<int>::iterator iter = cell_types.begin();
	 iter < cell_types.end(); ++iter)
    {
	meshfile >> *iter;
	--(*iter);
    }
    getline(meshfile, dummyString);

    meshfile >> dummyString >> ncell_flag_types;
    Insist(dummyString == "ncell_flag_types",
	   "Invalid mesh file: Dimension block missing ncell_flag_types");
    getline(meshfile, dummyString);

    ncell_flags.redim(ncell_flag_types);
    meshfile >> dummyString;
    Insist(dummyString == "ncell_flags",
	   "Invalid mesh file: Dimension block missing ncell_flags");
    for (dsxx::Mat1<int>::iterator iter = ncell_flags.begin();
	 iter < ncell_flags.end(); ++iter)
	meshfile >> *iter;
    getline(meshfile, dummyString);

    meshfile >> dummyString >> ncell_data;
    Insist(dummyString == "ncell_data",
	   "Invalid mesh file: Dimension block missing ncell_data");
    getline(meshfile, dummyString);
}

void RTT_Format::Dims::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_dims",
	   "Invalid mesh file: Dimension block missing end_dims");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::NodeFlags::readNodeFlags(ifstream & meshfile)
{
    readKeyword(meshfile);
    readFlagTypes(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::NodeFlags::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "node_flags",
	   "Invalid mesh file: node_flags block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::NodeFlags::readFlagTypes(ifstream & meshfile)
{
    int flagTypeNum;
    string dummyString;

    for (int i = 0; i < dims.get_nnode_flag_types(); ++i)
    {
	meshfile >> flagTypeNum >> dummyString;
	// Ignore plurals in node flag definitions
	if (dummyString[dummyString.size()-1] == 's' || 
	    dummyString[dummyString.size()-1] == 'S')
	    dummyString.resize(dummyString.size()-1);
	Insist(flagTypeNum == i+1,
	       "Invalid mesh file: node flag type out of order");
	flagTypes(i) = new Flags(dims.get_nnode_flags(i), dummyString);
	getline(meshfile, dummyString);
	flagTypes(i)->readFlags(meshfile);
    }
}

void RTT_Format::NodeFlags::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_node_flags",
	   "Invalid mesh file: node_flags block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::SideFlags::readSideFlags(ifstream & meshfile)
{
    readKeyword(meshfile);
    readFlagTypes(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::SideFlags::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "side_flags",
	   "Invalid mesh file: side_flags block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::SideFlags::readFlagTypes(ifstream & meshfile)
{
    int flagTypeNum;
    string dummyString;

    for (int i = 0; i < dims.get_nside_flag_types(); ++i)
    {
	meshfile >> flagTypeNum >> dummyString;
	// Ignore plurals in side flag definitions
	if (dummyString[dummyString.size()-1] == 's' || 
	    dummyString[dummyString.size()-1] == 'S')
	    dummyString.resize(dummyString.size()-1);
	Insist(flagTypeNum == i+1,
	       "Invalid mesh file: side flag type out of order");
	flagTypes(i) = new Flags(dims.get_nside_flags(i), dummyString);
	getline(meshfile, dummyString);
	flagTypes(i)->readFlags(meshfile);
    }
}

void RTT_Format::SideFlags::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_side_flags",
	   "Invalid mesh file: side_flags block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::CellFlags::readCellFlags(ifstream & meshfile)
{
    readKeyword(meshfile);
    readFlagTypes(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::CellFlags::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "cell_flags",
	   "Invalid mesh file: cell_flags block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::CellFlags::readFlagTypes(ifstream & meshfile)
{
    int flagTypeNum;
    string dummyString;

    for (int i = 0; i < dims.get_ncell_flag_types(); ++i)
    {
	meshfile >> flagTypeNum >> dummyString;
	// Ignore plurals in cell flag definitions
	if (dummyString[dummyString.size()-1] == 's' || 
	    dummyString[dummyString.size()-1] == 'S')
	    dummyString.resize(dummyString.size()-1);
	Insist(flagTypeNum == i+1,
	       "Invalid mesh file: cell flag type out of order");
	flagTypes(i) = new Flags(dims.get_ncell_flags(i), dummyString);
	getline(meshfile, dummyString);
	flagTypes(i)->readFlags(meshfile);
    }
}

void RTT_Format::CellFlags::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_cell_flags",
	   "Invalid mesh file: cell_flags block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::Flags::readFlags(ifstream & meshfile)
{
    string dummyString;

    for (int i = 0; i < nflags; ++i)
    {
	meshfile >> flag_nums(i) >> flag_names(i);
	getline(meshfile, dummyString);
    }
}

void RTT_Format::NodeDataIDs::readDataIDs(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::NodeDataIDs::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "node_data_ids",
	   "Invalid mesh file: node_data_ids block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::NodeDataIDs::readData(ifstream & meshfile)
{
    int dataIDNum;
    string dummyString;

    for (int i = 0; i < dims.get_nnode_data(); ++i)
    {
	meshfile >> dataIDNum >> names(i) >> units(i);
	Insist(dataIDNum == i+1,
	       "Invalid mesh file: node data ID out of order");
	getline(meshfile, dummyString);
    }
}

void RTT_Format::NodeDataIDs::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_node_data_ids",
	   "Invalid mesh file: node_data_ids block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::SideDataIDs::readDataIDs(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::SideDataIDs::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "side_data_ids",
	   "Invalid mesh file: side_data_ids block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::SideDataIDs::readData(ifstream & meshfile)
{
    int dataIDNum;
    string dummyString;

    for (int i = 0; i < dims.get_nside_data(); ++i)
    {
	meshfile >> dataIDNum >> names(i) >> units(i);
	Insist(dataIDNum == i+1,
	       "Invalid mesh file: side data ID out of order");
	getline(meshfile, dummyString);
    }
}

void RTT_Format::SideDataIDs::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_side_data_ids",
	   "Invalid mesh file: side_data_ids block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::CellDataIDs::readDataIDs(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::CellDataIDs::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "cell_data_ids",
	   "Invalid mesh file: cell_data_ids block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::CellDataIDs::readData(ifstream & meshfile)
{
    int dataIDNum;
    string dummyString;

    for (int i = 0; i < dims.get_ncell_data(); ++i)
    {
	meshfile >> dataIDNum >> names(i) >> units(i);
	Insist(dataIDNum == i+1,
	       "Invalid mesh file: cell data ID out of order");
	getline(meshfile, dummyString);
    }
}

void RTT_Format::CellDataIDs::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_cell_data_ids",
	   "Invalid mesh file: cell_data_ids block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::CellDefs::readCellDefs(ifstream & meshfile)
{
    readKeyword(meshfile);
    readDefs(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::CellDefs::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "cell_defs",
	   "Invalid mesh file: cell_defs block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::CellDefs::readDefs(ifstream & meshfile)
{
    int cellDefNum;
    string dummyString;

    for (int i = 0; i < dims.get_ncell_defs(); ++i)
    {
	meshfile >> cellDefNum >> dummyString;
	Insist(cellDefNum == i+1,
	       "Invalid mesh file: cell def out of order");
	// Ignore plurals in cell definitions
	if (dummyString[dummyString.size()-1] == 's')
	    dummyString.resize(dummyString.size()-1);
	defs(i) = new CellDef(*this, dummyString);
	getline(meshfile, dummyString);
	defs(i)->readDef(meshfile);
    }
    if (dims.get_renumber())
        sortData();
}

void RTT_Format::CellDefs::sortData()
{
    for (int i = 0; i < dims.get_ncell_defs(); ++i)
    {
	defs(i)->sortData();
    }
}

void RTT_Format::CellDefs::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_cell_defs",
	   "Invalid mesh file: cell_defs block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::CellDef::readDef(ifstream & meshfile)
{
    string dummyString;

    meshfile >> nnodes >> nsides;
    side_types.redim(nsides);
    sides.redim(nsides);
    ordered_sides.resize(nsides);
    getline(meshfile, dummyString);

    for (int i = 0; i < nsides; ++i)
    {
	meshfile >> side_types(i);
	--side_types(i);
    }
    if (nsides > 0)
	getline(meshfile, dummyString);

    // note that this implementation does not preserve the "right hand rule"
    // of the cell definitions due to the use of a set container (which is 
    // sorted). It is slicker than snail snot when it comes time to implement 
    // the connectivity, however. The ordered_sides vector was added to allow
    // the original ordered data to be retained.
    int side;
    for (int i = 0; i < nsides; ++i)
    {
        int numb_nodes = cellDefs.get_cell_def(side_types(i)).get_nnodes();
        ordered_sides[i].resize(numb_nodes);
	for (int j = 0; j < numb_nodes; ++j)
	{
	    meshfile >> side;
	    --side;
	    sides(i).insert(side);
	    ordered_sides[i][j] = side;
	}
	if (sides(i).size() > 0)
	    getline(meshfile, dummyString);
    }
}

void RTT_Format::CellDef::sortData()
{
    // The cell definitions built into ICEM do not correspond to our node
    // numbering scheme. Define the alternative cell definitions herein, 
    // while retaining the right hand rule of cell definition. This section 
    // takes heavy advantage of the C++ rules of truncating integer division.
    // Note that our cell definition scheme has the distinct advantage that
    // the nodes are inherently in increasing integer order. This eliminates
    // any need to define a coordinate transform from the user-input cell
    // definitions.  

    ordered_sides.resize(nsides);
    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;

    for (int i = 0; i < nsides; ++i)
    {

        int numb_nodes = cellDefs.get_cell_def(side_types(i)).get_nnodes();
	sides(i).erase(sides(i).begin(),sides(i).end());
        ordered_sides[i].resize(numb_nodes);

	if (name == "line")
        {
	    ordered_sides[i][0] = i;
	    sides(i).insert(ordered_sides[i][0]);
	    if (debugging)
	    {
	        cout << name << endl;
		cout << i << " " << ordered_sides[i][0] << endl;
	    }
	}

	else if (name == "triangle")
	{
	    ordered_sides[i][0] = (i + 1)%2 + i/2;
	    ordered_sides[i][1] = (3 - i)%3;
	    sides(i).insert(ordered_sides[i][0]);
	    sides(i).insert(ordered_sides[i][1]);
	    if (debugging)
	    {
	        cout << name << endl;
		cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1] << endl;
	    }
	}

	else if (name == "quad")
	{
	    ordered_sides[i][0] = (i + 1)%2 + 2 * (i/2);
	    ordered_sides[i][1] = i%2 + (i + 1)/2;
	    sides(i).insert(ordered_sides[i][0]);
	    sides(i).insert(ordered_sides[i][1]);
	    if (debugging)
	    {
	        cout << name << endl;
		cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1] << endl;
	    }
	}

	else if (name == "tetrahedron")
	{
	    ordered_sides[i][0] = i/3;
	    ordered_sides[i][1] = 1 + (3 - i)%3 + 2 * (i/3);
	    ordered_sides[i][2] = (2 - i)%3 + 3 * (i/2);
	    sides(i).insert(ordered_sides[i][0]);
	    sides(i).insert(ordered_sides[i][1]);
	    sides(i).insert(ordered_sides[i][2]);
	    if (debugging)
	    {
	        cout << name << endl;
		cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2] << endl;
	    }
	}

	else if (name == "quad_pyr")
	{
	    if (i != 0)
	    {
	        // determine the index for the triangular side type
	        int j = 0;
		string side_name = cellDefs.get_cell_def(j).get_name();
		while (side_name != "triangle")
		{
		    j++;
		    Insist(j < cellDefs.get_ncell_defs(), 
			   "triangle side type not found for quad pyramid!")
		    side_name = 
		         cellDefs.get_cell_def(j).get_name();
		}
	        ordered_sides[i].resize(3);
		side_types(i) = j;
	    }
	    else
	    {
	        // determine the index for the quad side type
	        int j = 0;
		string side_name = cellDefs.get_cell_def(j).get_name();
		while (side_name != "quad")
		{
		    j++;
		    Insist(j < cellDefs.get_ncell_defs(), 
			   "quad side type not found for quad pyramid!")
		    side_name = 
		        cellDefs.get_cell_def(j).get_name();
		}
	        ordered_sides[i].resize(4);
		side_types(i) = j;
		ordered_sides[i][3] = 2;
	        sides(i).insert(ordered_sides[i][3]);
	    }
	    ordered_sides[i][0] = i/3 + i/4 ;
	    ordered_sides[i][1] = (i + 1)%2 * (1 + i/2) + 4 * (i%2);
	    ordered_sides[i][2] = (i + 1)%2 * (3 + i/2 - i/4) + i * (i%2);
	    sides(i).insert(ordered_sides[i][0]);
	    sides(i).insert(ordered_sides[i][1]);
	    sides(i).insert(ordered_sides[i][2]);
	    if (debugging)
	    {
	        cout << name << endl;
		cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2];
		if (i == 0)
		    cout << " " << ordered_sides[i][3];
		cout <<  endl;
	    }
	}

	else if (name == "tri_prism")
	{
	    if (i == 0 || i == nsides - 1)
	    {
	        // determine the index for the triangular side type
	        int j = 0;
		string side_name = cellDefs.get_cell_def(j).get_name();
		while (side_name != "triangle")
		{
		    j++;
		    Insist(j < cellDefs.get_ncell_defs(), 
			   "triangle side type not found for quad pyramid!")
		    side_name = 
		         cellDefs.get_cell_def(j).get_name();
		}
	        ordered_sides[i].resize(3);
		side_types(i) = j;
	    }
	    else
	    {
	        // determine the index for the quad side type
	        int j = 0;
		string side_name = cellDefs.get_cell_def(j).get_name();
		while (side_name != "quad")
		{
		    j++;
		    Insist(j < cellDefs.get_ncell_defs(), 
			   "quad side type not found for quad pyramid!")
		    side_name = 
		        cellDefs.get_cell_def(j).get_name();
		}
	        ordered_sides[i].resize(4);
		side_types(i) = j;
		ordered_sides[i][3] = i%4 + i/2 - 2 * (i/3);
	        sides(i).insert(ordered_sides[i][3]);
	    }
	    ordered_sides[i][0] = i/3 + 2 * (i/4) ;
	    ordered_sides[i][1] = (i + 1)%2 * (1 + i/2 + 2 * (i/4)) + 
	                          (i%2) * (3 + i/3);
	    ordered_sides[i][2] = (i + 1)%2 * (2 + 3 * (i/2) - 4 * (i/4)) +
	                          (i%2) * (4 + i/3);
	    sides(i).insert(ordered_sides[i][0]);
	    sides(i).insert(ordered_sides[i][1]);
	    sides(i).insert(ordered_sides[i][2]);
	    if (debugging)
	    {
	        cout << name << endl;
		cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2];
		if (i != 0 && i != nsides - 1)
		    cout << " " << ordered_sides[i][3];
		cout <<  endl;
	    }
	}

	else if (name == "hexahedron")
	{
	    ordered_sides[i][0] = i/3 + i/4 + 2 * (i/5);
	    ordered_sides[i][1] = ((i+1)%2) * (1 + i/2) + (i%2) * (4 + i/2);
	    ordered_sides[i][2] = ((i+1)%2) * (3 * (1 + i/2) - 2 * (i/4)) +
	                          (i%2) * (5 + 2 * (i/3));
	    ordered_sides[i][3] = ((i+1)%2) * (i + 2) + (i%2) * i;
	    sides(i).insert(ordered_sides[i][0]);
	    sides(i).insert(ordered_sides[i][1]);
	    sides(i).insert(ordered_sides[i][2]);
	    sides(i).insert(ordered_sides[i][3]);
	    if (debugging)
	    {
	        cout << name << endl;
		cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2]
		     << " " << ordered_sides[i][3] << endl;
	    }
	}
    }
}

void RTT_Format::Nodes::readNodes(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::Nodes::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "nodes",
	   "Invalid mesh file: nodes block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::Nodes::readData(ifstream & meshfile)
{
    string dummyString;
    int nodeNum;

    for (int i = 0; i < dims.get_nnodes(); ++i)
    {
	meshfile >> nodeNum;
	Insist(nodeNum == i+1,
	       "Invalid mesh file: node index out of order");
	for (int j = 0; j < dims.get_ndim(); ++j)
	    meshfile >> coords(i,j);
	meshfile >> parents(i);
	--parents(i);
	for (int j = 0; j < dims.get_nnode_flag_types(); ++j)
	{
	    meshfile >> flags(i,j);
	    Insist(nodeFlags.allowed_flag(j, flags(i,j)),
		   "Invalid mesh file: illegal node flag");
	}
	getline(meshfile, dummyString);
    }

    // renumber the nodes to increasing x,y,z coordinate ordering if this 
    // option was selected.
    if (dims.get_renumber())
        sortData();
}

bool RTT_Format::Nodes::sortXYZ(const vector<double> & low_val, 
				  const vector<double> & high_val)
{
    // require agreement to six significant figures for equality. Note that
    // Shawn's tet mesh only agrees to four significant digits.
    const double EPSILON = 1.0e-06;
    vector<double> epsilon(low_val.size());
    bool sorted = true;

    Insist(low_val.size() == high_val.size(),"Improper sort arguments!");
    int dim = low_val.size();
    for (int d = dim - 1; d >= 0; d--)
    {
        if (low_val[d] != 0 && high_val[d] != 0)
	    epsilon[d] = EPSILON * ((fabs(low_val[d]) + fabs(high_val[d]))/2.);
	else
	    epsilon[d] = EPSILON;
        // this strange looking logical operator will sort x,y,(z) coordinates
        // with x varying fastest, followed by y, and lastly z.
        if (high_val[d] < low_val[d] && 
	     fabs(low_val[d] - high_val[d]) > epsilon[d] && 
	       (d == dim-1 || 
	           (d == dim-2 && 
		     fabs(high_val[d+1]-low_val[d+1]) < epsilon[d+1]) ||
	                 (fabs(high_val[d+1]-low_val[d+1]) < epsilon[d+1] &&
		           fabs(high_val[d+2]-low_val[d+2]) < epsilon[d+2] )))
	{
	    sorted = false;
	    d = -1;
	}
    }
    return sorted;
}

void RTT_Format::Nodes::sortData()
{
    vector<vector<double> > sort_vector(dims.get_nnodes(),dims.get_ndim());
    vector<double> original(dims.get_ndim());
    vector<int> temp_parents(dims.get_nnodes());
    vector<vector<int> > temp_flags(dims.get_nnodes(),
				    dims.get_nnode_flag_types());
    sort_map.resize(dims.get_nnodes());

    for (int i = 0; i < dims.get_nnodes(); ++i)
    {
	for (int j = 0; j < dims.get_ndim(); ++j)
	    sort_vector[i][j] = coords(i,j);
	temp_parents[i] = parents(i);
	for (int f = 0; f < dims.get_nnode_flag_types(); ++f)
	    temp_flags[i][f] = flags(i,f);
    }
    sort(sort_vector.begin(),sort_vector.end(),sortXYZ);

    // establish the mapping between the old and new node numbers, and assign
    // the coordinates, parents, and flags with the new numbering.
    for (int i = 0; i < dims.get_nnodes(); ++i)
    {
	for (int d = 0; d < dims.get_ndim(); ++d)
	    original[d] = coords(i,d);
        bool mapped = false;
        int low_index  = 0;
        int high_index = dims.get_nnodes() - 1;
        int k = (high_index + low_index) / 2;
	while (!mapped)
	{
	    // Correlate the original and sorted node numbers/coordinates with
	    // a binary search
            while ((high_index - low_index) > 1)
            {
	        if (sortXYZ(original,sort_vector[k]))
	            high_index = k;
	        else
	            low_index  = k;
                k = (high_index + low_index) / 2;
            }
	    Insist(k < dims.get_nnodes(), "Overflow in node sort while loop!");

	    k = low_index;
	    // Check to make sure these nodes have the same coordinates
	    if (original == sort_vector[k])
	    {
	        mapped = true;
	        sort_map[i] = k;
		for (int d = 0; d < dims.get_ndim(); ++d)
		    coords(i,d) = sort_vector[i][d];
		parents(k) = temp_parents[i];
		for (int f = 0; f < dims.get_nnode_flag_types(); ++f)
		    flags(k,f) = temp_flags[i][f];
	    }
	    else if (k < high_index)
	        ++low_index;
	    else
	        Insist(0, "Error in node sort loop!");
	}
    }
    // finally, now that the nodes are presumably sorted, update the parent
    // nodes to the new numbering system. Note that this step cannot be done
    // in the previous loop because the node and parent numbers are not 
    // necessarily the same.
    for (int i = 0; i < dims.get_nnodes(); ++i)
        parents(i) =  sort_map[parents(i)];

    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;
    if (debugging)
    {
        for (int i = 0; i < dims.get_nnodes(); ++i)
	{
	    cout << i << " " << sort_map[i] << " " << parents(i) << " ";
	    for (int j = 0; j < dims.get_ndim(); ++j)
	        cout << coords(i,j) << " ";
	    for (int f = 0; f < dims.get_nnode_flag_types(); ++f)
	        cout << flags(i,f);
	    cout << endl;
	}
    }
    // free memory
    sort_vector.resize(0);
    original.resize(0);
    temp_parents.resize(0);
    temp_flags.resize(0);
}

void RTT_Format::Nodes::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_nodes",
	   "Invalid mesh file: nodes block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

int RTT_Format::Nodes::get_node(vector<double> node_coords) const
{
    const double EPSILON = 1.0e-06;
    int node_number = 0;
    bool found = false;
    int dim = dims.get_ndim();

    // Find the desired coordinates with a binary search if the nodes are
    // sorted, otherwise we are stuck with a linear search
    if (dims.get_renumber())
    {
        vector<double> local_coords(dim);
        int low_index  = 0;
        int high_index = dims.get_nnodes() - 1;
        node_number = (high_index + low_index) / 2;
        while ((high_index - low_index) > 1)
        {
	    local_coords = get_coords(node_number);
            if (sortXYZ(node_coords,local_coords))
	        high_index = node_number;
	    else
	        low_index  = node_number;
            node_number = (high_index + low_index) / 2;
        }
        Insist(node_number < dims.get_nnodes(), 
	       "Overflow in get_node binary search routine!");
        node_number = low_index;
    }

    // A full linear search if the nodes are not sorted or just a couple
    // of nodes otherwise.
    vector<double> epsilon(dim);
    while (node_number < dims.get_nnodes() && !found)
    {
	int true_count = 0;
        for (int d = 0; d < dim; d++)
	{
            if (node_coords[d] != 0 && get_coords(node_number,d) != 0)
	        epsilon[d] = EPSILON * ((fabs(node_coords[d]) + 
				         fabs(get_coords(node_number,d)))/2.);
	    else
	        epsilon[d] = EPSILON;

	    if (fabs(get_coords(node_number,d) - node_coords[d]) <= epsilon[d])
	        ++true_count;
	}
        if (true_count == dim)
	    found = true;
        else
      	    ++node_number;
    }
    Insist(node_number < dims.get_nnodes(),
     	   "Node number could not be found from its coordinates!");
    return node_number;
}

void RTT_Format::Sides::readSides(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::Sides::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "sides",
	   "Invalid mesh file: sides block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::Sides::readData(ifstream & meshfile)
{
    string dummyString;
    int sideNum;

    for (int i = 0; i < dims.get_nsides(); ++i)
    {
	meshfile >> sideNum;
	Insist(sideNum == i+1,
	       "Invalid mesh file: side index out of order");
	meshfile >> sideType(i);
	--sideType(i);
	Insist(dims.allowed_side_type(sideType(i)),
	       "Invalid mesh file: illegal side type");
	for (int j = 0; j < cellDefs.get_nnodes(sideType(i)); ++j)
	{
	    meshfile >> nodes(i,j);
	    --nodes(i,j);
	}
	for (int j = 0; j < dims.get_nside_flag_types(); ++j)
	{
	    meshfile >> flags(i,j);
	    Insist(sideFlags.allowed_flag(j, flags(i,j)),
		   "Invalid mesh file: illegal side flag");
	}
	getline(meshfile, dummyString);
    }
    // renumber the sides to increasing x,y,z coordinate ordering if this 
    // option was selected.
    if (dims.get_renumber())
        sortData();
}

bool RTT_Format::Sides::sortXYZ(vector<int> low_val, vector<int> high_val)
{
    // compare two integer vectors to determine which contains the lowest 
    // numbers. First sort the vectors into ascending order.
    sort(low_val.begin(),low_val.end());
    sort(high_val.begin(),high_val.end());

    // compare the vectors
    return (low_val < high_val);
}

void RTT_Format::Sides::sortData()
{
    vector<vector<int> > sort_vector(dims.get_nsides(),1);
    vector<int> original(1);
    vector<int> temp_sideType(dims.get_nsides());
    vector<vector<int> > temp_flags(dims.get_nsides(),
				    dims.get_nside_flag_types());

    for (int i = 0; i < dims.get_nsides(); ++i)
    {
        sort_vector[i].resize(cellDefs.get_nnodes(sideType(i)));
	for (int j = 0; j < cellDefs.get_nnodes(sideType(i)); ++j)
	{
	    // map the user-input node numbers to the sorted node numbers.
	    nodes(i,j) = ptrNodes.get_map(nodes(i,j));
	    sort_vector[i][j] = nodes(i,j);
	}

	temp_sideType[i] = sideType(i);

	for (int f = 0; f < dims.get_nside_flag_types(); ++f)
	    temp_flags[i][f] = flags(i,f);
    }
    sort(sort_vector.begin(),sort_vector.end(),sortXYZ);

    // establish the mapping between the old and new side numbers, and assign
    // the nodes, sidetypes, and flags with the new numbering.
    for (int i = 0; i < dims.get_nsides(); ++i)
    {
        if (original.size() != cellDefs.get_nnodes(sideType(i)))
	    original.resize(cellDefs.get_nnodes(sideType(i)));
	for (int j = 0; j < cellDefs.get_nnodes(sideType(i)); ++j)
	    original[j] = nodes(i,j);

        bool mapped = false;
        int low_index  = 0;
        int high_index = dims.get_nsides() - 1;
        int k = (high_index + low_index) / 2;
	while (!mapped)
	{
	    // Correlate the original and sorted side numbers with a binary 
	    // search
            while ((high_index - low_index) > 1)
            {
	        if (sortXYZ(original,sort_vector[k]))
	            high_index = k;
	        else
	            low_index  = k;
                k = (high_index + low_index) / 2;
            }
	    Insist(k < dims.get_nsides(), "Overflow in side sort while loop!");

	    k = low_index;
	    // Check to make sure these sides have the same nodes
	    if (original == sort_vector[k])
	    {
	        mapped = true;
		for (int j = 0; j < cellDefs.get_nnodes(sideType(i)); ++j)
		    nodes(i,j) = sort_vector[i][j];
		sideType(k) = temp_sideType[i];
		for (int f = 0; f < dims.get_nside_flag_types(); ++f)
		    flags(k,f) = temp_flags[i][f];
	    }
	    else if (k < high_index)
	        ++low_index;
	    else
	        Insist(0, "Error in side sort loop!");
	}
    }

    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;

    if (debugging)
    {
        for (int i = 0; i < dims.get_nsides(); ++i)
	{
	    cout << i << " " << sideType(i) << " ";
	    for (int j = 0; j < cellDefs.get_nnodes(sideType(i)); ++j)
	        cout << nodes(i,j) << " ";
	    for (int f = 0; f < dims.get_nside_flag_types(); ++f)
	        cout << flags(i,f);
	    cout << endl;
	}
    }
    // free memory
    sort_vector.resize(0);
    original.resize(0);
    temp_sideType.resize(0);
    temp_flags.resize(0);
}

void RTT_Format::Sides::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_sides",
	   "Invalid mesh file: sides block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::Cells::readCells(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}

void RTT_Format::Cells::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "cells",
	   "Invalid mesh file: cells block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::Cells::readData(ifstream & meshfile)
{
    string dummyString;
    int cellNum;

    for (int i = 0; i < dims.get_ncells(); ++i)
    {
	meshfile >> cellNum;
	Insist(cellNum == i+1,
	       "Invalid mesh file: cell index out of order");
	meshfile >> cellType(i);
	--cellType(i);
	Insist(dims.allowed_cell_type(cellType(i)),
	       "Invalid mesh file: illegal cell type");
	for (int j = 0; j < cellDefs.get_nnodes(cellType(i)); ++j)
	{
	    meshfile >> nodes(i,j);
	    --nodes(i,j);
	}
	for (int j = 0; j < dims.get_ncell_flag_types(); ++j)
	{
	    meshfile >> flags(i,j);
	    Insist(cellFlags.allowed_flag(j, flags(i,j)),
		   "Invalid mesh file: illegal cell flag");
	}
	getline(meshfile, dummyString);
    }
    // renumber the cells to increasing x,y,z coordinate ordering if this 
    // option was selected.
    if (dims.get_renumber())
        sortData();
}

void RTT_Format::Cells::sortData()
{
    vector<vector<int> > sort_vector(dims.get_ncells(),1);
    vector<int> original(1);
    vector<int> temp_cellType(dims.get_ncells());
    vector<vector<int> > temp_flags(dims.get_ncells(),
				    dims.get_ncell_flag_types());

    for (int i = 0; i < dims.get_ncells(); ++i)
    {
        sort_vector[i].resize(cellDefs.get_nnodes(cellType(i)));
	for (int j = 0; j < cellDefs.get_nnodes(cellType(i)); ++j)
	{
	    // map the user-input node numbers to the sorted node numbers.
	    nodes(i,j) = ptrNodes.get_map(nodes(i,j));
	    sort_vector[i][j] = nodes(i,j);
	}
	sort(sort_vector[i].begin(),sort_vector[i].end());

	temp_cellType[i] = cellType(i);

	for (int f = 0; f < dims.get_ncell_flag_types(); ++f)
	    temp_flags[i][f] = flags(i,f);
    }
    sort(sort_vector.begin(),sort_vector.end());

    // establish the mapping between the old and new cell numbers, and assign
    // the nodes, celltypes, and flags with the new numbering.
    for (int i = 0; i < dims.get_ncells(); ++i)
    {
        if (original.size() != cellDefs.get_nnodes(cellType(i)))
	    original.resize(cellDefs.get_nnodes(cellType(i)));
	for (int j = 0; j < cellDefs.get_nnodes(cellType(i)); ++j)
	    original[j] = nodes(i,j);
	sort(original.begin(),original.end());

        bool mapped = false;
        int low_index  = 0;
        int high_index = dims.get_ncells() - 1;
        int k = (high_index + low_index) / 2;
	while (!mapped)
	{
	    // Correlate the original and sorted cell numbers with a binary 
	    // search
            while ((high_index - low_index) > 1)
            {
	        if (original < sort_vector[k])
	            high_index = k;
	        else
	            low_index  = k;
                k = (high_index + low_index) / 2;
            }
	    Insist(k < dims.get_ncells(), "Overflow in cell sort while loop!");

	    k = low_index;
	    // Check to make sure these cells have the same nodes
	    if (original == sort_vector[k])
	    {
	        mapped = true;
		for (int j = 0; j < cellDefs.get_nnodes(cellType(i)); ++j)
		    nodes(i,j) = sort_vector[i][j];
		cellType(k) = temp_cellType[i];
		for (int f = 0; f < dims.get_ncell_flag_types(); ++f)
		    flags(k,f) = temp_flags[i][f];
	    }
	    else if (k < high_index)
	        ++low_index;
	    else
	        Insist(0, "Error in cell sort loop!");
	}
    }

    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;

    if (debugging)
    {
        for (int i = 0; i < dims.get_ncells(); ++i)
	{
	    cout << i << " " << cellType(i) << " ";
	    for (int j = 0; j < cellDefs.get_nnodes(cellType(i)); ++j)
	        cout << nodes(i,j) << " ";
	    for (int f = 0; f < dims.get_ncell_flag_types(); ++f)
	        cout << flags(i,f);
	    cout << endl;
	}
    }
    // free memory
    sort_vector.resize(0);
    original.resize(0);
    temp_cellType.resize(0);
    temp_flags.resize(0);
}

void RTT_Format::Cells::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_cells",
	   "Invalid mesh file: cells block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::NodeData::readNodeData(ifstream & meshfile)
{
    readKeyword(meshfile);
    if (dims.get_nnode_data() > 0)
    {
        if (!dims.get_renumber())
	    readData(meshfile);
        else
	    sortData(meshfile);
    }
    readEndKeyword(meshfile);
}

void RTT_Format::NodeData::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "nodedat",
	   "Invalid mesh file: nodedat block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::NodeData::readData(ifstream & meshfile)
{
    string dummyString;
    int nodeNum;

    for (int i = 0; i < dims.get_nnodes(); ++i)
    {
	meshfile >> nodeNum;
	Insist(nodeNum == i+1,
	       "Invalid mesh file: node data index out of order");
	for (int j = 0; j < dims.get_nnode_data(); ++j)
	    meshfile >> data(i,j);
	getline(meshfile, dummyString);
    }
}

void RTT_Format::NodeData::sortData(ifstream & meshfile)
{
    string dummyString;
    int nodeNum;

    for (int i = 0; i < dims.get_nnodes(); ++i)
    {
	meshfile >> nodeNum;
	Insist(nodeNum == i+1,
	       "Invalid mesh file: node data index out of order");
	for (int j = 0; j < dims.get_nnode_data(); ++j)
	    meshfile >> data(ptrNodes.get_map(i),j);
	getline(meshfile, dummyString);
    }
}

void RTT_Format::NodeData::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_nodedat",
	   "Invalid mesh file: nodedat block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::SideData::readSideData(ifstream & meshfile)
{
    readKeyword(meshfile);
    if (dims.get_nside_data() > 0)
    {
        if (!dims.get_renumber())
	    readData(meshfile);
        else
	    sortData(meshfile);
    }
    readEndKeyword(meshfile);
}

void RTT_Format::SideData::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "sidedat",
	   "Invalid mesh file: sidedat block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::SideData::readData(ifstream & meshfile)
{
    string dummyString;
    int sideNum;

    for (int i = 0; i < dims.get_nsides(); ++i)
    {
	meshfile >> sideNum;
	Insist(sideNum == i+1,
	       "Invalid mesh file: side data index out of order");
	for (int j = 0; j < dims.get_nside_data(); ++j)
	    meshfile >> data(i,j);
	getline(meshfile, dummyString);
    }
}

void RTT_Format::SideData::sortData(ifstream & meshfile)
{
    string dummyString;
    int sideNum;

    for (int i = 0; i < dims.get_nsides(); ++i)
    {
	meshfile >> sideNum;
	Insist(sideNum == i+1,
	       "Invalid mesh file: side data index out of order");
	for (int j = 0; j < dims.get_nside_data(); ++j)
	    meshfile >> data(ptrNodes.get_map(i),j);
	getline(meshfile, dummyString);
    }
}

void RTT_Format::SideData::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_sidedat",
	   "Invalid mesh file: sidedat block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

void RTT_Format::CellData::readCellData(ifstream & meshfile)
{
    readKeyword(meshfile);
    if (dims.get_ncell_data() > 0)
    {
        if (!dims.get_renumber())
	    readData(meshfile);
        else
	    sortData(meshfile);
    }
    readEndKeyword(meshfile);
}

void RTT_Format::CellData::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "celldat",
	   "Invalid mesh file: celldat block missing");
    getline(meshfile, dummyString);
}

void RTT_Format::CellData::readData(ifstream & meshfile)
{
    string dummyString;
    int cellNum;

    for (int i = 0; i < dims.get_ncells(); ++i)
    {
	meshfile >> cellNum;
	Insist(cellNum == i+1,
	       "Invalid mesh file: cell data index out of order");
	for (int j = 0; j < dims.get_ncell_data(); ++j)
	    meshfile >> data(i,j);
	getline(meshfile, dummyString);
    }
}

void RTT_Format::CellData::sortData(ifstream & meshfile)
{
    string dummyString;
    int cellNum;

    for (int i = 0; i < dims.get_ncells(); ++i)
    {
	meshfile >> cellNum;
	Insist(cellNum == i+1,
	       "Invalid mesh file: cell data index out of order");
	for (int j = 0; j < dims.get_ncell_data(); ++j)
	    meshfile >> data(ptrNodes.get_map(i),j);
	getline(meshfile, dummyString);
    }
}

void RTT_Format::CellData::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_celldat",
	   "Invalid mesh file: celldat block missing end");
    getline(meshfile, dummyString);            // read and discard blank line.
}

RTT_Format::Connectivity::Connectivity(const Dims & dims_, 
    const CellDefs & cellDefs_, const Cells & cells_, const Sides & sides_,
    const Nodes & nodes_)
    : dims(dims_), cellDefs(cellDefs_), cells(cells_),sides(sides_),
      nodes(nodes_), adjCell(dims.get_ncells(),
      vector<vector<int> >(dims.get_nsides_max(),vector<int> (1)))
{
    calcAdjacentCells();
}

void RTT_Format::Connectivity::calcAdjacentCells()
{
    multimap<set<int>, int> faceCells;
    multimap<set<int>, int> sideCells;
    multimap<vector<int>, vector<int> > unstCells;
    multimap<double, vector<int>, greater<double> > circum;
    multimap<int, vector<int> > unstNodes;
    vector<int> cell_faces(2);

    // create a multimap of the node sets that define the cell sides so that 
    // the nodes are in ascending order. Pretty smart Shawn!
    for (int c = 0; c < dims.get_ncells(); ++c)
    {
        // But we have to decide the cell type to determine the number of
        // nodes used in the cell definition first for the general case.
        int cellType = cells.get_type(c);
	int nsides = cellDefs.get_nsides(cellType);
	if (adjCell[c].size() != nsides)
	    adjCell[c].resize(nsides);
	const CellDef & cellDef = cellDefs.get_cell_def(cellType);

	for (int s = 0; s < nsides ; ++s)
	{
	    set<int> face;
	    const set<int> & faceNodes = cellDef.get_side(s);
	    for (set<int>::const_iterator mmiter = faceNodes.begin();
		 mmiter != faceNodes.end(); ++mmiter)
	    {
	         int node_index = * mmiter;
	         face.insert(cells.get_nodes(c,node_index));
	    }
	    faceCells.insert(make_pair(face, c));
	}
    }

    // Do the same for the sides (but all of the side nodes are included in
    // the sets).
    for (int s = 0; s < dims.get_nsides(); ++s)
    {
        // Decide the side type to determine the number of  nodes used in the 
        // side definition.
	int nnodes = cellDefs.get_nnodes(sides.get_type(s));
	set<int> face;

	for (int node_index = 0; node_index < nnodes ; ++node_index)
	    face.insert(sides.get_nodes(s,node_index));

	sideCells.insert(make_pair(face, s));
    }
    // find the side flag number that contains the boundary conditions
    int boundary = sides.get_boundary_flag_number();

    // create a separate iterator for the sides.
    multimap<set<int>, int>::iterator sitr = sideCells.begin();

    // loop over all of the cell faces. Iteration index is also incremented
    // once within the loop.
    for (multimap<set<int>, int>::iterator mmiter = faceCells.begin();
	 mmiter != faceCells.end(); ++mmiter)
    {
	const set<int> & globalFaceNodes = mmiter->first;
	int cell = mmiter->second;
	int otherCell;
	set<int> localFaceNodes;
	set<int> otherLocalFaceNodes;
	int cellType = cells.get_type(cell);
	int nnodes = cellDefs.get_nnodes(cellType);
	int nsides = cellDefs.get_nsides(cellType);
	int faceNum;
	int otherFaceNum;

	// Determine which of this cell's nodes are tied to this face and put
	// this (normalized) info into the localFaceNodes set.
	for (int n = 0; n < nnodes; ++n)
	{
	    int globalNode = cells.get_nodes(cell,n);
	    if (globalFaceNodes.count(globalNode) == 1)
		localFaceNodes.insert(n);
	}

	// Correlate the localFaceNodes set with the user input face numbering
	// scheme and assign the face number to faceNum.
	const CellDef & cellDef = cellDefs.get_cell_def(cellType);
	for (int f = 0; f < nsides; ++f)
	{
	    const set<int> & faceNodes = cellDef.get_side(f);
	    if (localFaceNodes == faceNodes)
		faceNum = f;
	}

	// If the face only appears once in the definitions of the face cells
	// we must be on a boundary. For a structured mesh these faces are the
	// sides, while for an unstructured mesh these faces can also be next
	// to an irregular node. We can use this fact to resolve the refinement
	// level of a continuous adaptive mesh. Since the sides are sorted in 
	// the same order as the cell faces, we only have to look at the next 
	// element in the side multimap.
	if (faceCells.count(globalFaceNodes) == 1)
	{
	    // given my current hind-sight, having a complete list of all
	    // of the cell faces that are either a boundary between cells
	    // with different generation levels or on the physical problem
	    // boundary will be really useful for assigning generations 
	    // later (since I ended up having to reconstruct this exact 
	    // data outside this routine with yet another linear search of
	    // every cell in the mesh). Store the pairs with the face as 
	    // the multimap key so that the faces will be grouped according
	    // to direction.
	    bndryFaces.insert(make_pair(faceNum,cell));

	    const set<int> & globalSideNodes = sitr->first;
	    int side = sitr->second;
	    // if the nodes that comprise this face also occur in the nodes
	    // that define the sides, assign the negative of the boundary side 
	    // flag as the adjacent cell.
	    if (globalSideNodes == globalFaceNodes)
	    {
		adjCell[cell][faceNum][0] = - sides.get_flags(side,boundary);
		++sitr;
	    }
	    // this face is the boundary between two unstructured cells. 
	    else
	    {
		int faceType = cellDef.get_side_types(faceNum);
		const CellDef & faceDef = cellDefs.get_cell_def(faceType);
		vector<int> unstGlobalNodes(faceDef.get_nnodes());
		set<int> lines;
		cell_faces[0] = cell;
		cell_faces[1] = faceNum;

		// create an ordered set of the cell face nodes to calculate
		// the face circumference length. Add all of these nodes to
		// the unstNodes multimap.
		const vector<int> & unstFaceNodes = 
		    cellDef.get_ordered_side(faceNum);
		for (int n = 0; n < faceDef.get_nnodes(); ++n)
		{
		    unstGlobalNodes[n] = 
		        cells.get_nodes(cell,unstFaceNodes[n]);
		    unstNodes.insert(make_pair(unstGlobalNodes[n],cell_faces));
		}
		// create a multimap of the unconnected cells, faces, and 
		// nodes.
	        unstCells.insert(make_pair(cell_faces,unstGlobalNodes));

		// build the face lines from the "side cell type" definition.
		// calculate the face circumference length and create a 
		// multimap of this information.
		double length = 0.;
		for (int s = 0; s < faceDef.get_nsides(); ++s)
		{
		    set<int> lineNodes = faceDef.get_side(s);

		    set<int>::iterator bitr = lineNodes.begin();
		    set<int>::reverse_iterator eitr = lineNodes.rbegin();
		    double line_length = 0;
		    for (int d = 0; d < dims.get_ndim(); d++)
		    {
		        double diff = 
			    nodes.get_coords(unstGlobalNodes[* bitr],d) -
			    nodes.get_coords(unstGlobalNodes[* eitr],d);
		            line_length += pow(diff,2.0);
		    }
		    length += sqrt(line_length);
		}
		circum.insert(make_pair(length,cell_faces));
	    }
	}
	// this face is the junction to another structured cell (the next 
	// face in the multimap).
	else
	{
	    // temporary work-around to ICEM bug - 16 Aug 99.
	    const set<int> & globalSideNodes = sitr->first;
	    int side = sitr->second;
	    if (globalSideNodes == globalFaceNodes)
	    {
	        ++mmiter;
		otherCell = mmiter->second;	        
		cout << "Warning: Face shared by cells " << cell + 1 << " & "
		     << otherCell + 1 << " also exists as boundary side "
		     << side + 1 << endl;
		++sitr;
	        --mmiter;
	    }
	    ++mmiter;
	    const set<int> & globalFaceNodes = mmiter->first;
	    otherCell = mmiter->second;
	    int otherCellType = cells.get_type(otherCell);
	    int otherNnodes = cellDefs.get_nnodes(otherCellType);
	    int otherNsides = cellDefs.get_nsides(otherCellType);
	    for (int n = 0; n < otherNnodes; ++n)
	    {
	        // Determine if the next cell shares these nodes and store
	        // the (normalized) data set in otherLocalFaceNodes.
		int globalNode = cells.get_nodes(otherCell,n);
		if (globalFaceNodes.count(globalNode) == 1)
		    otherLocalFaceNodes.insert(n);
	    }
	    // Correlate the otherlocalFaceNodes set with the user input face
	    // numbering scheme and assign to otherfaceNum.
	    for (int f = 0; f < otherNsides; ++f)
	    {
		const set<int> & faceNodes = cellDef.get_side(f);
		if (otherLocalFaceNodes == faceNodes)
		    otherFaceNum = f;
	    }

	    adjCell[cell][faceNum][0] = otherCell;
	    adjCell[otherCell][otherFaceNum][0] = cell;
	}
    }
    // make sure all of the sides were correctly identified and assigned to 
    // the corresponding cell faces, then clear memory.
    Insist(sitr == sideCells.end(),"Side/Cell Face correspondance not found!");
    sideCells.clear();
    faceCells.clear();

    // Done with all of the structured cells. This essentially leaves a three-
    // dimensional jigsaw puzzle that is composed of blocks formed by cells 
    // that share a common geometry. The unassigned faces represent the
    // interlocks between these various pieces. Some of these faces will be 
    // "whole" in that groups of cells need to be combined into "subfaces" to
    // form the corresponding mating face. The whole faces can be distinguished
    // from the subfaces by size. The circum multimap contains the cell face 
    // circumferences in descending order. For the case of a continuous AMR
    // mesh, the center of each face corresponds to the location of a single 
    // node that is common to each of the subgrouped cells (and this node will
    // not occur anywhere else). Use the data in the unstCells  multimap to 
    // calculate the coordinates of the center of the face, find the node
    // that lives there, and assign each of the cell faces that contain this
    // node (as identified in the unstNode multimap) to this whole face. Then
    // remove all of these cell faces from the unstCells multimap. Note that 
    // the trick of locating the grouped subfaces by finding the whole face
    // center is unique to the continuous adaptive refinement mesh. Another
    // technique would probably be required for another mesh topology. The fact
    // that these "junction nodes" appear less frequently (because they do not
    // actually exist on the whole face) might be of use. The fact that the
    // "junction lines" between two grouped subfaces appear exactly twice could
    // also be valuable.

    multimap<double, vector<int>, greater<double> >::iterator cfiter = 
        circum.begin();

    // this is a no-op if the mesh is fully structured - a countdown otherwise.
    while (!unstCells.empty())
    {
        // skip over grouped cell faces in circum that have already been 
        // assigned to whole faces and erased from the unstCells multimap.
	while (unstCells.count(cfiter->second) == 0)
	    ++cfiter;

	cell_faces = cfiter->second;	    
	int cell = cell_faces[0];
	int faceNum = cell_faces[1];

	// calculate the coordinates of the face center point.
	vector<int> & faceNodes = unstCells.find(cell_faces)->second;
	vector<double> center(dims.get_ndim());
	for (int n = 0; n < faceNodes.size(); n++)
	{
	    for (int d = 0; d < dims.get_ndim(); d++)
	        center[d] += nodes.get_coords(faceNodes[n],d);
	}
	for (int d = 0; d < dims.get_ndim(); d++)
	    center[d] /= static_cast<double>(faceNodes.size());

	// find a node with these coordinates, and set up a corresponding 
	// iterator in the unstNodes multimap.
	int centerNode = nodes.get_node(center);
	multimap<int,vector<int> >::iterator fiter=unstNodes.find(centerNode);

	// remove the null value from the end of adjCell vector.
	if (!adjCell[cell][faceNum].empty())
	    adjCell[cell][faceNum].pop_back();

	// assign all of the cell faces that contain this node to this whole
	// face, and assign this whole face to each of the subfaces.
	while (fiter !=  unstNodes.end() && fiter->first == centerNode)
	{
	    vector<int> other_cell_faces = fiter->second;
	    int otherCell = other_cell_faces[0];
	    int otherFaceNum = other_cell_faces[1];

	    adjCell[cell][faceNum].push_back(otherCell);
	    adjCell[otherCell][otherFaceNum][0] = cell;
	    unstCells.erase(other_cell_faces);
	    ++fiter;
	}
	unstCells.erase(cell_faces);
	++cfiter;
    }
    // clean house.
    circum.clear();
    unstNodes.clear();
    cell_faces.resize(0);

    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;

    if (debugging)
    {
        for (int c = 0; c < dims.get_ncells(); ++c)
	{
	    int cellType = cells.get_type(c);
	    int nsides = cellDefs.get_nsides(cellType);
	    for (int f = 0; f < nsides; f++)
	    {
	        cout << "cell " << c << " face " << f << " adjacent cell(s) ";
	        for (int af = 0; af < adjCell[c][f].size(); af++)
	            cout << get_adjCell(c,f,af) << " ";
		cout << endl;
	    }
	}
    }
}

} // end namespace rtt_format

//---------------------------------------------------------------------------//
//                              end of RTT_Format.cc
//---------------------------------------------------------------------------//
