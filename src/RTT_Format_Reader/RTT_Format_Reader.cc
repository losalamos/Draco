//----------------------------------*-C++-*--------------------------------//
// RTT_Format_Reader.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/RTT_Format_Reader.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader library.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "RTT_Format_Reader.hh"
#include "ds++/Assert.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Constructs an RTT_Format_Reader object and parses the mesh data.
 * \param RTT_File Mesh file name.
 * \param renumber Turns the option to reassign the node, side, and cell 
 *        numbers based upon coordinates in ascending order (x, y, and then z)
 *        on and off (default is no renumbering).
 */
RTT_Format_Reader::RTT_Format_Reader(const string & RTT_File, 
				     const bool & renumber)
{
    readMesh(RTT_File,renumber);
    calculateConnectivity();
}
/*!
 * \brief Instantiates a Connectivity and determines the mesh connectivity.
 */
void RTT_Format_Reader::calculateConnectivity()
{
    spConnectivity = new Connectivity(dims, * spCellDefs, * spCells, 
				      * spSides, * spNodes);
}
/*!
 * \brief Parses the mesh file data via calls to the member data class objects
 *        public "read" member functions.
 * \param RTT_File Mesh file name.
 * \param renumber Turns the option to reassign the node, side, and cell 
 *        numbers based upon coordinates in ascending order (x, y, and then z)
 *        on and off (defaults to no renumbering).
 */
void RTT_Format_Reader::readMesh(const string & RTT_File, 
				 const bool & renumber)
{
    const char * file = RTT_File.c_str();
    ifstream meshfile(file, std::ios::in);
    if (!meshfile)
    {
        std::cout << "File could not be opened\n";
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
    catch (rtt_dsxx::assertion as)
    {
        std::cout << "Assertion thrown: " << as.what() << std::endl;
        Insist(false, as.what());
    }

    if (renumber)
    {
        spCellDefs->sortData();
        spNodes->sortData();
        spSides->sortData();
        spCells->sortData();
        spNodeData->sortData();
        spSideData->sortData();
        spCellData->sortData();
    }
}
/*!
 * \brief Reads and validates the magic cookie at the beginning of the mesh 
 *        file.
 * \param meshfile Mesh file name.
 */
void RTT_Format_Reader::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "rtt_ascii", "Invalid mesh file: Not an RTT file");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}

/*!
 * \brief Instantiates the majority of the RTT_Format_Reader member data class
 *        objects.
 */
void RTT_Format_Reader::createMembers()
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
    spSideData = new SideData(dims, * spSides);
    spCellData = new CellData(dims, * spCells);
}

/*!
 * \brief Reads the node, side, and cell flag blocks from the mesh file.
 * \param meshfile Mesh file name.
 */
void RTT_Format_Reader::readFlagBlocks(ifstream & meshfile)
{
    spNodeFlags->readNodeFlags(meshfile);
    spSideFlags->readSideFlags(meshfile);
    spCellFlags->readCellFlags(meshfile);
}
/*!
 * \brief Reads the node, side, and cell data_id blocks from the mesh file.
 * \param meshfile Mesh file name.
 */
void RTT_Format_Reader::readDataIDs(ifstream & meshfile)
{
    spNodeDataIds->readDataIDs(meshfile);
    spSideDataIds->readDataIDs(meshfile);
    spCellDataIds->readDataIDs(meshfile);
}
/*!
 * \brief Reads and validates the end_rtt_mesh keyword at the end of the mesh 
 *        file.
 * \param meshfile Mesh file name.
 */
void RTT_Format_Reader::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_rtt_mesh",
	   "Invalid mesh file: RTT file missing end");
    std::getline(meshfile, dummyString);
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                            end of RTT_Format_Reader.cc
//---------------------------------------------------------------------------//
