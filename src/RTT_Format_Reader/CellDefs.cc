//----------------------------------*-C++-*--------------------------------//
// CellDefs.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/CellDefs.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/CellDefs class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "CellDefs.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Parses the cell_defs (cell definitions) data block from the mesh 
 *        file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void CellDefs::readCellDefs(ifstream & meshfile)
{
    readKeyword(meshfile);
    readDefs(meshfile);
    readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the cell_defs block (cell definitions) keyword.
 * \param meshfile Mesh file name.
 */
void CellDefs::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "cell_defs",
	   "Invalid mesh file: cell_defs block missing");
    std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the cell_defs (cell definitions) block data.
 * \param meshfile Mesh file name.
 */
void CellDefs::readDefs(ifstream & meshfile)
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
	defs[i] = new CellDef(*this, dummyString);
	std::getline(meshfile, dummyString);
	defs[i]->readDef(meshfile);
    }
}
/*!
 * \brief Reads and validates the end_cell_defs block keyword.
 * \param meshfile Mesh file name.
 */
void CellDefs::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_cell_defs",
	   "Invalid mesh file: cell_defs block missing end");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}
/*!
 * \brief Used by the CellDefs class objects to parse the number of nodes and
 *        sides per cell, the side type indices, and the nodes for each side.
 * \param meshfile Mesh file name.
 */
void CellDef::readDef(ifstream & meshfile)
{
    string dummyString;

    meshfile >> nnodes >> nsides;
    side_types.resize(nsides);
    sides.resize(nsides);
    ordered_sides.resize(nsides);
    std::getline(meshfile, dummyString);

    for (int i = 0; i < nsides; ++i)
    {
	meshfile >> side_types[i];
	--side_types[i];
    }
    if (nsides > 0)
	std::getline(meshfile, dummyString);

    // note that this implementation does not preserve the "right hand rule"
    // of the cell definitions due to the use of a set container (which is 
    // sorted). It is slicker than snail snot when it comes time to implement 
    // the connectivity, however. The ordered_sides vector was added to allow
    // the original ordered data to be retained.
    int side;
    for (int i = 0; i < nsides; ++i)
    {
        int numb_nodes = cellDefs.get_cell_def(side_types[i]).get_nnodes();
        ordered_sides[i].resize(numb_nodes);
	for (int j = 0; j < numb_nodes; ++j)
	{
	    meshfile >> side;
	    --side;
	    sides[i].insert(side);
	    ordered_sides[i][j] = side;
	}
	if (sides[i].size() > 0)
	    std::getline(meshfile, dummyString);
    }
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                    end of RTT_Format_Reader/CellDefs.cc
//---------------------------------------------------------------------------//
