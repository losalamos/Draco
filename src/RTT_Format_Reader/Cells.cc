//----------------------------------*-C++-*--------------------------------//
// Cells.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/Cells.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Cells class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Cells.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Parses the cells block data from the mesh file via calls to 
 *        private member functions.
 * \param meshfile Mesh file name.
 */
void Cells::readCells(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the cells block keyword.
 * \param meshfile Mesh file name.
 */
void Cells::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "cells", "Invalid mesh file: cells block missing");
    std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the cells block data.
 * \param meshfile Mesh file name.
 */
void Cells::readData(ifstream & meshfile)
{
    string dummyString;
    int cellNum;

    for (int i = 0; i < dims.get_ncells(); ++i)
    {
	meshfile >> cellNum;
	Insist(cellNum == i+1, "Invalid mesh file: cell index out of order");
	meshfile >> cellType[i];
	--cellType[i];
	Insist(dims.allowed_cell_type(cellType[i]),
	       "Invalid mesh file: illegal cell type");
	nodes[i].resize(cellDefs.get_nnodes(cellType[i]));
	for (int j = 0; j < cellDefs.get_nnodes(cellType[i]); ++j)
	{
	    meshfile >> nodes[i][j];
	    --nodes[i][j];
	}
	for (int j = 0; j < dims.get_ncell_flag_types(); ++j)
	{
	    meshfile >> flags[i][j];
	    Insist(cellFlags.allowed_flag(j, flags[i][j]),
		   "Invalid mesh file: illegal cell flag");
	}
	std::getline(meshfile, dummyString);
    }
}
/*!
 * \brief Reads and validates the end_cells block keyword.
 * \param meshfile Mesh file name.
 */
void Cells::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_cells",
	   "Invalid mesh file: cells block missing end");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                    end of RTT_Format_Reader/Cells.cc
//---------------------------------------------------------------------------//
