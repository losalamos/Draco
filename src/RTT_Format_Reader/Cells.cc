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
/*!
 * \brief Renumbers the cells in ascending order based upon their node numbers.
 */
void Cells::sortData()
{
    vector_vector_int sort_vector(dims.get_ncells(),1);
    vector_int original(1);
    vector_int temp_cellType = cellType;
    vector_vector_int temp_flags = flags;
    sort_map.resize(dims.get_ncells());

    for (int i = 0; i < dims.get_ncells(); ++i)
    {
        sort_vector[i].resize(cellDefs.get_nnodes(cellType[i]));
	for (int j = 0; j < cellDefs.get_nnodes(cellType[i]); ++j)
	    // map the user-input node numbers to the sorted node numbers.
	    nodes[i][j] = nodesClass.get_map(nodes[i][j]);
	sort_vector[i] = nodes[i];
	std::sort(sort_vector[i].begin(),sort_vector[i].end());
    }
    std::sort(sort_vector.begin(),sort_vector.end());

    // establish the mapping between the old and new cell numbers, and assign
    // the nodes, celltypes, and flags with the new numbering.
    for (int i = 0; i < dims.get_ncells(); ++i)
    {
        if (original.size() != cellDefs.get_nnodes(cellType[i]))
	    original.resize(cellDefs.get_nnodes(cellType[i]));
	original = nodes[i];
	std::sort(original.begin(),original.end());

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
		sort_map[i] = k;
		nodes[i] = sort_vector[i];
		cellType[k] = temp_cellType[i];
		flags[k] = temp_flags[i];
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
	    std::cout << i << " " << sort_map[i] << " " << cellType[i] << " ";
	    for (int j = 0; j < cellDefs.get_nnodes(cellType[i]); ++j)
	        std::cout << nodes[i][j] << " ";
	    for (int f = 0; f < dims.get_ncell_flag_types(); ++f)
		std::cout << flags[i][f];
	    std::cout << std::endl;
	}
    }
    // free memory
    sort_vector.resize(0);
    original.resize(0);
    temp_cellType.resize(0);
    temp_flags.resize(0);
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                    end of RTT_Format_Reader/Cells.cc
//---------------------------------------------------------------------------//
