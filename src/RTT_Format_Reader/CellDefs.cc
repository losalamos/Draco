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
 * \brief Redefines the cells with the nodes and sides numbered in ascending 
 *        order based upon their coordinates (x, y, and then z).
 * \param meshfile Mesh file name.
 */
void CellDefs::sortData()
{
    for (int i = 0; i < dims.get_ncell_defs(); ++i)
    {
	defs[i]->sortData();
    }
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
/*!
 * \brief Redefines the cells with the nodes and sides numbered in ascending 
 *        order based upon their coordinates (x, y, and then z). 
 *
 *        This cell definition scheme has the distinct advantage that the 
 *        cell nodes are inherently in increasing integer order. This 
 *        eliminates the need to derive a coordinate transform from the 
 *        user-input cell definitions.
 */
void CellDef::sortData()
{
    // The cell definitions built into ICEM do not correspond to the AMR mesh
    // node numbering scheme. Define the alternative cell definitions herein, 
    // while retaining the right hand rule of cell definition. This section 
    // takes heavy advantage of the C++ rules of truncating integer division.
    // Note that our cell definition scheme has the distinct advantage that
    // the cell nodes are inherently in increasing integer order. This 
    // eliminates any need to derive a coordinate transform from the 
    // user-input cell definitions.  

    ordered_sides.resize(nsides);
    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;

    for (int i = 0; i < nsides; ++i)
    {

        int numb_nodes = cellDefs.get_cell_def(side_types[i]).get_nnodes();
	sides[i].erase(sides[i].begin(),sides[i].end());
        ordered_sides[i].resize(numb_nodes);

	if (name == "line")
        {
	    ordered_sides[i][0] = i;
	    sides[i].insert(ordered_sides[i][0]);
	    if (debugging)
	    {
	        std::cout << name << std::endl;
		std::cout << i << " " << ordered_sides[i][0] << std::endl;
	    }
	}

	else if (name == "triangle")
	{
	    ordered_sides[i][0] = (i + 1)%2 + i/2;
	    ordered_sides[i][1] = (3 - i)%3;
	    sides[i].insert(ordered_sides[i][0]);
	    sides[i].insert(ordered_sides[i][1]);
	    if (debugging)
	    {
	        std::cout << name << std::endl;
		std::cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1] << std::endl;
	    }
	}

	else if (name == "quad")
	{
	    ordered_sides[i][0] = (i + 1)%2 + 2 * (i/2);
	    ordered_sides[i][1] = i%2 + (i + 1)/2;
	    sides[i].insert(ordered_sides[i][0]);
	    sides[i].insert(ordered_sides[i][1]);
	    if (debugging)
	    {
	        std::cout << name << std::endl;
		std::cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1] << std::endl;
	    }
	}

	else if (name == "tetrahedron")
	{
	    ordered_sides[i][0] = i/3;
	    ordered_sides[i][1] = 1 + (3 - i)%3 + 2 * (i/3);
	    ordered_sides[i][2] = (2 - i)%3 + 3 * (i/2);
	    sides[i].insert(ordered_sides[i][0]);
	    sides[i].insert(ordered_sides[i][1]);
	    sides[i].insert(ordered_sides[i][2]);
	    if (debugging)
	    {
	        std::cout << name << std::endl;
		std::cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2] << std::endl;
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
		side_types[i] = j;
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
		side_types[i] = j;
		ordered_sides[i][3] = 2;
	        sides[i].insert(ordered_sides[i][3]);
	    }
	    ordered_sides[i][0] = i/3 + i/4 ;
	    ordered_sides[i][1] = (i + 1)%2 * (1 + i/2) + 4 * (i%2);
	    ordered_sides[i][2] = (i + 1)%2 * (3 + i/2 - i/4) + i * (i%2);
	    sides[i].insert(ordered_sides[i][0]);
	    sides[i].insert(ordered_sides[i][1]);
	    sides[i].insert(ordered_sides[i][2]);
	    if (debugging)
	    {
	        std::cout << name << std::endl;
		std::cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2];
		if (i == 0)
		    std::cout << " " << ordered_sides[i][3];
		std::cout <<  std::endl;
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
		side_types[i] = j;
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
		side_types[i] = j;
		ordered_sides[i][3] = i%4 + i/2 - 2 * (i/3);
	        sides[i].insert(ordered_sides[i][3]);
	    }
	    ordered_sides[i][0] = i/3 + 2 * (i/4) ;
	    ordered_sides[i][1] = (i + 1)%2 * (1 + i/2 + 2 * (i/4)) + 
	                          (i%2) * (3 + i/3);
	    ordered_sides[i][2] = (i + 1)%2 * (2 + 3 * (i/2) - 4 * (i/4)) +
	                          (i%2) * (4 + i/3);
	    sides[i].insert(ordered_sides[i][0]);
	    sides[i].insert(ordered_sides[i][1]);
	    sides[i].insert(ordered_sides[i][2]);
	    if (debugging)
	    {
	        std::cout << name << std::endl;
		std::cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2];
		if (i != 0 && i != nsides - 1)
		    std::cout << " " << ordered_sides[i][3];
		std::cout <<  std::endl;
	    }
	}

	else if (name == "hexahedron")
	{
	    ordered_sides[i][0] = i/3 + i/4 + 2 * (i/5);
	    ordered_sides[i][1] = ((i+1)%2) * (1 + i/2) + (i%2) * (4 + i/2);
	    ordered_sides[i][2] = ((i+1)%2) * (3 * (1 + i/2) - 2 * (i/4)) +
	                          (i%2) * (5 + 2 * (i/3));
	    ordered_sides[i][3] = ((i+1)%2) * (i + 2) + (i%2) * i;
	    sides[i].insert(ordered_sides[i][0]);
	    sides[i].insert(ordered_sides[i][1]);
	    sides[i].insert(ordered_sides[i][2]);
	    sides[i].insert(ordered_sides[i][3]);
	    if (debugging)
	    {
	        std::cout << name << std::endl;
		std::cout << i << " " << ordered_sides[i][0]  
		     << " " << ordered_sides[i][1]
		     << " " << ordered_sides[i][2]
		     << " " << ordered_sides[i][3] << std::endl;
	    }
	}
    }
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                    end of RTT_Format_Reader/CellDefs.cc
//---------------------------------------------------------------------------//
