//----------------------------------*-C++-*--------------------------------//
// Connect.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   meshReaders_Services/Connect.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for meshReaders_Services/Connect 
 *         class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Connect.hh"
#include "ds++/Assert.hh"
#include <utility>
#include <iterator>
#include <algorithm>
#include <iostream>

namespace rtt_meshReaders_Services
{
/*!
 * \brief Loads internal data structures from the mesh reader data.
 * \param mesh_reader_ meshReaders class object.
 */
void Connect::organizeData(rtt_dsxx::SP<Mesh_Reader> mesh_reader)
{
    nodes_coords = mesh_reader->get_node_coords();
    nnodes = nodes_coords.size();
    int ndim = nodes_coords[0].size();
    std::vector<Element_Type> element_types = mesh_reader->get_element_types();
    vector_vector_int element_nodes = mesh_reader->get_element_nodes();
    elem_sets = mesh_reader->get_element_sets();
    
    rtt_dsxx::SP<Element_Definition> element_definition;  
    for (int e = 0; e < element_types.size(); e++)
    {
        std::map<Element_Type, rtt_dsxx::SP<Element_Definition> >::iterator
	    elem_iter = elem_defs.find(element_types[e]);
        if (elem_iter == elem_defs.end())
	{
	    element_definition = new Element_Definition(element_types[e]);
	    elem_defs.insert(std::make_pair(element_types[e], 
			     element_definition));
	    elem_iter = elem_defs.find(element_types[e]);
	}
	if (elem_iter->second->get_dimension() == ndim)
	{
	    cells_types.push_back(element_types[e]);
	    cells_nodes.push_back(element_nodes[e]);
	}
	else
	{
	    sides_types.push_back(element_types[e]);
	    sides_nodes.push_back(element_nodes[e]);
	    for (int f = 0; f < bndry_flags.size(); f++)
	    {
	        std::map<std::string, std::set<int> >::iterator esiter = 
		    elem_sets.find(bndry_flags[f]);
		if (esiter->second.count(e) != 0)
		    sides_bndry.push_back(flag_numbs[f]);
	    }	        
	}
    }
    nsides = sides_nodes.size();
    ncells = cells_nodes.size();
    adjCell.resize(ncells);
}
/*!
 * \brief Generates a multimap containing the set of nodes (key) comprising 
 *        all of the sides.
 * \return The sideNodes multimap.
 */
std::multimap<std::set<int>, int> Connect::generateSideNodes()
{
    std::multimap<std::set<int>, int> sideNodes;
    for (int s = 0; s < nsides; ++s)
    {
        // Decide the side type to determine the number of nodes used in the 
        // side definition.
	int side_nnodes = 
	    elem_defs.find(sides_types[s])->second->get_number_of_nodes();

	set_int face;
	for (int node_index = 0; node_index < side_nnodes ; ++node_index)
	    face.insert(sides_nodes[s][node_index]);

	sideNodes.insert(std::make_pair(face, s));
    }
    return sideNodes;
}
/*!
 * \brief Generates a multimap containing the set of nodes (key) comprising 
 *        all of the cell faces.
 * \return The cellFaceNodes multimap.
 */
std::multimap<std::set<int>, int> Connect::generateCellFaceNodes()
{
    std::multimap<std::set<int>, int> cellFaceNodes;
    for (int c = 0; c < ncells; ++c)
    {
        // We have to decide the cell type to determine the number of
        // nodes used in the cell definition first for the general case.
	adjCell[c].resize(elem_defs.find(cells_types[c])->second->
			  get_number_of_sides());
	for (int s = 0; s < adjCell[c].size(); ++s)
	{
	    set_int face;
	    const vector_int & faceNodes = 
	        elem_defs.find(cells_types[c])->second->get_side_nodes(s);
	    for (int n = 0; n < faceNodes.size(); n++)
	         face.insert(cells_nodes[c][faceNodes[n]]);
	    cellFaceNodes.insert(std::make_pair(face, c));
	}
    }
    return cellFaceNodes;
}
/*!
 * \brief Generates data needed to treat the connectivity for complex cell 
 *        faces.
 * \param cell_faces The cell and face number.
 * \return The Connectivity class private data members cmplxNodes, cmplxCells,
 *         and circum are populated.
 */
void Connect::treatComplexFace(const vector_int & cell_faces)
{
    // this face is the boundary between complex cells. 
    Element_Definition faceDef = elem_defs.find(cells_types[cell_faces[0]])->
        second->get_side_type(cell_faces[1]);
    set_int lines;

    // create an ordered set of the cell face nodes to calculate
    // the face circumference length. Add all of these nodes to
    // the cmplxNodes multimap.
    const vector_int & cmplxFaceNodes = 
        elem_defs.find(cells_types[cell_faces[0]])->second->
        get_side_nodes(cell_faces[1]);
    vector_int cmplxGlobalNodes(cmplxFaceNodes.size());
    for (int n = 0; n < cmplxFaceNodes.size(); n++)
    {
       	cmplxGlobalNodes[n] = cells_nodes[cell_faces[0]][cmplxFaceNodes[n]];
       	cmplxNodes.insert(std::make_pair(cmplxGlobalNodes[n], cell_faces));
    }
    // create a multimap of the unconnected cells, faces, and 
    // nodes.
    cmplxCells.insert(std::make_pair(cell_faces, cmplxGlobalNodes));

    // build the face lines from the "side cell type" definition.
    // calculate the face circumference length and create a 
    // multimap of this information.
    double length = 0.;
    for (int s = 0; s < faceDef.get_number_of_nodes(); ++s)
    {
        vector_int lineVector = faceDef.get_side_nodes(s);
       	set_int lineNodes;
	for (int n = 0; n < lineVector.size(); n++)
	    lineNodes.insert(lineVector[n]);

        set_int::iterator bitr = lineNodes.begin();
    	set_int::reverse_iterator eitr = lineNodes.rbegin();
       	double line_length = 0;
        for (int d = 0; d < nodes_coords[0].size(); d++)
     	{
            double diff = nodes_coords[cmplxGlobalNodes[* bitr]][d] -
    		          nodes_coords[cmplxGlobalNodes[* eitr]][d];
	    line_length += std::pow(diff,2.0);
        }
        length += std::sqrt(line_length);
    }
    circum.insert(std::make_pair(length, cell_faces));
}
/*!
 * \brief Handles the connectivity for simple cell faces.
 * \param cell The cell number.
 * \param faceNum The cell face number.
 * \param mmiter Iterator to the cellFaceNodes multimap (modified).
 * \return The Connectivity class private data member adjCell is populated.
 */
void Connect::treatSimpleFace(const int & cell, const int & faceNum,
			      std::multimap<set_int, int>::iterator & mmiter)
{
    ++mmiter;
    const set_int & globalFaceNodes = mmiter->first;
    int otherCell = mmiter->second;
    int otherNnodes = elem_defs.find(cells_types[otherCell])->second->
        get_number_of_nodes();
    int otherNsides = elem_defs.find(cells_types[otherCell])->second->
        get_number_of_sides();
    set_int otherLocalFaceNodes;
    int otherFaceNum;
    for (int n = 0; n < otherNnodes; ++n)
    {
        // Determine if the next cell shares these nodes and store
        // the (normalized) data set in otherLocalFaceNodes.
     	int globalNode = cells_nodes[otherCell][n];
	if (globalFaceNodes.count(globalNode) == 1)
       	    otherLocalFaceNodes.insert(n);
    }
    // Correlate the otherlocalFaceNodes set with the user input face
    // numbering scheme and assign to otherfaceNum.
    for (int f = 0; f < otherNsides; ++f)
    {
	const vector_int & faceVector = elem_defs.find(cells_types[cell])->
	    second->get_side_nodes(f);
	set_int faceNodes;
	for (int n = 0; n < faceVector.size(); n++)
	    faceNodes.insert(faceVector[n]);
	if (otherLocalFaceNodes == faceNodes)
	    otherFaceNum = f;
    }
    adjCell[cell][faceNum].push_back(otherCell);
    adjCell[otherCell][otherFaceNum].push_back(cell);
}
/*!
 * \brief Determines the mesh connectivity for cells that share common faces.
 * \return The Connectivity class private data members cmplxNodes, cmplxCells,
 *         and circum are generated if complex faces (e.g., AMR) are found.
 */
void Connect::connectMutualFaces()
{
    vector_int cell_faces(2);
    // create a multimap of the node sets that define the sides so that 
    // the nodes are in ascending order.
    sideNodes = generateSideNodes();
    // create an iterator for the sides multimap.
    std::multimap<set_int, int>::iterator sitr = sideNodes.begin();

    // create a multimap of the node sets that define the cell faces so that 
    // the nodes are in ascending order.
    cellFaceNodes = generateCellFaceNodes();

    // loop over all of the cell faces. Iteration index is also incremented
    // once within the loop.
    for (std::multimap<set_int, int>::iterator mmiter = cellFaceNodes.begin();
	 mmiter != cellFaceNodes.end(); ++mmiter)
    {
	const set_int & globalFaceNodes = mmiter->first;
	int cell = mmiter->second;
	int otherCell;
	set_int localFaceNodes;
	int cell_nnodes = elem_defs.find(cells_types[cell])->second->
	    get_number_of_nodes();
	int cell_nsides = elem_defs.find(cells_types[cell])->second->
	    get_number_of_sides();
	int faceNum;

	// Determine which of this cell's nodes are tied to this face and put
	// this (normalized) info into the localFaceNodes set.
	for (int n = 0; n < cell_nnodes; ++n)
	{
	    int globalNode = cells_nodes[cell][n];
	    if (globalFaceNodes.count(globalNode) == 1)
		localFaceNodes.insert(n);
	}

	// Correlate the localFaceNodes set with the user input face numbering
	// scheme and assign the face number to faceNum.
	for (int f = 0; f < cell_nsides; ++f)
	{
	    const vector_int & faceVector = elem_defs.find(cells_types[cell])->
	        second->get_side_nodes(f);
	    set_int faceNodes;
	    for (int n = 0; n < faceVector.size(); n++)
	        faceNodes.insert(faceVector[n]);
	    if (localFaceNodes == faceNodes)
		faceNum = f;
	}
	// If the face only appears once in the definitions of the face cells
	// we must be on a boundary. For a simple mesh these faces are the
	// sides, while for an complex mesh these faces can also be next
	// to a hanging node. We can use this fact to resolve the refinement
	// level of a continuous adaptive mesh. Since the sides are sorted in 
	// the same order as the cell faces, we only have to look at the next 
	// element in the side multimap.
	if (cellFaceNodes.count(globalFaceNodes) == 1)
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
	    bndryFaces.insert(std::make_pair(faceNum,cell));
	    cell_faces[0] = cell;
	    cell_faces[1] = faceNum; 

	    const set_int & globalSideNodes = sitr->first;
	    int side = sitr->second;
	    // if the nodes that comprise this face also occur in the nodes
	    // that define the sides, assign the negative of the boundary side 
	    // flag as the adjacent cell. Correlate the cell face number to the
	    // corresponding side number in Cell_Faces_to_Sides.
	    if (globalSideNodes == globalFaceNodes)
	    {
	        adjCell[cell][faceNum].push_back(- sides_bndry[side]);
		++sitr;
		Side_to_Cell_Face.insert(std::make_pair(side, cell_faces));
	    }
	    // this face is the boundary between complex cells. 
	    else
	        treatComplexFace(cell_faces);
	}
	// this face is the junction to another simple cell (the next face 
	// in the multimap).
	else
	{
	    // temporary work-around to ICEM bug - 16 Aug 99.
	    const set_int & globalSideNodes = sitr->first;
	    int side = sitr->second;
	    if (globalSideNodes == globalFaceNodes)
	    {
	        ++mmiter;
		otherCell = mmiter->second;	        
		std::cout << "Warning: Face shared by cells " << cell + 1 
			  << " & " << otherCell + 1 
			  << " also exists as boundary side " << side + 1 
			  << std::endl;
		++sitr;
	        --mmiter;
	    }
	    treatSimpleFace(cell, faceNum, mmiter);
	}
    }
    // make sure all of the sides were correctly identified and assigned to 
    // the corresponding cell faces, then clear memory.
    Insist(sitr == sideNodes.end(),"Side/Cell Face correspondance not found!");
    sideNodes.clear();
    cellFaceNodes.clear();
}
/*!
 * \brief Prints the adjacent cell list.
 */
void Connect::printAdjacentCells() const
{
    for (int c = 0; c < ncells; ++c)
    {
        int cell_nsides = elem_defs.find(cells_types[c])->second->
	    get_number_of_sides();
        for (int f = 0; f < cell_nsides; f++)
        {
            std::cout << "cell " << c << " face " << f << " adjacent cell(s) ";
	    for (int af = 0; af < adjCell[c][f].size(); af++)
	        std::cout << get_adjCell(c,f,af) << " ";
	    std::cout << std::endl;
	}
    }
}
/*!
 * \brief Returns the cells that have boundary faces (i.e., faces that are 
 *        either on the outer boundary of the problem geometry or a connection
 *        between cells with different refinement levels in an AMR mesh) with 
 *        the specified face number
 * \param face Face number.
 * \return The cells with these boundary faces.
 */
std::set<int> Connect::get_bndryCells(int face) const
{ 
    set_int bndryCells;
    std::multimap<int, int>::const_iterator citer = bndryFaces.find(face);
    int cell_index = 0;
    int face_index = citer->first;

    while (citer !=  bndryFaces.end() && citer->first == face_index)
    {
        bndryCells.insert(citer->second);
        ++cell_index;
        ++citer;
    }
    return bndryCells;
}
/*!
 * \brief Returns true if the specified cell face is a boundary face (i.e., 
 *        a faces that is either on the outer boundary of the problem geometry
 *        or a connection between cells with different refinement levels in an
 *        AMR mesh).
 * \param cell Cell number.
 * \param face Face number.
 * \return Boundary face status.
 */
bool Connect::check_bndryFace(int cell, int face) const
{
    std::multimap<int, int>::const_iterator citer = bndryFaces.find(face);
    int face_index = citer->first;
    bool bndryFace = false;
	
    while (citer !=  bndryFaces.end() && citer->first == face_index &&
           !bndryFace)
    {
        if (citer->second == cell)
	    bndryFace = true;
	else
	    ++citer;
    }
    return bndryFace;
}
/*!
 * \brief Returns the cell number associated with the specified side number.
 * \param side Side number.
 * \return The cell number.
 */
int Connect::get_Cell_from_Side(int side) const
{
    std::multimap<int, vector_int >::const_iterator sitr = 
        Side_to_Cell_Face.find(side);

    if (sitr != Side_to_Cell_Face.end())
    {
        vector_int Cell_Face = sitr->second;
        return Cell_Face[0];
    }
    else
        return -1;
}
/*!
 * \brief Returns the cell face number associated with the specified side 
 *        number.
 * \param side The side number.
 * \return The face number.
 */
int Connect::get_Cell_Face_from_Side(int side) const
{
    std::multimap<int, vector_int >::const_iterator sitr = 
        Side_to_Cell_Face.find(side);

    if (sitr != Side_to_Cell_Face.end())
    {
        vector_int Cell_Face = sitr->second;
        return Cell_Face[1];
    }
    else
        return -1;
}

} // end namespace rtt_meshReaders_Services

//---------------------------------------------------------------------------//
//                 end of meshReaders_Services/Connect.cc
//---------------------------------------------------------------------------//
