//----------------------------------*-C++-*--------------------------------//
// RTT_Format_Connect.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/RTT_Format_Connect.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/RTT_Format_Connect class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "RTT_Format_Connect.hh"
#include "ds++/Assert.hh"
#include <utility>
#include <iterator>
#include <cmath>

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Constructs a Connectivity class object which determines the 
 *        mesh connectivity via a call to the calcAdjacentCells private
 *        memember function.
 * \param dims_ RTT_Format_Reader Dims class object.
 * \param cellDefs_ CellDefs class object.
 * \param cells_ Cells class object.
 * \param sides_ Sides class object.
 * \param nodes_ Nodes class object.
 */
Connectivity::Connectivity(const Dims & dims_, const CellDefs & cellDefs_, 
			   const Cells & cells_, const Sides & sides_,
			   const Nodes & nodes_) : dims(dims_), 
    cellDefs(cellDefs_), cells(cells_), sides(sides_), nodes(nodes_), 
    adjCell(dims.get_ncells(), vector_vector_int(dims.get_nsides_max(),
						 vector_int (1)))
{
    calcAdjacentCells();
}

/*!
 * \brief Determines the mesh connectivity using existing and accessible 
 *        Dims, CellDefs, Nodes, Sides, and Cells class objects.
 */
void Connectivity::calcAdjacentCells()
{
    std::multimap<set_int, int> faceCells;
    std::multimap<set_int, int> sideCells;
    std::multimap<vector_int, vector_int > unstCells;
    std::multimap<double, vector_int, std::greater<double> > circum;
    std::multimap<int, vector_int > unstNodes;
    vector_int cell_faces(2);

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
	    set_int face;
	    const set_int & faceNodes = cellDef.get_side(s);
	    for (set_int::const_iterator mmiter = faceNodes.begin();
		 mmiter != faceNodes.end(); ++mmiter)
	    {
	         int node_index = * mmiter;
	         face.insert(cells.get_nodes(c,node_index));
	    }
	    faceCells.insert(std::make_pair(face, c));
	}
    }

    // Do the same for the sides (but all of the side nodes are included in
    // the sets).
    for (int s = 0; s < dims.get_nsides(); ++s)
    {
        // Decide the side type to determine the number of  nodes used in the 
        // side definition.
	int nnodes = cellDefs.get_nnodes(sides.get_type(s));
	set_int face;

	for (int node_index = 0; node_index < nnodes ; ++node_index)
	    face.insert(sides.get_nodes(s,node_index));

	sideCells.insert(std::make_pair(face, s));
    }
    // find the side flag number that contains the boundary conditions
    int boundary = sides.get_boundary_flag_number();

    // create a separate iterator for the sides.
    std::multimap<set_int, int>::iterator sitr = sideCells.begin();

    // loop over all of the cell faces. Iteration index is also incremented
    // once within the loop.
    for (std::multimap<set_int, int>::iterator mmiter = faceCells.begin();
	 mmiter != faceCells.end(); ++mmiter)
    {
	const set_int & globalFaceNodes = mmiter->first;
	int cell = mmiter->second;
	int otherCell;
	set_int localFaceNodes;
	set_int otherLocalFaceNodes;
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
	    const set_int & faceNodes = cellDef.get_side(f);
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
	    bndryFaces.insert(std::make_pair(faceNum,cell));

	    const set_int & globalSideNodes = sitr->first;
	    int side = sitr->second;
	    // if the nodes that comprise this face also occur in the nodes
	    // that define the sides, assign the negative of the boundary side 
	    // flag as the adjacent cell. Correlate the cell face number to the
	    // corresponding side number in Cell_Faces_to_Sides.
	    if (globalSideNodes == globalFaceNodes)
	    {
		adjCell[cell][faceNum][0] = - sides.get_flags(side,boundary);
		++sitr;
		cell_faces[0] = cell;
		cell_faces[1] = faceNum;
		Side_to_Cell_Face.insert(std::make_pair(side, cell_faces));
	    }
	    // this face is the boundary between unstructured cells. 
	    else
	    {
		int faceType = cellDef.get_side_types(faceNum);
		const CellDef & faceDef = cellDefs.get_cell_def(faceType);
		vector_int unstGlobalNodes(faceDef.get_nnodes());
		set_int lines;
		cell_faces[0] = cell;
		cell_faces[1] = faceNum;

		// create an ordered set of the cell face nodes to calculate
		// the face circumference length. Add all of these nodes to
		// the unstNodes multimap.
		const vector_int & unstFaceNodes = 
		    cellDef.get_ordered_side(faceNum);
		for (int n = 0; n < faceDef.get_nnodes(); ++n)
		{
		    unstGlobalNodes[n] = 
		        cells.get_nodes(cell,unstFaceNodes[n]);
		    unstNodes.insert(std::make_pair(unstGlobalNodes[n],
						    cell_faces));
		}
		// create a multimap of the unconnected cells, faces, and 
		// nodes.
	        unstCells.insert(std::make_pair(cell_faces,unstGlobalNodes));

		// build the face lines from the "side cell type" definition.
		// calculate the face circumference length and create a 
		// multimap of this information.
		double length = 0.;
		for (int s = 0; s < faceDef.get_nsides(); ++s)
		{
		    set_int lineNodes = faceDef.get_side(s);

		    set_int::iterator bitr = lineNodes.begin();
		    set_int::reverse_iterator eitr = lineNodes.rbegin();
		    double line_length = 0;
		    for (int d = 0; d < dims.get_ndim(); d++)
		    {
		        double diff = 
			    nodes.get_coords(unstGlobalNodes[* bitr],d) -
			    nodes.get_coords(unstGlobalNodes[* eitr],d);
		            line_length += std::pow(diff,2.0);
		    }
		    length += std::sqrt(line_length);
		}
		circum.insert(std::make_pair(length,cell_faces));
	    }
	}
	// this face is the junction to another structured cell (the next 
	// face in the multimap).
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
	    ++mmiter;
	    const set_int & globalFaceNodes = mmiter->first;
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
		const set_int & faceNodes = cellDef.get_side(f);
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

    std::multimap<double, vector_int, std::greater<double> >::iterator cfiter= 
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
	vector_int & faceNodes = unstCells.find(cell_faces)->second;
	vector_dbl center(dims.get_ndim());
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
	std::multimap<int,vector_int >::iterator fiter = 
	    unstNodes.find(centerNode);

	// remove the null value from the end of adjCell vector.
	if (!adjCell[cell][faceNum].empty())
	    adjCell[cell][faceNum].pop_back();

	// assign all of the cell faces that contain this node to this whole
	// face, and assign this whole face to each of the subfaces.
	while (fiter !=  unstNodes.end() && fiter->first == centerNode)
	{
	    vector_int other_cell_faces = fiter->second;
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
	        std::cout << "cell " << c << " face " << f 
			  << " adjacent cell(s) ";
	        for (int af = 0; af < adjCell[c][f].size(); af++)
	            std::cout << get_adjCell(c,f,af) << " ";
		std::cout << std::endl;
	    }
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
std::set<int> Connectivity::get_bndryCells(int face) const
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
bool Connectivity::check_bndryFace(int cell, int face) const
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
int Connectivity::get_Cell_from_Side(int side) const
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
 * \param side Side number.
 * \return The face number.
 */
int Connectivity::get_Cell_Face_from_Side(int side) const
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

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                        end of RTT_Format_Connect.cc
//---------------------------------------------------------------------------//
